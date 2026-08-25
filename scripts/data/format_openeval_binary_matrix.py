#!/usr/bin/env python3
"""Format OpenEval benchmark responses as an LLM-by-item binary matrix.

The OpenEval dataset stores benchmark-level item files and sharded response
files on Hugging Face:

    item/<benchmark>-00000-of-00001.parquet
    response/<benchmark>-00000-of-0000K.parquet

Each response has a nested model descriptor and a nested score object.  This
script extracts model names, item identifiers, and numeric scores, then writes a
binary response matrix suitable for the probit factor-model scripts.

Examples:

    python3 format_openeval_binary_matrix.py \
      --benchmarks ifeval \
      --out-dir openeval_ifeval_formatted_uncapped \
      --min-item-response-prop 0.25 \
      --min-model-response-prop 0.25

The command above keeps the full retained IFEval item set after coverage and
constant-item filtering.  Do not use an item cap for the IFEval analysis in
this repository.  The optional cap is only for separate exploratory
multi-benchmark skill-battery subsets, e.g.

    python3 format_openeval_binary_matrix.py \
      --benchmarks ifeval,gpqa,mmlu_pro \
      --out-dir openeval_llm_skill_battery_formatted \
      --min-item-response-prop 0.25 \
      --min-model-response-prop 0.25 \
      --max-items-per-benchmark 500

Matrix entries are not read from a preexisting strict-accuracy column.  The
script extracts numeric values from the nested OpenEval score object, averages
multiple numeric values for a model-item pair, and codes the entry as 1 when
the mean score is at least --binary-threshold, which defaults to 0.5.
"""

from __future__ import annotations

import argparse
import json
import os
import re
from collections import defaultdict
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd

try:
    from huggingface_hub import hf_hub_download, list_repo_files
except ImportError:  # local cached snapshots can be used without this package
    hf_hub_download = None
    list_repo_files = None


REPO_ID = "human-centered-eval/OpenEval"
REPO_TYPE = "dataset"
ITEM_ID_RE = re.compile(r"^(?P<benchmark>.+?)_(?P<timestamp>\d{8}T\d{6}Z)_(?P<item_index>\d+)")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Build binary model-by-item matrices from OpenEval responses."
    )
    parser.add_argument(
        "--benchmarks",
        default="ifeval,gpqa,mmlu_pro",
        help="Comma-separated benchmark names. Use --list-benchmarks to inspect options.",
    )
    parser.add_argument(
        "--out-dir",
        default="openeval_llm_skill_battery_formatted",
        help="Directory for CSV/JSON outputs.",
    )
    parser.add_argument(
        "--cache-dir",
        default=None,
        help="Optional local cache/download directory.",
    )
    parser.add_argument(
        "--local-snapshot-dir",
        default=None,
        help=(
            "Optional path to a local OpenEval snapshot containing item/ and "
            "response/ directories. When supplied, no Hugging Face network "
            "calls are needed."
        ),
    )
    parser.add_argument(
        "--token",
        default=None,
        help="Optional Hugging Face token. If omitted, HF_TOKEN/HUGGINGFACE_HUB_TOKEN/cache is used.",
    )
    parser.add_argument(
        "--list-benchmarks",
        action="store_true",
        help="Print available response benchmark names and exit.",
    )
    parser.add_argument(
        "--search",
        default=None,
        help="When listing benchmarks, only print names containing this substring.",
    )
    parser.add_argument(
        "--binary-threshold",
        type=float,
        default=0.5,
        help="Mean numeric score >= this value is coded as correct.",
    )
    parser.add_argument(
        "--min-item-response-prop",
        type=float,
        default=0.50,
        help="Drop items answered by less than this fraction of retained models.",
    )
    parser.add_argument(
        "--min-model-response-prop",
        type=float,
        default=0.50,
        help="Drop models answering less than this fraction of retained items.",
    )
    parser.add_argument(
        "--drop-constant-items",
        action=argparse.BooleanOptionalAction,
        default=True,
        help="Drop items that are all correct or all incorrect among observed responses.",
    )
    parser.add_argument(
        "--max-items-per-benchmark",
        type=int,
        default=None,
        help="Optional reproducible cap on retained items per benchmark.",
    )
    parser.add_argument(
        "--seed",
        type=int,
        default=20260727,
        help="Seed used only for item subsampling.",
    )
    return parser.parse_args()


def token_from_args(args: argparse.Namespace) -> str | None:
    return (
        args.token
        or os.environ.get("HF_TOKEN")
        or os.environ.get("HUGGINGFACE_HUB_TOKEN")
        or None
    )


def repo_files(token: str | None) -> list[str]:
    if list_repo_files is None:
        raise ImportError(
            "huggingface_hub is required unless --local-snapshot-dir is supplied."
        )
    return list_repo_files(REPO_ID, repo_type=REPO_TYPE, token=token)


def local_repo_files(snapshot_dir: str) -> list[str]:
    root = Path(snapshot_dir)
    return sorted(str(path.relative_to(root)) for path in root.rglob("*") if path.is_file())


def available_response_benchmarks(files: list[str]) -> list[str]:
    benches = set()
    for path in files:
        match = re.match(r"^response/(.+)-\d{5}-of-\d{5}\.parquet$", path)
        if match:
            benches.add(match.group(1))
    return sorted(benches)


def response_shards(files: list[str], benchmark: str) -> list[str]:
    prefix = f"response/{benchmark}-"
    return sorted(
        path
        for path in files
        if path.startswith(prefix) and path.endswith(".parquet")
    )


def download_file(path: str, args: argparse.Namespace, token: str | None) -> Path:
    if args.local_snapshot_dir is not None:
        local_path = Path(args.local_snapshot_dir) / path
        if not local_path.exists():
            raise FileNotFoundError(f"Local OpenEval snapshot file not found: {local_path}")
        return local_path
    if hf_hub_download is None:
        raise ImportError(
            "huggingface_hub is required unless --local-snapshot-dir is supplied."
        )
    kwargs: dict[str, Any] = {
        "repo_id": REPO_ID,
        "filename": path,
        "repo_type": REPO_TYPE,
        "token": token,
    }
    if args.cache_dir is not None:
        kwargs["local_dir"] = args.cache_dir
    return Path(hf_hub_download(**kwargs))


def to_plain(value: Any) -> Any:
    """Convert nested numpy/pandas containers to JSON-friendly Python values."""
    if value is None or value is pd.NA:
        return None
    if isinstance(value, np.generic):
        return value.item()
    if isinstance(value, np.ndarray):
        return [to_plain(v) for v in value.tolist()]
    if isinstance(value, dict):
        return {str(k): to_plain(v) for k, v in value.items()}
    if isinstance(value, (list, tuple)):
        return [to_plain(v) for v in value]
    return value


def as_json_text(value: Any, max_chars: int = 8000) -> str:
    try:
        text = json.dumps(to_plain(value), ensure_ascii=True, sort_keys=True)
    except TypeError:
        text = str(value)
    return text[:max_chars]


def parse_json_maybe(value: Any) -> Any:
    """Parse JSON strings inside OpenEval payloads when possible."""
    value = to_plain(value)
    if isinstance(value, str):
        try:
            return json.loads(value)
        except json.JSONDecodeError:
            return value
    return value


def first_nonempty(values: list[Any]) -> str:
    for value in values:
        if value is None:
            continue
        text = str(value)
        if text:
            return text
    return ""


def extract_ifeval_item_details(item_content: Any) -> dict[str, Any]:
    """Extract prompt and instruction metadata from an OpenEval IFEval item.

    OpenEval stores IFEval prompts in item_content["input"] as JSON strings.
    Each JSON object contains the original prompt, the instruction_id_list used
    by the strict evaluator, and the instruction-specific kwargs.  Keeping
    these columns explicit makes it possible to join item loadings back to the
    raw instruction checks.
    """
    content = parse_json_maybe(item_content)
    if not isinstance(content, dict):
        return {
            "prompt": "",
            "instruction_ids": "",
            "instruction_families": "",
            "n_instructions": 0,
            "n_unique_instruction_ids": 0,
            "instruction_kwargs": "[]",
        }

    raw_inputs = content.get("input", [])
    if not isinstance(raw_inputs, list):
        raw_inputs = [raw_inputs]

    prompts: list[str] = []
    instruction_ids: list[str] = []
    kwargs_all: list[Any] = []
    for raw_input in raw_inputs:
        parsed = parse_json_maybe(raw_input)
        if not isinstance(parsed, dict):
            continue
        prompts.append(str(parsed.get("prompt", "")))
        ids = parsed.get("instruction_id_list", [])
        if not isinstance(ids, list):
            ids = [ids]
        instruction_ids.extend(str(x) for x in ids if x is not None)
        kwargs = parsed.get("kwargs", [])
        if not isinstance(kwargs, list):
            kwargs = [kwargs]
        kwargs_all.extend(kwargs)

    families = [x.split(":", 1)[0] for x in instruction_ids]
    return {
        "prompt": first_nonempty(prompts),
        "instruction_ids": "|".join(instruction_ids),
        "instruction_families": "|".join(families),
        "n_instructions": len(instruction_ids),
        "n_unique_instruction_ids": len(set(instruction_ids)),
        "instruction_kwargs": as_json_text(kwargs_all, max_chars=12000),
    }


def item_instruction_long(item_meta: pd.DataFrame, retained_items: set[str]) -> pd.DataFrame:
    """Build one row per retained item-instruction pair."""
    if item_meta.empty or "item_id" not in item_meta.columns:
        return pd.DataFrame()
    rows: list[dict[str, Any]] = []
    meta = item_meta[item_meta["item_id"].isin(retained_items)].copy()
    for _, row in meta.iterrows():
        ids = [x for x in str(row.get("instruction_ids", "")).split("|") if x]
        families = [x for x in str(row.get("instruction_families", "")).split("|") if x]
        for k, instruction_id in enumerate(ids, start=1):
            family = families[k - 1] if k - 1 < len(families) else instruction_id.split(":", 1)[0]
            rows.append(
                {
                    "benchmark": row.get("benchmark", ""),
                    "item_id": row.get("item_id", ""),
                    "instruction_position": k,
                    "instruction_id": instruction_id,
                    "instruction_family": family,
                    "prompt": row.get("prompt", ""),
                }
            )
    return pd.DataFrame(rows)


def extract_model_name(model_obj: Any) -> str | None:
    model_obj = to_plain(model_obj)
    if isinstance(model_obj, dict):
        return model_obj.get("name") or model_obj.get("model_name") or model_obj.get("id")
    if model_obj is None:
        return None
    return str(model_obj)


def numeric_values(value: Any) -> list[float]:
    """Recursively pull finite numeric values from a nested score object."""
    value = to_plain(value)
    if value is None:
        return []
    if isinstance(value, bool):
        return [float(value)]
    if isinstance(value, (int, float)):
        val = float(value)
        return [val] if np.isfinite(val) else []
    if isinstance(value, dict):
        if "value" in value:
            vals = numeric_values(value["value"])
            if vals:
                return vals
        vals: list[float] = []
        for nested in value.values():
            vals.extend(numeric_values(nested))
        return vals
    if isinstance(value, (list, tuple)):
        vals = []
        for nested in value:
            vals.extend(numeric_values(nested))
        return vals
    return []


def extract_score(score_obj: Any) -> float:
    vals = numeric_values(score_obj)
    if not vals:
        return np.nan
    return float(np.mean(vals))


def extract_metric_names(score_obj: Any) -> str:
    score_obj = to_plain(score_obj)
    metrics = []
    if isinstance(score_obj, dict) and "metric" in score_obj:
        metric_obj = score_obj["metric"]
        if isinstance(metric_obj, list):
            for metric in metric_obj:
                if isinstance(metric, dict) and metric.get("name"):
                    metrics.append(str(metric["name"]))
        elif isinstance(metric_obj, dict) and metric_obj.get("name"):
            metrics.append(str(metric_obj["name"]))
    return "|".join(sorted(set(metrics)))


def item_id_from_response_id(response_id: Any) -> str | None:
    if response_id is None:
        return None
    response_id = str(response_id)
    match = ITEM_ID_RE.match(response_id)
    if not match:
        return None
    return f"{match.group('benchmark')}_{match.group('timestamp')}_{match.group('item_index')}"


def read_item_metadata(benchmark: str, args: argparse.Namespace, token: str | None) -> pd.DataFrame:
    path = f"item/{benchmark}-00000-of-00001.parquet"
    try:
        item_path = download_file(path, args, token)
    except Exception:
        return pd.DataFrame({"benchmark": [], "item_id": []})

    item_df = pd.read_parquet(item_path).copy()
    if "item_id" not in item_df.columns:
        item_df["item_id"] = [f"{benchmark}_{i}" for i in range(len(item_df))]
    if "item_content" in item_df.columns:
        details = pd.DataFrame(
            item_df["item_content"].map(extract_ifeval_item_details).tolist(),
            index=item_df.index,
        )
        item_df = pd.concat([item_df, details], axis=1)
    item_df.insert(0, "benchmark", benchmark)

    for col in item_df.columns:
        if col not in {
            "benchmark",
            "item_id",
            "prompt",
            "instruction_ids",
            "instruction_families",
            "n_instructions",
            "n_unique_instruction_ids",
        }:
            item_df[col] = item_df[col].map(as_json_text)
    return item_df


def read_response_long(
    benchmark: str,
    shard_paths: list[str],
    args: argparse.Namespace,
    token: str | None,
) -> pd.DataFrame:
    parts = []
    for shard in shard_paths:
        parquet_path = download_file(shard, args, token)
        df = pd.read_parquet(parquet_path, columns=["response_id", "model", "scores"])
        part = pd.DataFrame(
            {
                "benchmark": benchmark,
                "response_id": df["response_id"].astype(str),
                "item_id": df["response_id"].map(item_id_from_response_id),
                "model_name": df["model"].map(extract_model_name),
                "score": df["scores"].map(extract_score),
                "metric_names": df["scores"].map(extract_metric_names),
            }
        )
        parts.append(part)

    out = pd.concat(parts, ignore_index=True)
    out = out.dropna(subset=["item_id", "model_name", "score"])
    out["score"] = pd.to_numeric(out["score"], errors="coerce")
    out = out.dropna(subset=["score"])
    return out


def filter_matrix(
    matrix: pd.DataFrame,
    args: argparse.Namespace,
) -> pd.DataFrame:
    """Iteratively apply item/model coverage filters and optional item cap."""
    mat = matrix.copy()

    for _ in range(10):
        old_shape = mat.shape
        if mat.shape[0] > 0 and mat.shape[1] > 0:
            item_keep = mat.notna().mean(axis=0) >= args.min_item_response_prop
            mat = mat.loc[:, item_keep]
        if mat.shape[0] > 0 and mat.shape[1] > 0:
            model_keep = mat.notna().mean(axis=1) >= args.min_model_response_prop
            mat = mat.loc[model_keep, :]
        if mat.shape == old_shape:
            break

    if args.drop_constant_items and mat.shape[1] > 0:
        observed_min = mat.min(axis=0, skipna=True)
        observed_max = mat.max(axis=0, skipna=True)
        mat = mat.loc[:, observed_min != observed_max]

    if args.max_items_per_benchmark is not None and mat.shape[1] > 0:
        rng = np.random.default_rng(args.seed)
        by_bench: dict[str, list[str]] = defaultdict(list)
        for item_id in mat.columns:
            bench = item_id.split("_", 1)[0]
            if item_id.startswith("mmlu_pro_"):
                bench = "mmlu_pro"
            elif item_id.startswith("omni_math_"):
                bench = "omni_math"
            elif item_id.startswith("hi_tom_"):
                bench = "hi_tom"
            elif item_id.startswith("do_not_answer_"):
                bench = "do_not_answer"
            by_bench[bench].append(item_id)

        keep = []
        for items in by_bench.values():
            items = sorted(items)
            if len(items) > args.max_items_per_benchmark:
                items = sorted(rng.choice(items, size=args.max_items_per_benchmark, replace=False))
            keep.extend(items)
        mat = mat.loc[:, sorted(keep)]

    return mat


def write_outputs(
    response_long: pd.DataFrame,
    retained_pairs: pd.DataFrame,
    item_meta: pd.DataFrame,
    matrix_raw: pd.DataFrame,
    matrix_with_missing: pd.DataFrame,
    args: argparse.Namespace,
) -> None:
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    complete = matrix_with_missing.dropna(axis=0, how="any").dropna(axis=1, how="any")
    if complete.shape[1] > 0 and args.drop_constant_items:
        complete = complete.loc[:, complete.min(axis=0) != complete.max(axis=0)]

    response_long.to_csv(out_dir / "openeval_response_long_scores.csv", index=False)
    retained_pairs.to_csv(out_dir / "openeval_response_long_binary.csv", index=False)
    matrix_raw.to_csv(out_dir / "openeval_binary_matrix_raw.csv")
    matrix_with_missing.to_csv(out_dir / "openeval_binary_matrix_with_missing.csv")
    complete.to_csv(out_dir / "openeval_binary_matrix_complete.csv")
    complete_nonconstant = complete.loc[:, complete.min(axis=0) != complete.max(axis=0)] if complete.shape[1] else complete
    complete_nonconstant.to_csv(out_dir / "openeval_binary_matrix_complete_nonconstant.csv")

    retained_items = set(matrix_with_missing.columns)
    item_meta = item_meta[item_meta["item_id"].isin(retained_items)].copy()
    item_meta.to_csv(out_dir / "openeval_item_metadata.csv", index=False)
    instructions = item_instruction_long(item_meta, retained_items)
    instructions.to_csv(out_dir / "openeval_item_instruction_metadata_long.csv", index=False)

    model_meta = pd.DataFrame({"model_name": matrix_with_missing.index})
    model_meta["n_observed_items"] = matrix_with_missing.notna().sum(axis=1).to_numpy()
    model_meta["mean_correct_observed"] = matrix_with_missing.mean(axis=1, skipna=True).to_numpy()
    model_meta.to_csv(out_dir / "openeval_model_metadata.csv", index=False)

    item_summary = pd.DataFrame(
        {
            "item_id": matrix_with_missing.columns,
            "n_observed_models": matrix_with_missing.notna().sum(axis=0).to_numpy(),
            "mean_correct_observed": matrix_with_missing.mean(axis=0, skipna=True).to_numpy(),
        }
    )
    item_summary["benchmark"] = item_summary["item_id"].map(lambda x: x.rsplit("_", 2)[0])
    item_summary.to_csv(out_dir / "openeval_item_summary.csv", index=False)

    benchmark_summary = (
        response_long.groupby("benchmark")
        .agg(
            raw_rows=("score", "size"),
            raw_models=("model_name", "nunique"),
            raw_items=("item_id", "nunique"),
            mean_score=("score", "mean"),
        )
        .reset_index()
    )
    benchmark_summary.to_csv(out_dir / "openeval_benchmark_summary.csv", index=False)

    summary = {
        "repo_id": REPO_ID,
        "benchmarks": sorted(response_long["benchmark"].unique().tolist()),
        "binary_threshold": args.binary_threshold,
        "min_item_response_prop": args.min_item_response_prop,
        "min_model_response_prop": args.min_model_response_prop,
        "drop_constant_items": args.drop_constant_items,
        "max_items_per_benchmark": args.max_items_per_benchmark,
        "raw_long_rows": int(len(response_long)),
        "retained_unique_instruction_ids": (
            int(instructions["instruction_id"].nunique()) if not instructions.empty else 0
        ),
        "retained_unique_instruction_families": (
            int(instructions["instruction_family"].nunique()) if not instructions.empty else 0
        ),
        "raw_matrix_shape": list(matrix_raw.shape),
        "matrix_with_missing_shape": list(matrix_with_missing.shape),
        "matrix_with_missing_mean_correct": float(np.nanmean(matrix_with_missing.to_numpy())),
        "complete_matrix_shape": list(complete.shape),
        "complete_nonconstant_shape": list(complete_nonconstant.shape),
        "complete_nonconstant_mean_correct": (
            float(np.nanmean(complete_nonconstant.to_numpy()))
            if complete_nonconstant.size
            else None
        ),
        "output_dir": str(out_dir.resolve()),
    }
    with open(out_dir / "openeval_format_summary.json", "w", encoding="utf-8") as handle:
        json.dump(summary, handle, indent=2, sort_keys=True)
    pd.DataFrame([summary]).to_csv(out_dir / "openeval_format_summary.csv", index=False)

    print(json.dumps(summary, indent=2, sort_keys=True))


def main() -> None:
    args = parse_args()
    token = token_from_args(args)
    files = local_repo_files(args.local_snapshot_dir) if args.local_snapshot_dir else repo_files(token)

    if args.list_benchmarks:
        benches = available_response_benchmarks(files)
        if args.search:
            benches = [b for b in benches if args.search.lower() in b.lower()]
        print("\n".join(benches))
        return

    requested = [b.strip() for b in args.benchmarks.split(",") if b.strip()]
    all_benchmarks = set(available_response_benchmarks(files))
    missing = sorted(set(requested) - all_benchmarks)
    if missing:
        raise ValueError(f"Benchmarks not found in response files: {missing}")

    response_parts = []
    item_parts = []
    for benchmark in requested:
        shards = response_shards(files, benchmark)
        if not shards:
            raise ValueError(f"No response shards found for benchmark: {benchmark}")
        print(f"Reading {benchmark}: {len(shards)} response shard(s)")
        response_parts.append(read_response_long(benchmark, shards, args, token))
        item_parts.append(read_item_metadata(benchmark, args, token))

    response_long = pd.concat(response_parts, ignore_index=True)
    item_meta = pd.concat(item_parts, ignore_index=True) if item_parts else pd.DataFrame()

    collapsed = (
        response_long.groupby(["model_name", "item_id"], as_index=False)
        .agg(
            benchmark=("benchmark", "first"),
            score=("score", "mean"),
            n_attempts=("score", "size"),
            metric_names=("metric_names", lambda x: "|".join(sorted(set("|".join(x).split("|"))))),
        )
    )
    collapsed["correct"] = (collapsed["score"] >= args.binary_threshold).astype(int)
    matrix_raw = collapsed.pivot(index="model_name", columns="item_id", values="correct").sort_index()
    matrix = filter_matrix(matrix_raw, args)

    retained_pairs = collapsed[
        collapsed["model_name"].isin(matrix.index) & collapsed["item_id"].isin(matrix.columns)
    ].copy()
    write_outputs(response_long, retained_pairs, item_meta, matrix_raw, matrix, args)


if __name__ == "__main__":
    main()
