# Tests

This directory is reserved for the audit pass.

Suggested first tests:

- truncated-normal augmentation respects binary response signs
- intercept centering leaves `alpha + Lambda f` invariant after factor-location normalization
- permutation/sign alignment recovers a known synthetic loading matrix
- flattened parameter correlation is unchanged by a signed permutation after alignment
- Gibbs joint-profile indexing has exactly `G^H` profiles
