"""Small-system brute-force enumeration helpers for cross-check tests."""

from __future__ import annotations

from itertools import combinations
from math import comb


def analytic_binary_labelings(n_sites: int, n_a: int) -> int:
    """Return the raw number of A/B labelings before symmetry reduction."""
    if not 0 <= n_a <= n_sites:
        raise ValueError("n_a must be between 0 and n_sites")
    return comb(n_sites, n_a)


def brute_force_binary_patterns(n_sites: int, n_a: int) -> list[tuple[int, ...]]:
    """Generate raw binary occupation patterns for small tests.

    A site value of 1 means species A; 0 means species B/vacancy.
    """
    if n_sites > 24:
        raise ValueError("brute-force helper is limited to <=24 sites")
    patterns = []
    for occ in combinations(range(n_sites), n_a):
        row = [0] * n_sites
        for i in occ:
            row[i] = 1
        patterns.append(tuple(row))
    return patterns
