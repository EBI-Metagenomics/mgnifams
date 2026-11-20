#!/usr/bin/env python3
"""Self-check: select_nonredundant_families.py. Run: python3 bin/test_select_nonredundant_families.py"""

import tempfile
from pathlib import Path

from select_nonredundant_families import main

EDGES_HEADER = "family_1,family_2,shared_proteins,protein_jaccard,residue_jaccard,residue_containment\n"
METADATA_HEADER = (
    "family_id,converged,seed_msa_size,full_msa_size,"
    "rep_protein,rep_region,rep_length,consensus_length,rep_sequence,consensus_sequence\n"
)


def run(edges, sizes, score="0.95", residue="0.95"):
    """`sizes`: {chunk: {local id: full_msa_size}}; returns the redundant names."""
    with tempfile.TemporaryDirectory() as tmp:
        d = Path(tmp)
        (d / "edges.csv").write_text(EDGES_HEADER + "".join(f"{a},{b},1,{p},{r},1.000000\n" for a, b, p, r in edges))
        metadata = []
        for chunk, families in sizes.items():
            path = d / f"{chunk}_metadata.csv"
            rows = "".join(f'{i},True,1,{n},"1",1-10,10,10,AAA,AAA\n' for i, n in families.items())
            path.write_text(METADATA_HEADER + rows)
            metadata.append(str(path))
        main(
            ["--edges", str(d / "edges.csv"), "--metadata", *metadata, "--score_threshold", score]
            + ["--residue_threshold", residue, "--out", str(d / "out.txt")]
        )
        return (d / "out.txt").read_text().split()


# Chain A > B > C with A~B and B~C: A and C kept
assert run([("c_1", "c_2", "1.0", "1.0"), ("c_2", "c_3", "1.0", "1.0")], {"c": {1: 30, 2: 20, 3: 10}}) == ["c_2"]
# Equal sizes: the name that sorts first is kept
assert run([("c_2", "c_1", "0.96", "0.96")], {"c": {1: 10, 2: 10}}) == ["c_2"]
# A pair failing either threshold is not redundant; values at the threshold are
assert run([("c_1", "c_2", "0.949999", "1.0")], {"c": {1: 10, 2: 5}}) == []
assert run([("c_1", "c_2", "1.0", "0.949999")], {"c": {1: 10, 2: 5}}) == []
assert run([("c_1", "c_2", "0.950000", "0.950000")], {"c": {1: 10, 2: 5}}) == ["c_2"]
# Containment-only edge (residue Jaccard below threshold): not redundant
assert run([("c_1", "c_2", "1.0", "0.300000")], {"c": {1: 100, 2: 30}}) == []
# Empty edge file
assert run([], {"c": {1: 10}}) == []
# Two chunks with local family 1 stay two families
assert run([("000_1", "001_1", "1.0", "1.0")], {"000": {1: 5}, "001": {1: 9}}) == ["000_1"]
# Thresholds are user-set
assert run([("c_1", "c_2", "0.5", "0.5")], {"c": {1: 10, 2: 5}}, score="0.5", residue="0.5") == ["c_2"]
# A family in the edges but in no metadata fails loudly
try:
    run([("c_1", "x_9", "1.0", "1.0")], {"c": {1: 10}})
    raise AssertionError("expected ValueError")
except ValueError as error:
    assert "x_9" in str(error)

print("test_select_nonredundant_families.py: all passed")
