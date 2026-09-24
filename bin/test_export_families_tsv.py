#!/usr/bin/env python3
"""Self-check for export_families_tsv.py. Run: python3 bin/test_export_families_tsv.py"""

import os
import tempfile

from export_families_tsv import read_delta, update_rows

row = {
    "family_id": "",
    "status": "new",
    "not_updated_reason": "",
    "first_release": "1.0",
    "model_release": "1.0",
    "model_proteins_release": "2024_04",
    "members_release": "1.0",
    "members_proteins_release": "2024_04",
    "full_size": "10",
    "seed_size": "5",
    "rep_id": "MGYP000000000001",
    "rep_region": "1-100",
}
previous = [row | {"family_id": "1"}, row | {"family_id": "2"}]

with tempfile.TemporaryDirectory() as tmp:
    delta = os.path.join(tmp, "updated_delta.csv")
    with open(delta, "w") as handle:
        handle.write(
            "family_id,model_length_before,model_length_after,round1_recruits,full_msa_size,retention,rounds_run,converged,outcome\n"
            "1,122,122,43,47,1.0,1,False,successful\n"
            "2,114,,,,,1,False,internal error during a, b\n"
        )
    outcomes = read_delta(delta)
assert outcomes == {"1": "successful", "2": "internal error during a, b"}, outcomes

metadata = {"1": {"full_msa_size": "47", "protein": "1814953751", "region": "541-660"}}
updated, not_updated = update_rows(previous, outcomes, metadata, "1.1", "2026_07")
assert updated == row | {
    "family_id": "1",
    "status": "updated",
    "members_release": "1.1",
    "members_proteins_release": "2026_07",
    "full_size": "47",
    "rep_id": "MGYP001814953751",
    "rep_region": "541-660",
}, updated
assert not_updated == previous[1] | {"status": "not_updated", "not_updated_reason": "internal error during a, b"}


def fails(*args):
    try:
        update_rows(*args)
    except SystemExit:
        return True
    return False


assert fails(previous, {"1": "successful"}, metadata, "1.1", "2026_07")  # previous family missing from the delta
assert fails(previous, outcomes | {"3": "successful"}, metadata, "1.1", "2026_07")  # unknown family
assert fails(previous, outcomes, {}, "1.1", "2026_07")  # successful family without metadata
assert fails(previous, outcomes, metadata, "1.0", "2026_07")  # previous file already from this release
print("ok")
