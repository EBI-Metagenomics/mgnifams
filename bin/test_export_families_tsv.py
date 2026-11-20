#!/usr/bin/env python3
"""Self-check for export_families_tsv.py. Run: python3 bin/test_export_families_tsv.py"""

import os
import tempfile

from export_families_tsv import COLUMNS, new_rows, read_delta, update_rows

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

metadata = {"1": {"full_msa_size": "47", "rep_protein": "1814953751", "rep_region": "541-660"}}
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
# New families: one `new` row each, every release column from this release, whole protein as `-`
new_metadata = {
    "35460": {"full_msa_size": "31", "seed_msa_size": "22", "rep_protein": "782510898", "rep_region": "278-392"},
    "35461": {"full_msa_size": "25", "seed_msa_size": "20", "rep_protein": "376297851", "rep_region": "-"},
}
new = new_rows(new_metadata, "2.0", "2026_07")
assert [list(r) for r in new] == [COLUMNS, COLUMNS]
assert new[0] == {
    "family_id": "35460",
    "status": "new",
    "not_updated_reason": "",
    "first_release": "2.0",
    "model_release": "2.0",
    "model_proteins_release": "2026_07",
    "members_release": "2.0",
    "members_proteins_release": "2026_07",
    "full_size": "31",
    "seed_size": "22",
    "rep_id": "MGYP000782510898",
    "rep_region": "278-392",
}, new[0]
assert new[1]["rep_region"] == "-" and new[1]["rep_id"] == "MGYP000376297851"
assert list(previous[0]) == COLUMNS  # the update path keeps the same columns
print("ok")
