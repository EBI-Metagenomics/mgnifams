#!/usr/bin/env python3
"""Self-check: pool_nonredundant_families.py on real mgnifam 4.1.0 generate_families outputs of two chunks that both
hold local families 1 and 2 (modules/local/pool_nonredundant_families/tests/fixtures/mgnifam_4.1.0.tar.gz, made once with
`mgnifam generate_families --chunk_id 000|001` on two halves of the nf-core test-datasets
proteomics/proteinfamilies/mgnify_proteins_trimmed_clustering.tsv, the .ssi, log and stats removed).
Run: python3 bin/test_pool_nonredundant_families.py"""

import gzip
import json
import shutil
import tarfile
import tempfile
from pathlib import Path

from pool_nonredundant_families import METADATA_HEADER, main

FIXTURE_TAR = (
    Path(__file__).resolve().parent.parent
    / "modules/local/pool_nonredundant_families/tests/fixtures/mgnifam_4.1.0.tar.gz"
)
# POOL_NONREDUNDANT_FAMILIES stageAs folder <- mgnifam output (per chunk)
STAGING = {
    "seed_msa_sto": "seed_msa/*.sto.gz",
    "full_msa_sto": "full_msa/*.sto.gz",
    "hmm": "hmm/*.hmm.gz",
    "rf": "rf/*.txt",
    "refined_families": "*_families.tsv",
    "discarded_clusters": "*_discarded.csv",
    "successful_clusters": "*_successful.txt",
    "converged_families": "*_converged.txt",
    "family_metadata": "*_metadata.csv",
    "family_reps": "*_reps.fasta.gz",
}

with tempfile.TemporaryDirectory() as tmp:
    d = Path(tmp)
    with tarfile.open(FIXTURE_TAR) as tar:
        tar.extractall(d, filter="data")
    FIXTURE = d / "mgnifam_4.1.0"
    for folder, pattern in STAGING.items():
        (d / "in" / folder).mkdir(parents=True)
        for chunk in ("000", "001"):
            for f in (FIXTURE / chunk).glob(pattern):
                shutil.copy(f, d / "in" / folder / f.name)
    (d / "redundant.txt").write_text("001_2\n")
    main(
        [
            "-d",
            str(d / "in"),
            "-o",
            str(d / "out"),
            "-r",
            str(d / "redundant.txt"),
            "-i",
            "5",
        ]
    )
    out = d / "out"

    # Keys sort by name; the redundant family is dropped before numbering from --iter
    assert json.loads((out / "family_to_id.json").read_text()) == {"000_1": 5, "000_2": 6, "001_1": 7}

    # Artifacts: decompressed, in the pooled names; HMM NAME = pooled id
    for folder, ext in (("hmm", ".hmm"), ("seed_msa_sto", ".sto"), ("full_msa_sto", ".sto"), ("rf", ".txt")):
        assert sorted(p.name for p in (out / folder).iterdir()) == [f"{i}{ext}" for i in (5, 6, 7)], folder
    assert (out / "hmm/7.hmm").read_text().splitlines()[1] == "NAME  7"
    assert (out / "seed_msa_sto/5.sto").read_text().startswith("# STOCKHOLM 1.0")
    with gzip.open(FIXTURE / "001/full_msa/001_1.sto.gz", "rt") as f:
        assert (out / "full_msa_sto/7.sto").read_text() == f.read()

    # Metadata: mgnifam's header, ids rewritten, every other field as mgnifam wrote it
    rows = (out / "family_metadata.csv").read_text().splitlines()
    assert rows[0] == METADATA_HEADER
    src = {
        f"{c}_{line.split(',', 1)[0]}": line.split(",", 1)[1]
        for c in ("000", "001")
        for line in (FIXTURE / c / f"{c}_metadata.csv").read_text().splitlines()[1:]
    }
    assert rows[1:] == [f"5,{src['000_1']}", f"6,{src['000_2']}", f"7,{src['001_1']}"], rows

    # Members: every row of the kept families, none of 001_2
    members = [line.split("\t") for line in (out / "refined_families.tsv").read_text().splitlines()]
    n = {
        k: sum(
            1
            for line in (FIXTURE / k[:3] / f"{k[:3]}_families.tsv").read_text().splitlines()
            if line.split("\t")[0] == k[4:]
        )
        for k in ("000_1", "000_2", "001_1")
    }
    assert [m[0] for m in members] == ["5"] * n["000_1"] + ["6"] * n["000_2"] + ["7"] * n["001_1"]

    # Converged: 000 lists 2, 001 lists 1 and 2 (2 redundant)
    assert (out / "converged_families.txt").read_text().split() == ["6", "7"]

    # Representatives: plain FASTA, family column = pooled id; ids-only FASTA for annotation
    reps = (out / "family_reps.fasta").read_text().splitlines()
    assert [h for h in reps if h.startswith(">")] == [
        ">782510898/278-392\t5",
        ">376297851/617-724\t6",
        ">5761513631/864-968\t7",
    ], reps
    assert (out / "family_ids.fasta").read_text().splitlines()[::2] == [">5", ">6", ">7"]

    # Clusters: headers dropped, rows kept
    assert (out / "discarded_clusters.txt").read_text() == ""
    assert len((out / "successful_clusters.txt").read_text().split()) == 4

    # A chunk in another format fails instead of being read under the wrong column names
    bad = d / "in/family_metadata/000_metadata.csv"
    bad.write_text("family_id,full_msa_size,protein,region,length,sequence,consensus,converged\n")
    try:
        main(
            [
                "-d",
                str(d / "in"),
                "-o",
                str(d / "out2"),
                "-r",
                str(d / "redundant.txt"),
                "-i",
                "5",
            ]
        )
    except ValueError as e:
        assert "header is not" in str(e), e
    else:
        raise AssertionError("legacy metadata header accepted")

print("ok")
