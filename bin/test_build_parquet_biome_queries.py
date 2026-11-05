#!/usr/bin/env python3
"""Self-check for build_parquet_biome_queries.py and the parse_biomes.py lineage mode.
Run: python3 bin/test_build_parquet_biome_queries.py"""

import os
import subprocess
import sys
import tempfile

import pyarrow as pa
import pyarrow.parquet as pq
from build_parquet_biome_queries import main

A, B = "root:Host-associated:Human", "root:Host-associated:Human:Skin"
with tempfile.TemporaryDirectory() as d:
    clusters = os.path.join(d, "clusters.parquet")
    pq.write_table(
        pa.table({"cluster_rep": [1, 2, 3, 4], "cluster_rep_biomes": [f"{A};{B};{B}", B, None, ""]}), clusters
    )
    families = os.path.join(d, "families.tsv")
    with open(families, "w") as fh:
        # member 1 twice (two regions): counted once; 3 and 4 have no biome
        fh.write("7\t1/1-50\n7\t1/60-120\n7\t2/5-80\n8\t3/1-40\n8\t4_10_90/2-30\n")
    out = os.path.join(d, "query_results")
    main(["--families", families, "--clusters", clusters, "--output_dir", out])
    assert open(f"{out}/7.tsv").read() == f"1\t{A};{B};{B}\n2\t{B}\n"
    assert open(f"{out}/8.tsv").read() == "3\t\n4\t\n"

    # A member that is not a cluster_rep fails
    with open(families, "a") as fh:
        fh.write("8\t5/1-40\n")
    try:
        main(["--families", families, "--clusters", clusters, "--output_dir", out])
        raise AssertionError("member 5 is not a cluster_rep and must fail")
    except SystemExit as e:
        assert "1 members are not cluster_rep rows" in str(e), e

    # Lineage mode: 1 per member per distinct lineage, then prefix sums
    res = os.path.join(d, "biome_results")
    script = os.path.join(os.path.dirname(os.path.abspath(__file__)), "parse_biomes.py")
    subprocess.run([sys.executable, script, "--query_results", out, "--output_dir", res], check=True)
    rows = open(f"{res}/7.csv").read().splitlines()
    assert rows[0] == "ids,labels,parents,counts"
    assert sorted(rows[1:]) == sorted(
        [
            "root,root,,3",
            "root:Host-associated,Host-associated,root,3",
            f"{A},Human,root:Host-associated,3",
            f"{B},Skin,{A},2",
        ]
    ), rows
    assert open(f"{res}/8.csv").read() == "ids,labels,parents,counts\n"
print("ok")
