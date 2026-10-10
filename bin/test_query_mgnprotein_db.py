#!/usr/bin/env python3
"""Self-check for query_mgnprotein_db.py against a fake psycopg2. Run: python3 bin/test_query_mgnprotein_db.py"""

import sys
import tempfile
import types

# mgyp -> (metadata->'b', metadata->'p'); 3 has no 'p' key, 9 is not in the DB
DB = {1: ([10, 11], [["PF00001", 1e-5]]), 2: ([12], [["PF00002", 1e-3]]), 3: ([13], None)}
queries = []


class FakeCursor:
    def execute(self, sql, params):
        queries.append(params[0])
        self.rows = [(int(s), *DB[int(s)]) for s in params[0] if int(s) in DB]

    def fetchall(self):
        return self.rows

    def close(self):
        pass


class FakeConnection:
    def cursor(self):
        return FakeCursor()

    def close(self):
        pass


sys.modules["psycopg2"] = types.SimpleNamespace(connect=lambda **_: FakeConnection())
import query_mgnprotein_db  # noqa: E402

query_mgnprotein_db.BATCH_SIZE = 2
with tempfile.TemporaryDirectory() as tmp:
    with open(f"{tmp}/families.tsv", "w") as f:
        f.write("famA\t1/1-50\nfamA\t1_5_60/2-30\nfamA\t2/3-40\nfamB\t3/1-20\nfamB\t9/1-20\nfamC\t9/1-20\n")
    query_mgnprotein_db.query_sequence_explorer_protein({}, [f"{tmp}/families.tsv"], tmp, threads=2)

    def rows(family):
        return sorted(open(f"{tmp}/{family}.tsv").read().splitlines())

    assert rows("famA") == ["1\t[10, 11]\t[['PF00001', 1e-05]]", "2\t[12]\t[['PF00002', 0.001]]"], rows("famA")
    assert rows("famB") == ["3\t[13]\t"], rows("famB")
    assert rows("famC") == [], rows("famC")
# Whole families per chunk (famA, famB, famC each reach BATCH_SIZE rows), each member queried once per chunk
assert sorted(map(sorted, queries)) == [["1", "2"], ["3", "9"], ["9"]], queries
print("ok")
