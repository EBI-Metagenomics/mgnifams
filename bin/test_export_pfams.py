#!/usr/bin/env python3
"""Self-check: export_pfams.py template-coverage filter (#69). Run: python3 bin/test_export_pfams.py"""

import tempfile
from fractions import Fraction
from pathlib import Path

from export_pfams import parse_summary_block, template_coverage

# hhsearch summary block, fixed columns as HH-suite writes them
HHR = """Query         1
 No Hit                             Prob E-value P-value  Score    SS Cols Query HMM  Template HMM
  1 PF00001.24 ; 7tm_1 ; 7 transm 100.0 1.1E-30 1.2E-35  200.0   0.0   80    1-80     11-90 (100)
  2 PF00002.29 ; 7tm_2 ; 7 transm  99.9 2.2E-20 2.3E-25  150.0   0.0   79    1-79     11-89 (100)

No 1
"""

assert template_coverage("11-90 (100)") == Fraction(4, 5)
assert template_coverage("58-470 (471)") == Fraction(413, 471)

with tempfile.TemporaryDirectory() as tmp:
    hhr = Path(tmp) / "x.hhr"
    hhr.write_text(HHR)
    every = parse_summary_block(hhr, {})
    assert [r["pfam"] for r in every] == ["PF00001", "PF00002"], every
    assert every[0]["template_hmm"] == "11-90 (100)"
    # 80% of the template is kept at 0.8, 79% is not
    kept = parse_summary_block(hhr, {}, Fraction("0.8"))
    assert [r["pfam"] for r in kept] == ["PF00001"], kept

print("test_export_pfams.py: all passed")
