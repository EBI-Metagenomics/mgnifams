#!/usr/bin/env python3
"""Self-check for parse_domains.py (#58, #59). Run: python3 bin/test_parse_domains.py"""

from parse_domains import calculate_mgnifam_start, construct_domain_architecture, extract_mgyp

cases = {
    "250671917": 1,
    "250671917/1_50": 1,
    "250671917_50_150": 50,
    "250671917_50_200/2_34": 51,
    "3387826881/357-470": 357,
    "250671917_50_200/2-34": 51,
}
for name, start in cases.items():
    assert calculate_mgnifam_start(name) == start, (name, calculate_mgnifam_start(name))
    assert extract_mgyp(name) == name.split("/")[0].split("_")[0]

# Pfam tuples are [acc, i_evalue, score, hmm_from, hmm_to, env_from, env_to]; order by env_from, not hmm_from
pfams = "[['PF00001', 1e-5, 30.0, 1, 50, 200, 260], ['PF00002', 1e-5, 30.0, 100, 150, 10, 60]]"
assert construct_domain_architecture(pfams, 7, [100]) == "PF00002\t7\tPF00001"
print("ok")
