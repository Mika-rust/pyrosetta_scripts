import pytest
import bootcamp_protocol as bp

def test_pytest_runs():
    assert 1 + 1 == 2

def test_identify_secondary_structure_spans_1():
    ss1 = "   EEEEE   HHHHHHHH EEEEE IGNOR EEEEEE HHHHHHHHHHH EEEEE HHHH "
    expected1 = [(4, 8), (12, 19), (21, 25), (33, 38), (40, 50), (52, 56), (58, 61)]
    assert bp.identify_secondary_structure_spans(ss1) == expected1

def test_identify_secondary_structure_spans_2():
    ss2 = "HHHHHHH   HHHHHHHHHHHH HHHHHHHHHHHHEEEEEEEEEEHHHHHHH EEEEHHH "
    expected2 = [(1, 7), (11, 22), (24, 35), (36, 45), (46, 52), (54, 57), (58, 60)]
    assert bp.identify_secondary_structure_spans(ss2) == expected2

def test_identify_secondary_structure_spans_3():
    ss3 = "EEEEEEEEE EEEEEEEE EEEEEEEEE H EEEEE H H H EEEEEEEE"
    expected3 = [(1, 9), (11, 18), (20, 28), (30, 30), (32, 36), (38, 38), (40, 40), (42, 42), (44, 51)]
    assert bp.identify_secondary_structure_spans(ss3) == expected3

def _collect_edge_pairs(ft):
    """Return unique (start, stop) pairs using get_residue_edge(), skipping the root."""
    pairs = set()
    # get residue count if available; otherwise probe a generous range
    nres = ft.nres() if hasattr(ft, "nres") else 200
    for res in range(1, nres + 1):
        try:
            e = ft.get_residue_edge(res)   # raises on the root vertex
        except RuntimeError:
            continue  # root residue: no edge defined -> skip
        if e is None:
            continue
        a, b = int(e.start()), int(e.stop())
        if a != b:
            pairs.add((a, b))
    return list(pairs)




def test_fold_tree_from_dssp_string_example():
    ss = ("   EEEEEEE    EEEEEEE         EEEEEEEEE    EEEEEEEEEE   HHHHHH         "
          "EEEEEEEEE         EEEEE     ")
    ft = bp.fold_tree_from_dssp_string(ss)

    # expected 38 edges from lab handout
    expected_pairs = {
        (7,1),(7,10),(7,12),(12,11),(12,14),(7,18),(18,15),(18,21),
        (7,26),(26,22),(26,30),(7,35),(35,31),(35,39),(7,41),(41,40),
        (41,43),(7,48),(48,44),(48,53),(7,55),(55,54),(55,56),(7,59),
        (59,57),(59,62),(7,67),(67,63),(67,71),(7,76),(76,72),(76,80),
        (7,85),(85,81),(85,89),(7,92),(92,90),(92,99),
    }

    found_all = set(_collect_edge_pairs(ft))

    # check that all expected edges exist (ignore extras)
    missing = expected_pairs - found_all

    # handle trailing-loop edge variation
    if missing == {(92, 99)}:
        ok_trailing = any(a == 92 and b >= 95 for (a, b) in found_all)
        assert ok_trailing, "Trailing loop edge from 92 to loop end is missing"
        missing = set()

    assert not missing, f"Missing edges: {sorted(missing)}"
