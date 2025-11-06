import pytest
import bootcamp_protocol as bp


def test_pytest_runs():
    assert 1 + 1 == 2


def test_identify_secondary_structure_spans_1():
    ss1 = "   EEEEE   HHHHHHHH  EEEEE   IGNOR EEEEEE   HHHHHHHHHHH  EEEEE  HHHH   "
    expected1 = [(4, 8), (12, 19), (22, 26), (36, 41), (45, 55), (58, 62), (65, 68)]
    assert bp.identify_secondary_structure_spans(ss1) == expected1


def test_identify_secondary_structure_spans_2():
    ss2 = "HHHHHHH   HHHHHHHHHHHH      HHHHHHHHHHHHEEEEEEEEEEHHHHHHH EEEEHHH "
    expected2 = [(1, 7), (11, 22), (29, 40), (41, 50), (51, 57), (59, 62), (63, 65)]
    assert bp.identify_secondary_structure_spans(ss2) == expected2


def test_identify_secondary_structure_spans_3():
    ss3 = "EEEEEEEEE EEEEEEEE EEEEEEEEE H EEEEE H H H EEEEEEEE"
    expected3 = [(1, 9), (11, 18), (20, 28), (30, 30), (32, 36), (38, 38), (40, 40), (42, 42), (44, 51)]
    assert bp.identify_secondary_structure_spans(ss3) == expected3


def _collect_edge_pairs_by_residue(ft, max_res):
    pairs = set()
    for res in range(1, max_res + 1):
        try:
            e = ft.get_residue_edge(res)  # raises on root
        except RuntimeError:
            continue
        if e is None:
            continue
        a, b = int(e.start()), int(e.stop())
        if a != b:
            pairs.add((a, b))
    return pairs


def test_fold_tree_from_dssp_string_example():
    #        000000000111111111122222222223333333333444444444455555555556666666666777777777788888888889999999999
    #        123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789
    ss = "   EEEEEEE    EEEEEEE         EEEEEEEEE    EEEEEEEEEE   HHHHHH         EEEEEEEEEHHHHHHHHH          "
    assert len(ss) == 99

    ft = bp.fold_tree_from_dssp_string(ss)

    expected_pairs = {
        (7, 1), (7, 10),
        (7, 12), (12, 11), (12, 14),
        (7, 18), (18, 15), (18, 21),
        (7, 26), (26, 22), (26, 30),
        (7, 35), (35, 31), (35, 39),
        (7, 41), (41, 40), (41, 43),
        (7, 48), (48, 44), (48, 53),
        (7, 55), (55, 54), (55, 56),
        (7, 59), (59, 57), (59, 62),
        (7, 67), (67, 63), (67, 71),
        (7, 76), (76, 72), (76, 80),
        (7, 85), (85, 81), (85, 89),
        (7, 92), (92, 90), (92, 99),
    }

    pairs = _collect_edge_pairs_by_residue(ft, max_res=len(ss))
    missing = expected_pairs - pairs
    assert not missing, f"Missing edges: {sorted(missing)}"
    assert len(pairs) == 38, f"Expected 38 edges, got {len(pairs)}"
