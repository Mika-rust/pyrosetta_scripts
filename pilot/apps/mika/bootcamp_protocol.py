# bootcamp_protocol.py
from pyrosetta.rosetta.core.kinematics import FoldTree
from pyrosetta.rosetta.core.scoring.dssp import Dssp

def identify_secondary_structure_spans(ss):
    elements = []
    i = 0
    n = len(ss)
    while i < n:
        if ss[i] not in ('H', 'E'):
            i += 1
            continue
        start = i + 1
        t = ss[i]
        j = i
        while j < n and ss[j] == t:
            j += 1
        elements.append((start, j))
        i = j
    return elements

def fold_tree_from_dssp_string(ss):
    def mid(a, b): 
        return (a + b) // 2

    ss_spans = identify_secondary_structure_spans(ss)
    nres = len(ss)
    ft = FoldTree()
    if not ss_spans:
        return ft

    # build inter-SS loops
    inter_loops = []
    for (s1, e1), (s2, e2) in zip(ss_spans, ss_spans[1:]):
        if e1 + 1 <= s2 - 1:
            inter_loops.append((e1 + 1, s2 - 1))

    # trailing loop after last SS
    last_start, last_end = ss_spans[-1]
    trailing_loop = (last_end + 1, nres) if last_end + 1 <= nres else None

    # root at midpoint of first SS
    r_start, r_end = ss_spans[0]
    root = mid(r_start, r_end)

    peptide_edge = -1
    jump_id = 1
    seen = set()

    def add_edge(u, v, t):
        if u == v:
            return
        key = (u, v, t)
        if key in seen:
            return
        ft.add_edge(u, v, t)
        seen.add(key)

    # per handout: root -> 1 (back), root -> end(first SS) (forward)
    add_edge(root, 1, peptide_edge)
    add_edge(root, r_end, peptide_edge)

    # other SS midpoints: jump + fan peptide edges
    for a, b in ss_spans[1:]:
        m = mid(a, b)
        add_edge(root, m, jump_id); jump_id += 1
        add_edge(m, a, peptide_edge)
        add_edge(m, b, peptide_edge)

    # inter-SS loops
    for a, b in inter_loops:
        m = mid(a, b)
        add_edge(root, m, jump_id); jump_id += 1
        add_edge(m, a, peptide_edge)
        add_edge(m, b, peptide_edge)

    # trailing loop
    if trailing_loop:
        a, b = trailing_loop
        m = mid(a, b)
        add_edge(root, m, jump_id); jump_id += 1
        add_edge(m, a, peptide_edge)
        add_edge(m, b, peptide_edge)

    return ft

def fold_tree_from_pose(pose):
    dssp = Dssp(pose)
    ss = dssp.get_dssp_secstruct()
    return fold_tree_from_dssp_string(ss)

