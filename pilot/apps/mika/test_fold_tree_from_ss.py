import pytest

def test_pytest_runs():
    assert 1 + 1 == 2

def identify_secondary_structure_spans(ss):
    elements = []
    i = 0
    n = len(ss)
    
    while i < n:
        # Skip non-H/E characters
        if ss[i] not in ['H', 'E']:
            i += 1
            continue
            
        # Found start of SS element
        start = i + 1  # Convert to 1-based indexing
        ss_type = ss[i]
        
        # Find the end of this SS element
        j = i
        while j < n and ss[j] == ss_type:
            j += 1
            
        end = j  # j is now at the first non-matching character
        
        # Add to elements
        elements.append((start, end))
        
        # Move i to the end of this element
        i = j
    
    return elements





def test_identify_secondary_structure_spans_1():
    ss1 = "   EEEEE   HHHHHHHH EEEEE IGNOR EEEEEE HHHHHHHHHHH EEEEE HHHH "
    expected1 = [(4, 8), (12, 19), (21, 25), (33, 38), (40, 50), (52, 56), (58, 61)]
    result = identify_secondary_structure_spans(ss1)
    assert result == expected1

def test_identify_secondary_structure_spans_2():
    ss2 = "HHHHHHH   HHHHHHHHHHHH HHHHHHHHHHHHEEEEEEEEEEHHHHHHH EEEEHHH "
    expected2 = [(1, 7), (11, 22), (24, 35), (36, 45), (46, 52), (54, 57), (58, 60)]
    result = identify_secondary_structure_spans(ss2)
    assert result == expected2

def test_identify_secondary_structure_spans_3():
    ss3 = "EEEEEEEEE EEEEEEEE EEEEEEEEE H EEEEE H H H EEEEEEEE"
    expected3 = [(1, 9), (11, 18), (20, 28), (30, 30), (32, 36), (38, 38), (40, 40), (42, 42), (44, 51)]
    result = identify_secondary_structure_spans(ss3)
    assert result == expected3

