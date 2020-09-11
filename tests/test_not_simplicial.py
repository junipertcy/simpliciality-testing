from SimplicialTest import *


def test_not_simplicial_1():
    degree_list, size_list = ([2], [2])
    st = SimplicialTest(degree_list, size_list)
    assert st.is_simplicial(greedy=True) is False

