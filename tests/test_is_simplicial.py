from SimplicialTest import *


def test_simple_1():
    size_list, degree_list = ([6, 6, 3, 3, 2], [3, 2, 2, 2, 3, 3, 2, 3])
    st = SimplicialTest(degree_list, size_list)
    assert st.is_simplicial(greedy=True) is True


def test_simple_2():
    size_list, degree_list = ([5, 4, 3, 3, 2], [4, 3, 3, 2, 2, 1, 1, 1])
    st = SimplicialTest(degree_list, size_list)
    assert st.is_simplicial(greedy=True) is True


def test_simple_3():
    size_list, degree_list =  ([6, 4, 4], [3, 2, 3, 1, 1, 2, 2])
    st = SimplicialTest(degree_list, size_list)
    assert st.is_simplicial(greedy=True) is True


def test_simple_4():
    degree_list, size_list = ([1, 3, 1, 3, 4, 1, 3, 1], [5, 4, 3, 3, 2])
    st = SimplicialTest(degree_list, size_list)
    assert st.is_simplicial(greedy=True) is True


def test_simple_5():
    degree_list, size_list = ([3, 3, 3, 1, 1, 1, 1, 1, 1], [6, 3, 3, 3])
    st = SimplicialTest(degree_list, size_list)
    assert st.is_simplicial(greedy=True) is True


def test_simple_6():
    degree_list, size_list = ([2, 2, 2, 1, 2, 1], [5, 3, 2])
    st = SimplicialTest(degree_list, size_list)
    assert st.is_simplicial(greedy=True) is True


def test_simple_7():
    size_list, degree_list = ([8, 6, 6, 5, 4], [4, 4, 4, 3, 3, 3, 3, 2, 2, 1])
    st = SimplicialTest(degree_list, size_list)
    assert st.is_simplicial(greedy=True) is True


def test_simple_8():
    degree_list, size_list = ([4, 4, 3, 3, 3, 2, 2], [6, 5, 4, 4, 2])
    st = SimplicialTest(degree_list, size_list)
    assert st.is_simplicial(greedy=True) is True


def test_simple_9():
    degree_list = [3, 2, 2, 2, 2]
    size_list = [4, 3, 2, 2]
    st = SimplicialTest(degree_list, size_list)
    assert st.is_simplicial(greedy=True) is True


def test_simple_10():
    degree_list = [4, 3, 3, 2, 2, 2, 2, 1]
    size_list = [6, 5, 4, 2, 2]
    st = SimplicialTest(degree_list, size_list)
    assert st.is_simplicial(greedy=True) is True


def test_simple_11():
    degree_list = [4, 3, 2, 2, 2, 2, 2, 2]
    size_list = [5, 5, 4, 3, 2]
    st = SimplicialTest(degree_list, size_list)
    assert st.is_simplicial(greedy=True) is True


def test_simple_12():
    degree_list = [11, 8, 8, 8, 7, 7, 6, 6, 6, 5, 5, 5, 5, 4, 4]
    size_list = [14, 12, 11, 9, 8, 8, 8, 8, 5, 5, 4, 3]
    st = SimplicialTest(degree_list, size_list)
    assert st.is_simplicial(greedy=True) is True


def test_simple_12_reduced():
    """TODO: a new kind of reduction yet to be written"""
    degree_list = [7, 7, 7, 6, 6, 5, 5, 5, 4, 4, 4, 4, 3, 3]
    size_list = [11, 10, 8, 7, 7, 7, 7, 4, 4, 3, 2]
    st = SimplicialTest(degree_list, size_list)
    assert st.is_simplicial(greedy=True) is True


def test_simple_13():
    degree_list = [8, 7, 7, 7, 6, 6, 5, 5, 5, 5, 5, 4, 4]
    size_list = [8, 7, 7, 7, 7, 7, 6, 6, 6, 5, 5, 3]
    st = SimplicialTest(degree_list, size_list)
    assert st.is_simplicial(greedy=True) is True


def test_simple_14():
    degree_list = [3, 2, 2, 2, 2, 1]
    size_list = [4, 4, 2, 2]
    st = SimplicialTest(degree_list, size_list)
    assert st.is_simplicial(greedy=True) is True


def test_cond_3():
    """Testing that cond_3 works."""
    degree_list, size_list = ([3, 3, 3, 3, 3, 3, 3, 2, 2], [6, 6, 6, 5, 2])
    st = SimplicialTest(degree_list, size_list)
    assert st.is_simplicial(greedy=True) is True


def test_sample_icebreaker_q_be_0():
    """Testing that sample_icebreaker should enforce q=0 works."""
    degree_list, size_list = ([4, 4, 3, 3, 3, 2, 2, 1, 1, 1], [8, 6, 4, 3, 3])
    st = SimplicialTest(degree_list, size_list)
    assert st.is_simplicial(greedy=True) is True


def test_cond_4():
    """Testing that cond_4 works."""
    degree_list, size_list = ([4, 3, 3, 3, 3, 3, 3, 2, 2], [6, 6, 6, 6, 2])
    st = SimplicialTest(degree_list, size_list)
    assert st.is_simplicial(greedy=True) is True


def test_preprocessing():
    """Testing that preprocessing works."""
    # TODO: one should also check that without preprocessing, the algorithm should still work.
    degree_list, size_list = ([5, 4, 3, 2, 2, 2, 2, 2, 1], [7, 6, 4, 4, 2])
    st = SimplicialTest(degree_list, size_list)
    assert st.is_simplicial(greedy=True, preprocess=True) is True


def test_crime_dataset():
    degree_list = list(map(int, "1 1 1 2 1 1 1 1 1 1 1 1 1 5 1 1 1 1 1 1 2 1 4 2 1 1 1 1 1 1 1 2 1 2 1 1 1 1 8 3 1 1 2 1 2 6 5 2 2 1 1 3 2 1 2 1 1 2 3 1 1 1 1 1 1 1 1 2 1 2 2 4 1 1 4 1 3 1 1 1 4 2 1 1 1 1 1 1 1 3 2 1 1 3 7 2 1 2 1 1 2 2 1 2 1 1 2 2 2 14 1 3 1 1 2 1 4 1 3 1 1 1 1 2 2 2 1 1 2 2 2 2 2 2 2 1 1 1 1 3 1 1 1 1 1 3 1 2 1 1 2 1 2 1 2 2 1 2 1 3 1 1 1 1 1 1 2 1 3 1 1 1 3 2 1 1 2 1 1 1 1 1 2 1 3 1 1 2 5 1 1 1 2 2 2 3 2 2 2 2 1 1 1 2 1 1 1 1 1 1 1 1 1 1 1 1 1 3 1 2 1 2 2 1 1 1 1 1 1 1 1 2 1 1 2 1 1 1 1 2 1 1 1 1 1 1 1 2 2 1 2 1 1 2 1 2 1 1 2 1 1 2 4 2 1 3 2 1 1 1 1 1 1 3 1 1 3 2 1 2 1 2 1 3 1 1 2 1 1 1 1 2 1 1 1 1 1 1 1 1 1 1 2 1 1 1 2 2 4 2 1 1 2 2 1 2 1 1 1 2 1 3 1 1 1 1 1 1 1 1 1 1 2 2 1 1 1 1 1 1 1 1 1 1 1 1 2 2 1 1 2 3 3 1 1 2 1 2 1 2 1 2 1 2 2 1 1 2 1 1 1 1 1 1 1 2 1 1 1 1 1 1 2 2 1 1 2 2 1 1 1 1 1 1 1 2 1 1 1 1 2 2 1 2 2 1 1 1 1 2 1 1 2 1 2 3 4 2 4 1 1 1 1 2 1 1 1 1 1 1 1 5 1 1 1 1 1 1 1 1 1 2 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 2 3 1 1 1 1 1 3 1 2 1 1 1 1 2 1 1 1 1 2 3 2 1 1 1 1 1 1 1 2 2 1 1 1 1 2 1 1 1 1 1 1 1 1 1 1 1 1 2 1 1 1 1 2 2 2 1 2 1 2 1 1 1 1 1 2 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1".split(" ")))
    size_list = list(map(int, "25 22 18 17 14 12 11 11 11 10 10 9 9 9 9 9 9 9 8 8 8 8 8 7 7 7 7 7 7 7 7 7 7 6 6 6 6 6 6 6 6 6 6 6 6 5 5 5 5 5 5 5 5 5 5 5 5 5 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 1 1 1 1 1 1 1 1 1 1 1 1".split(" ")))
    st = SimplicialTest(degree_list, size_list)
    assert st.is_simplicial(greedy=True) is True

