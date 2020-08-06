from SimplicialTest import *


def test_simple_1():
    size_list, degree_list = ([6, 6, 3, 3, 2], [3, 2, 2, 2, 3, 3, 2, 3])
    st = SimplicialTest(degree_list, size_list)
    assert st.is_simplicial(greedy=True) is True


def test_simple_2():
    size_list, degree_list = ([5, 4, 3, 3, 2], [4, 3, 3, 2, 2, 1, 1, 1])  # originally needs "special" "greedy" strategy
    st = SimplicialTest(degree_list, size_list)
    assert st.is_simplicial(greedy=True) is True


def test_simple_3():
    size_list, degree_list =  ([6, 4, 4], [3, 2, 3, 1, 1, 2, 2])  # something wrong with the forward strategyv (recursive case)
    st = SimplicialTest(degree_list, size_list)
    assert st.is_simplicial(greedy=True) is True


def test_simple_4():
    degree_list, size_list = ([1, 3, 1, 3, 4, 1, 3, 1], [5, 4, 3, 3, 2])  # another recursive case
    st = SimplicialTest(degree_list, size_list)
    assert st.is_simplicial(greedy=True) is True


def test_simple_5():
    degree_list, size_list = ([3, 3, 3, 1, 1, 1, 1, 1, 1], [6, 3, 3, 3])  # another recursive case
    st = SimplicialTest(degree_list, size_list)
    assert st.is_simplicial(greedy=True) is True


def test_simple_6():
    degree_list, size_list = ([2, 2, 2, 1, 2, 1], [5, 3, 2])
    st = SimplicialTest(degree_list, size_list)
    assert st.is_simplicial(greedy=True) is True


def test_simple_7():
    size_list, degree_list = ([8, 6, 6, 5, 4], [4, 4, 4, 3, 3, 3, 3, 2, 2, 1])  # TODO: RecursionError: maximum recursion depth exceeded while calling a Python object
    st = SimplicialTest(degree_list, size_list)
    assert st.is_simplicial(greedy=True) is True


def test_simple_8():
    degree_list, size_list = ([4, 4, 3, 3, 3, 2, 2], [6, 5, 4, 4, 2])  # NotImplementedError
    st = SimplicialTest(degree_list, size_list)
    assert st.is_simplicial(greedy=True) is True


def test_simple_9():
    degree_list = [3, 2, 2, 2, 2]
    size_list = [4, 3, 2, 2]
    st = SimplicialTest(degree_list, size_list)
    assert st.is_simplicial(greedy=True) is True
