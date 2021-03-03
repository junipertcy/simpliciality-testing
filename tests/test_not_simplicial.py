from simplicial_test import *


def test_not_simplicial_1():
    degree_list, size_list = ([2], [2])
    st = Test(degree_list, size_list)
    assert st.is_simplicial()[0] is False


def test_not_simplicial_2():
    degree_list, size_list = ([3, 3, 3, 2, 1, 1, 1], [6, 6, 2])
    st = Test(degree_list, size_list)
    assert st.is_simplicial()[0] is False


def test_not_simplicial_3():
    degree_list, size_list = ([3, 3, 3, 2, 1, 1, 1, 1, 1, 1], [6, 6, 2, 1, 1, 1])
    st = Test(degree_list, size_list)
    assert st.is_simplicial()[0] is False


def test_not_simplicial_4():
    degree_list, size_list = ([6, 3, 3, 3, 2, 1, 1, 1, 1, 1, 1], [7, 7, 3, 2, 2, 2])
    st = Test(degree_list, size_list)
    assert st.is_simplicial()[0] is False


def test_not_simplicial_5():
    degree_list, size_list = (
        [6, 3, 3, 3, 2, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1], [7, 7, 3, 2, 2, 2, 1, 1, 1, 1, 1, 1, 1])
    st = Test(degree_list, size_list)
    assert st.is_simplicial()[0] is False


def test_not_simplicial_unequal_seq_sum():
    degree_list, size_list = ([3, 3, 3, 2, 1, 1], [6, 6, 2])
    st = Test(degree_list, size_list)
    assert st.is_simplicial()[0] is False


def test_cutoff():
    degree_list, size_list = ([12, 8, 5, 3, 1, 1, 1, 1, 1, 1, 1, 1], [3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3])
    st = Test(degree_list, size_list, verbose=False, width=1e4, cutoff=3000)
    assert st.is_simplicial()[0] is False
