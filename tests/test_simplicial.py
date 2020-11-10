from simplicial_test import *


def test_simplicial_1():
    size_list, degree_list = ([6, 6, 3, 3, 2], [3, 2, 2, 2, 3, 3, 2, 3])
    st = Test(degree_list, size_list)
    is_simplicial, facets = st.is_simplicial()
    assert is_simplicial is True
    joint_seqs = compute_joint_seq(facets)
    assert if_facets_simplicial(facets) is True
    assert joint_seqs[0] == sorted(size_list, reverse=True)
    assert joint_seqs[1] == sorted(degree_list, reverse=True)


def test_simplicial_2():
    size_list, degree_list = ([5, 4, 3, 3, 2], [4, 3, 3, 2, 2, 1, 1, 1])
    st = Test(degree_list, size_list)
    is_simplicial, facets = st.is_simplicial()
    assert is_simplicial is True
    joint_seqs = compute_joint_seq(facets)
    assert if_facets_simplicial(facets) is True
    assert joint_seqs[0] == sorted(size_list, reverse=True)
    assert joint_seqs[1] == sorted(degree_list, reverse=True)


def test_simplicial_3():
    size_list, degree_list = ([6, 4, 4], [3, 2, 3, 1, 1, 2, 2])
    st = Test(degree_list, size_list)
    is_simplicial, facets = st.is_simplicial()
    assert is_simplicial is True
    joint_seqs = compute_joint_seq(facets)
    assert if_facets_simplicial(facets) is True
    assert joint_seqs[0] == sorted(size_list, reverse=True)
    assert joint_seqs[1] == sorted(degree_list, reverse=True)


def test_simplicial_4():
    degree_list, size_list = ([1, 3, 1, 3, 4, 1, 3, 1], [5, 4, 3, 3, 2])
    st = Test(degree_list, size_list)
    is_simplicial, facets = st.is_simplicial()
    assert is_simplicial is True
    joint_seqs = compute_joint_seq(facets)
    assert if_facets_simplicial(facets) is True
    assert joint_seqs[0] == sorted(size_list, reverse=True)
    assert joint_seqs[1] == sorted(degree_list, reverse=True)


def test_simplicial_5():
    degree_list, size_list = ([3, 3, 3, 1, 1, 1, 1, 1, 1], [6, 3, 3, 3])
    st = Test(degree_list, size_list)
    is_simplicial, facets = st.is_simplicial()
    assert is_simplicial is True
    joint_seqs = compute_joint_seq(facets)
    assert if_facets_simplicial(facets) is True
    assert joint_seqs[0] == sorted(size_list, reverse=True)
    assert joint_seqs[1] == sorted(degree_list, reverse=True)


def test_simplicial_6():
    degree_list, size_list = ([2, 2, 2, 1, 2, 1], [5, 3, 2])
    st = Test(degree_list, size_list)
    is_simplicial, facets = st.is_simplicial()
    assert is_simplicial is True
    joint_seqs = compute_joint_seq(facets)
    assert if_facets_simplicial(facets) is True
    assert joint_seqs[0] == sorted(size_list, reverse=True)
    assert joint_seqs[1] == sorted(degree_list, reverse=True)


def test_simplicial_7():
    size_list, degree_list = ([8, 6, 6, 5, 4], [4, 4, 4, 3, 3, 3, 3, 2, 2, 1])
    st = Test(degree_list, size_list)
    is_simplicial, facets = st.is_simplicial()
    assert is_simplicial is True
    joint_seqs = compute_joint_seq(facets)
    assert if_facets_simplicial(facets) is True
    assert joint_seqs[0] == sorted(size_list, reverse=True)
    assert joint_seqs[1] == sorted(degree_list, reverse=True)


def test_simplicial_8():
    degree_list, size_list = ([4, 4, 3, 3, 3, 2, 2], [6, 5, 4, 4, 2])
    st = Test(degree_list, size_list)
    is_simplicial, facets = st.is_simplicial()
    assert is_simplicial is True
    joint_seqs = compute_joint_seq(facets)
    assert if_facets_simplicial(facets) is True
    assert joint_seqs[0] == sorted(size_list, reverse=True)
    assert joint_seqs[1] == sorted(degree_list, reverse=True)


def test_simplicial_9():
    degree_list = [3, 2, 2, 2, 2]
    size_list = [4, 3, 2, 2]
    st = Test(degree_list, size_list)
    is_simplicial, facets = st.is_simplicial()
    assert is_simplicial is True
    joint_seqs = compute_joint_seq(facets)
    assert if_facets_simplicial(facets) is True
    assert joint_seqs[0] == sorted(size_list, reverse=True)
    assert joint_seqs[1] == sorted(degree_list, reverse=True)


def test_simplicial_10():
    degree_list = [4, 3, 3, 2, 2, 2, 2, 1]
    size_list = [6, 5, 4, 2, 2]
    st = Test(degree_list, size_list)
    is_simplicial, facets = st.is_simplicial()
    assert is_simplicial is True
    joint_seqs = compute_joint_seq(facets)
    assert if_facets_simplicial(facets) is True
    assert joint_seqs[0] == sorted(size_list, reverse=True)
    assert joint_seqs[1] == sorted(degree_list, reverse=True)


def test_simplicial_11():
    degree_list = [4, 3, 2, 2, 2, 2, 2, 2]
    size_list = [5, 5, 4, 3, 2]
    st = Test(degree_list, size_list)
    is_simplicial, facets = st.is_simplicial()
    assert is_simplicial is True
    joint_seqs = compute_joint_seq(facets)
    assert if_facets_simplicial(facets) is True
    assert joint_seqs[0] == sorted(size_list, reverse=True)
    assert joint_seqs[1] == sorted(degree_list, reverse=True)


def test_simplicial_12():
    degree_list = [8, 7, 7, 7, 6, 6, 5, 5, 5, 5, 5, 4, 4]
    size_list = [8, 7, 7, 7, 7, 7, 6, 6, 6, 5, 5, 3]
    st = Test(degree_list, size_list)
    is_simplicial, facets = st.is_simplicial()
    assert is_simplicial is True
    joint_seqs = compute_joint_seq(facets)
    assert if_facets_simplicial(facets) is True
    assert joint_seqs[0] == sorted(size_list, reverse=True)
    assert joint_seqs[1] == sorted(degree_list, reverse=True)


def test_simplicial_13():
    degree_list = [3, 2, 2, 2, 2, 1]
    size_list = [4, 4, 2, 2]
    st = Test(degree_list, size_list)
    is_simplicial, facets = st.is_simplicial()
    assert is_simplicial is True
    joint_seqs = compute_joint_seq(facets)
    assert if_facets_simplicial(facets) is True
    assert joint_seqs[0] == sorted(size_list, reverse=True)
    assert joint_seqs[1] == sorted(degree_list, reverse=True)


def test_simplicial_14():
    degree_list, size_list = ([3, 3, 3, 3, 3, 3, 3, 2, 2], [6, 6, 6, 5, 2])
    st = Test(degree_list, size_list)
    is_simplicial, facets = st.is_simplicial()
    assert is_simplicial is True
    joint_seqs = compute_joint_seq(facets)
    assert if_facets_simplicial(facets) is True
    assert joint_seqs[0] == sorted(size_list, reverse=True)
    assert joint_seqs[1] == sorted(degree_list, reverse=True)


def test_simplicial_15():
    degree_list, size_list = ([4, 4, 3, 3, 3, 2, 2, 1, 1, 1], [8, 6, 4, 3, 3])
    st = Test(degree_list, size_list)
    is_simplicial, facets = st.is_simplicial()
    assert is_simplicial is True
    joint_seqs = compute_joint_seq(facets)
    assert if_facets_simplicial(facets) is True
    assert joint_seqs[0] == sorted(size_list, reverse=True)
    assert joint_seqs[1] == sorted(degree_list, reverse=True)


def test_simplicial_16():
    degree_list, size_list = ([4, 3, 3, 3, 3, 3, 3, 2, 2], [6, 6, 6, 6, 2])
    st = Test(degree_list, size_list)
    is_simplicial, facets = st.is_simplicial()
    assert is_simplicial is True
    joint_seqs = compute_joint_seq(facets)
    assert if_facets_simplicial(facets) is True
    assert joint_seqs[0] == sorted(size_list, reverse=True)
    assert joint_seqs[1] == sorted(degree_list, reverse=True)


def test_simplicial_17():
    degree_list, size_list = ([5, 4, 3, 2, 2, 2, 2, 2, 1], [7, 6, 4, 4, 2])
    st = Test(degree_list, size_list)
    is_simplicial, facets = st.is_simplicial()
    assert is_simplicial is True
    joint_seqs = compute_joint_seq(facets)
    assert if_facets_simplicial(facets) is True
    assert joint_seqs[0] == sorted(size_list, reverse=True)
    assert joint_seqs[1] == sorted(degree_list, reverse=True)


def test_simplicial_18():
    degree_list, size_list = (
        [11, 9, 8, 8, 7, 7, 7, 7, 7, 7, 7, 6, 6, 6, 6, 6, 6, 6, 5, 4, 4, 4],
        [21, 14, 13, 13, 12, 12, 12, 11, 10, 9, 9, 8]
    )
    st = Test(degree_list, size_list)
    is_simplicial, facets = st.is_simplicial()
    assert is_simplicial is True
    joint_seqs = compute_joint_seq(facets)
    assert if_facets_simplicial(facets) is True
    assert joint_seqs[0] == sorted(size_list, reverse=True)
    assert joint_seqs[1] == sorted(degree_list, reverse=True)


def test_simplicial_19():
    degree_list, size_list = (
        [11, 11, 9, 8, 8, 8, 8, 8, 7, 7, 6, 6, 6, 6, 6, 5, 5, 5, 5, 5, 5, 5, 4, 3],
        [23, 17, 15, 15, 15, 14, 13, 11, 10, 9, 8, 7]
    )
    st = Test(degree_list, size_list)
    is_simplicial, facets = st.is_simplicial()
    assert is_simplicial is True
    joint_seqs = compute_joint_seq(facets)
    assert if_facets_simplicial(facets) is True
    assert joint_seqs[0] == sorted(size_list, reverse=True)
    assert joint_seqs[1] == sorted(degree_list, reverse=True)


def test_simplicial_20():
    degree_list, size_list = (
        [11, 9, 8, 8, 8, 8, 8, 8, 8, 7, 7, 7, 6, 6, 6, 5, 5],
        [16, 15, 13, 13, 11, 11, 10, 10, 9, 9, 4, 4]
    )
    st = Test(degree_list, size_list)
    is_simplicial, facets = st.is_simplicial()
    assert is_simplicial is True
    joint_seqs = compute_joint_seq(facets)
    assert if_facets_simplicial(facets) is True
    assert joint_seqs[0] == sorted(size_list, reverse=True)
    assert joint_seqs[1] == sorted(degree_list, reverse=True)


def test_simplicial_21():
    degree_list = [11, 8, 8, 8, 7, 7, 6, 6, 6, 5, 5, 5, 5, 4, 4]
    size_list = [14, 12, 11, 9, 8, 8, 8, 8, 5, 5, 4, 3]
    st = Test(degree_list, size_list)
    is_simplicial, facets = st.is_simplicial()
    assert is_simplicial is True
    joint_seqs = compute_joint_seq(facets)
    assert if_facets_simplicial(facets) is True
    assert joint_seqs[0] == sorted(size_list, reverse=True)
    assert joint_seqs[1] == sorted(degree_list, reverse=True)


def test_simplicial_22():
    degree_list = [7, 7, 7, 6, 6, 5, 5, 5, 4, 4, 4, 4, 3, 3]
    size_list = [11, 10, 8, 7, 7, 7, 7, 4, 4, 3, 2]
    st = Test(degree_list, size_list)
    is_simplicial, facets = st.is_simplicial()
    assert is_simplicial is True
    joint_seqs = compute_joint_seq(facets)
    assert if_facets_simplicial(facets) is True
    assert joint_seqs[0] == sorted(size_list, reverse=True)
    assert joint_seqs[1] == sorted(degree_list, reverse=True)
