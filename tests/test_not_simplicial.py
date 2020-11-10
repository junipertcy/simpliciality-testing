from simplicial_test import *


def test_not_simplicial_1():
    degree_list, size_list = ([2], [2])
    st = Test(degree_list, size_list)
    assert st.is_simplicial()[0] is False
