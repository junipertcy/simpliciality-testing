from simplicial_test import *
import pickle

many_small_cases = pickle.load(open("datasets/multiple_small_test_cases.pickle", "rb"))


def test_batch_dataset():
    for degs, sizes in many_small_cases:
        st = Test(degs, sizes)
        is_simplicial, facets = st.is_simplicial()
        assert is_simplicial is True
        joint_seqs = compute_joint_seq(facets)
        assert if_facets_simplicial(facets) is True
        assert joint_seqs[0] == sorted(sizes, reverse=True)
        assert joint_seqs[1] == sorted(degs, reverse=True)
