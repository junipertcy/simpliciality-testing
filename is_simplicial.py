import click
from simplicial_test.Test import *


@click.command()
@click.option('-k', '--degree_seq_file', 'degree_sequence', type=click.File('r'), help='Path to degree sequence file.')
@click.option('-s', '--size_seq_file', 'size_sequence', type=click.File('r'), help='Path to size sequence file.')
def is_simplicial(degree_sequence, size_sequence):
    degree_sequence = list(map(int, degree_sequence.read().replace("\n", "").split(" ")))
    size_sequence = list(map(int, size_sequence.read().replace("\n", "").split(" ")))
    is_simplicial, facets = Test(degree_sequence, size_sequence).is_simplicial()
    if is_simplicial is True:
        joint_seqs = compute_joint_seq(facets)
        assert if_facets_simplicial(facets) is True
        assert joint_seqs[0] == sorted(size_sequence, reverse=True)
        assert joint_seqs[1] == sorted(degree_sequence, reverse=True)
        print(f"Yes, the joint sequence is simplicial. \nThe complex is: {facets}")
    else:
        print("No, it cannot form a simplicial complex.")


if __name__ == "__main__":
    is_simplicial()
