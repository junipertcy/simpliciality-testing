import click
from SimplicialTest import *


@click.command()
@click.option('-k', '--degree_seq_file', 'degree_sequence', type=click.File('r'), help='Path to degree sequence file.')
@click.option('-s', '--size_seq_file', 'size_sequence', type=click.File('r'), help='Path to size sequence file.')
@click.option('--greedy/--no-greedy', default=False, help='Enable the Havel–Hakimi-type recursive algorithm.')
def is_simplicial(degree_sequence, size_sequence, greedy):
    degree_sequence = list(map(int, degree_sequence.read().replace("\n", "").split(" ")))
    size_sequence = list(map(int, size_sequence.read().replace("\n", "").split(" ")))
    st = SimplicialTest(degree_sequence, size_sequence)
    result = st.is_simplicial(greedy=greedy, preprocess=True)
    if result is True:
        print(f"Yes, the joint sequence is simplicial. \nThe complex is: {st.identifier2facets(st.identifier)}")
    else:
        print("No, it cannot form a simplicial complex.")


if __name__ == "__main__":
    is_simplicial()
