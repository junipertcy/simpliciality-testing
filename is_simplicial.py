import click
from SimplicialTest import SimplicialTest


@click.command()
@click.option('-k', '--degree_seq_file', 'degree_sequence', type=click.File('r'), help='Path to degree sequence file.')
@click.option('-s', '--size_seq_file', 'size_sequence', type=click.File('r'), help='Path to size sequence file.')
def is_simplicial(degree_sequence, size_sequence):
    degree_sequence = list(map(int, degree_sequence.read().replace("\n", "").split(" ")))
    size_sequence = list(map(int, size_sequence.read().replace("\n", "").split(" ")))
    st = SimplicialTest(degree_sequence, size_sequence)
    result = st.is_simplicial()
    if result is True:
        print(f"Yes, the joint sequence is simplicial. \nThe complex is: {st.identifier2facets(st.identifier)}")
    else:
        print("No, it cannot form a simplicial complex.")


if __name__ == "__main__":
    is_simplicial()
