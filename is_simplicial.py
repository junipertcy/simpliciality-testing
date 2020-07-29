import click
from SimplicialTest import SimplicialTest


@click.command()
@click.argument('degree_sequence', type=click.File('r'))
@click.argument('size_sequence', type=click.File('r'))
def is_simplicial(degree_sequence, size_sequence):
    degree_sequence = list(map(int, degree_sequence.read().replace("\n", "").split(" ")))
    size_sequence = list(map(int, size_sequence.read().replace("\n", "").split(" ")))
    st = SimplicialTest(degree_sequence, size_sequence)
    result = st.is_simplicial()
    if result is True:
        print("Yes, it forms a simplicial complex.")
        print(st.identifier2facets(st.identifier))
    else:
        print("No, it cannot form a simplicial complex.")


if __name__ == "__main__":
    is_simplicial()
