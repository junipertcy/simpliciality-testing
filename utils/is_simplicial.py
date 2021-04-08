#! /usr/bin/env python
# -*- coding: utf-8 -*-
#
# simplicial_test -- a python module to realize simplicial joint sequences
#
# Copyright (C) 2020-2021 Tzu-Chi Yen <tzuchi.yen@colorado.edu>
#
# This program is free software; you can redistribute it and/or modify it under
# the terms of the GNU Lesser General Public License as published by the Free
# Software Foundation; either version 3 of the License, or (at your option) any
# later version.
#
# This program is distributed in the hope that it will be useful, but WITHOUT
# ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
# FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public License for more
# details.
#
# You should have received a copy of the GNU Lesser General Public License
# along with this program. If not, see <http://www.gnu.org/licenses/>.


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
