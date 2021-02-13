#! /usr/bin/env python
# -*- coding: utf-8 -*-
#
# simplicial_test -- A greedy algorithm to check whether a joint sequence of degree/size distributions is simplicial.
#
# Copyright (C) 2020- Tzu-Chi Yen <tzuchi.yen@colorado.edu>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

"""``simplicial_test`` - Algorithm to check whether a joint sequence of degree/size distributions is simplicial
---------------------------------------------------------------------------------------------------------------

This module contains a greedy algorithm for checking simpliciality.

.. note::

   TODO.

"""
from simplicial_test.Test import *
from simplicial_test.utils import *
from simplicial_test.enumeration import *
from simplicial_test.sample import *
from simplicial_test.cm_count import *

__package__ = 'simplicial_test'
__title__ = 'simplicial_test: a recursive, deterministic algorithm to check whether a joint sequence of ' \
            'degree/size distributions is simplicial'
__description__ = ''
__copyright__ = 'Copyright 2020- Tzu-Chi Yen'
__license__ = "GPL version 3"
__author__ = """\n""".join([
    'Tzu-Chi Yen <tzuchi.yen@colorado.edu>',
])
__URL__ = "https://docs.netscied.tw/simplicialTest/index.html"
__version__ = '0.91.0'
__release__ = '0.91.0'

#  mark the class as not to be collected as pytest tests
Test.__test__ = False

__all__ = [
    "Test",
    "EnumRegistrar",
    "Enum",
    "get_hitting_sets",
    "compute_joint_seq",
    "flatten",
    "compute_dpv",
    "get_partition",
    "sort_facets",
    "get_relabeled_facets",
    "compute_entropy",
    "groupby_vtx_symm_class",
    "if_facets_simplicial",
    "gen_joint_sequence_from_sizes",
    "accel_asc",
    "NoMoreBalls",
    "paint_landscape",
    "paint_block",
    "__author__",
    "__URL__",
    "__version__",
    "__copyright__"
]
