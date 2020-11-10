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

"""``simplicial_test`` - Greedy algorithm to check whether a joint sequence of degree/size distributions is simplicial
---------------------------------------------------------------------------------------------------------------------

This module contains a greedy algorithm for checking simpliciality.

.. note::

   TODO.

"""
from simplicial_test.Test import *
from simplicial_test.utils import *
from simplicial_test.sample import *


__package__ = 'simplicial_test'
__title__ = 'simplicial_test: A greedy algorithm to check whether a joint sequence of degree/size distributions is simplicial'
__description__ = ''
__copyright__ = 'Copyright 2020- Tzu-Chi Yen'
__author__ = """\n""".join([
    'Tzu-Chi Yen <tzuchi.yen@colorado.edu>',
])
__URL__ = "https://docs.netscied.tw/simplicialTest/index.html"
__version__ = '0.90.0'
__release__ = '0.90.0'

#  mark the class as not to be collected as pytest tests
Test.__test__ = False

__all__ = [
    "Test",
    "compute_joint_seq",
    "if_facets_simplicial",
    "accel_asc",
    "__author__",
    "__URL__",
    "__version__",
    "__copyright__"
]
