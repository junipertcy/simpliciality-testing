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

from os.path import realpath, dirname, join
from setuptools import setup, find_packages
from simplicial_test import __version__

PROJECT_ROOT = dirname(realpath(__file__))

REQUIREMENTS_FILE = join(PROJECT_ROOT, 'requirements.txt')

with open(REQUIREMENTS_FILE) as f:
    install_reqs = f.read().splitlines()

install_reqs.append('setuptools')


def get_long_description():
    return ""


setup(
    name='simplicial-test',
    version=__version__,
    description='Python module to realize simplicial joint degree sequences',
    long_description=get_long_description(),
    long_description_content_type='text/x-rst',
    author='Tzu-Chi Yen',
    author_email='tzuchi@netscied.tw',
    license='LGPLv3',
    packages=find_packages(),
    package_data={
        'simplicial_test': [
            'COPYING',
            'COPYING.LESS',
            'README.rst',
            'requirements.txt'
        ]},
    include_package_data=True,
    install_requires=install_reqs,
    extras_require={
        'testing': ['pytest'],
        'docs': ['mkdocs']
    },
    url='https://github.com/junipertcy/simplicial-test',
    platforms='any',
    classifiers=[
        'Development Status :: 5 - Production/Stable',
        'Intended Audience :: Science/Research',
        'License :: OSI Approved :: GNU Lesser General Public License v3 (LGPLv3)',
        'Natural Language :: English',
        'Operating System :: Microsoft :: Windows',
        'Operating System :: POSIX',
        'Operating System :: Unix',
        'Operating System :: MacOS',
        'Programming Language :: Python :: 3.7',
        'Programming Language :: Python :: 3.8',
        'Programming Language :: Python :: 3.9',
        'Programming Language :: Python :: 3 :: Only',
        'Topic :: Scientific/Engineering',
        'Topic :: Software Development :: Libraries',
    ],
    keywords=[
        'simplicial-test',
        'simplicial complex',
        'degree sequence',
        'simplex',
        'topological data analysis',
        'topology',
        'recursive algorithm',
        'fixed-parameter tractable',
        'realizability',
    ],
)
