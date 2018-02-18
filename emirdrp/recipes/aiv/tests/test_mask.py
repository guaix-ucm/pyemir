#
# Copyright 2015 Universidad Complutense de Madrid
#
# This file is part of PyEmir
#
# PyEmir is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# PyEmir is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with PyEmir.  If not, see <http://www.gnu.org/licenses/>.
#

"""Test the AIV pinhole mask recipe"""

import pytest

from numina.user.cli import main
from emirdrp.tests.runrecipe import run_recipe

BASE_URL = 'http://guaix.fis.ucm.es/~spr/emir_test/'

@pytest.mark.remote
def test_mode_TEST6_set0(numinatpldir):

    run_recipe()


@pytest.mark.remote
def test_mode_TEST7_set0(numinatpldir):

    run_recipe()


@pytest.mark.remote
def test_mode_TEST8_set0(numinatpldir):

    run_recipe()
