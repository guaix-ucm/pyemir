import pytest

from tests.runrecipe import run_recipe


BASE_URL = "https://guaix.fis.ucm.es/data/pyemir/test/"


@pytest.mark.remote_data
def _test_mode_TEST6_set0(numinatpldir):
    run_recipe()


@pytest.mark.remote_data
def _test_mode_TEST7_set0(numinatpldir):
    run_recipe()


@pytest.mark.remote_data
def _test_mode_TEST8_set0(numinatpldir):
    run_recipe()
