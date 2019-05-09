
import pytest

from numina.user.baserun import run_reduce

# Variable for datamanager_remote
TEST_SET_FILE = "pyemir/pyemir-image-tutorial-v2.tar.gz"


@pytest.mark.result_compare
@pytest.mark.remote_data
@pytest.mark.parametrize("obsid", range(111, 118))
def test_block1(datamanager_remote, obsid):
    task = run_reduce(datamanager_remote, obsid)
    return task.result


@pytest.mark.result_compare
@pytest.mark.remote_data
@pytest.mark.parametrize("obsid", range(121, 128))
def test_block2(datamanager_remote, obsid):
    task = run_reduce(datamanager_remote, obsid)
    return task.result


@pytest.mark.result_compare
@pytest.mark.remote_data
@pytest.mark.parametrize("obsid", [20101])
def test_block3(datamanager_remote, obsid):
    task = run_reduce(datamanager_remote, obsid)
    return task.result


@pytest.mark.remote_data
def test_block4(datamanager_remote):
    from numina.util.context import working_directory

    with working_directory(datamanager_remote.basedir):
        with open('control_dump.yaml', 'w') as fp:
            datamanager_remote.backend.dump(fp)
        assert True
