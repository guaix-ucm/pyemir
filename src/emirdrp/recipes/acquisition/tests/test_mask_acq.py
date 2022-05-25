
import pytest

from numina.user.baserun import run_reduce

# Variable for datamanager_remote
TEST_SET_FILE = "pyemir/test/img_acq_set_20200702.tar.gz"


@pytest.mark.result_compare
@pytest.mark.remote_data
def test_img_tut_block1(datamanager_remote):
    task = run_reduce(datamanager_remote, "test1")
    return task.result


@pytest.mark.result_compare
@pytest.mark.remote_data
def test_test_img_tut_block2(datamanager_remote):
    task = run_reduce(datamanager_remote, "test2")
    return task.result


@pytest.mark.remote_data
def test_test_img_tut_block4(datamanager_remote):
    from numina.util.context import working_directory

    with working_directory(datamanager_remote.basedir):
        with open('control_dump.yaml', 'w') as fp:
            datamanager_remote.backend.dump(fp)
        assert True
