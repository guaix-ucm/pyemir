
import pytest

from numina.user.baserun import run_reduce

# Variable for datamanager_remote
TEST_SET_FILE = "pyemir/pyemir-image-tutorial-v2.tar.gz"


@pytest.mark.result_compare
@pytest.mark.remote_data
@pytest.mark.parametrize("obsid", range(111, 118))
def test_img_tut_block1(datamanager_remote, obsid):
    task = run_reduce(datamanager_remote, obsid)
    return task.result


@pytest.mark.result_compare
@pytest.mark.remote_data
@pytest.mark.parametrize("obsid", range(121, 128))
def test_test_img_tut_block2(datamanager_remote, obsid):
    task = run_reduce(datamanager_remote, obsid)
    return task.result


reduction_reqs1 = {
    'iterations': 1,
    'offsets': 'user_offsets.txt',
    'sky_images': 3
}


reduction_reqs2 = {
    'iterations': 1,
    'sky_images': 3
}

@pytest.mark.result_compare
@pytest.mark.remote_data
@pytest.mark.parametrize("obsid", [20101])
@pytest.mark.parametrize("requirements", [reduction_reqs1, reduction_reqs2])
def test_test_img_tut_block3(datamanager_remote, obsid, requirements):
    task = run_reduce(datamanager_remote, obsid, requirements=requirements)
    return task.result


@pytest.mark.remote_data
def test_test_img_tut_block4(datamanager_remote):
    from numina.util.context import working_directory

    with working_directory(datamanager_remote.basedir):
        with open('control_dump.yaml', 'w') as fp:
            datamanager_remote.backend.dump(fp)
        assert True
