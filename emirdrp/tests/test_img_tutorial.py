
import os
import tarfile

import pytest

from numina.tests.testcache import download_cache
from numina.user.baserun import run_job
from numina.user.helpers import create_datamanager


@pytest.fixture(scope='module')
def datamanagerfl(tmp_path_factory):
    from numina.util.context import working_directory

    basedir = tmp_path_factory.mktemp('manager')
    datadir = basedir / 'data'  # pathlib syntax
    reqfile  = basedir / 'control_v2.yaml'

    tarname = "pyemir/pyemir-image-tutorial-v2.tar.gz"
    base = "https://guaix.fis.ucm.es/data/"
    url = base + tarname

    # Download everything
    with working_directory(basedir):

        downloaded = download_cache(url)

        # Uncompress
        with tarfile.open(downloaded.name, mode="r:gz") as tar:
            tar.extractall()

        os.remove(downloaded.name)

    # Insert OBS in the control file....
    dm = create_datamanager(reqfile, basedir, datadir)

    # This is not really needed...
    # If everything is in the file already
    # with working_directory(basedir):
    #     obsresults = ['obs_ids.yaml']
    #     sessions, loaded_obs = load_observations(obsresults, is_session=False)
    #    dm.backend.add_obs(loaded_obs)

    return dm


@pytest.mark.remote_data
@pytest.mark.parametrize("obsid", range(111, 118))
def test_block1(datamanagerfl, obsid):
    run_job(datamanagerfl, obsid)
    assert True


@pytest.mark.remote_data
@pytest.mark.parametrize("obsid", range(121, 128))
def test_block2(datamanagerfl, obsid):
    print(datamanagerfl)
    run_job(datamanagerfl, obsid)
    assert True


@pytest.mark.remote_data
@pytest.mark.parametrize("obsid", [10101,20101])
def test_block3(datamanagerfl, obsid):
    run_job(datamanagerfl, obsid)
    assert True


@pytest.mark.remote_data
def test_block4(datamanagerfl):
    from numina.util.context import working_directory

    with working_directory(datamanagerfl.basedir):
        with open('control_dump.yaml', 'w') as fp:
            datamanagerfl.backend.dump(fp)
        assert True
