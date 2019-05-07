
from numina.tests.plugins import *

import os
import tarfile

import pytest

from numina.tests.testcache import download_cache
from numina.user.helpers import create_datamanager


@pytest.fixture(scope='module')
def datamanager_remote(tmp_path_factory, request):
    """Return a DataManager object create from a remote dataset"""
    from numina.util.context import working_directory

    req_base_default = "https://guaix.fis.ucm.es/data/"
    req_base = getattr(request.module, 'TEST_SET_HOST', req_base_default)
    req_tarname = getattr(request.module, 'TEST_SET_FILE')
    req_datadir = getattr(request.module, 'TEST_SET_DATADIR', 'data')
    req_control = getattr(request.module, 'TEST_SET_CONTROL', "control_v2.yaml")

    basedir = tmp_path_factory.mktemp('manager')

    datadir = basedir / req_datadir  # pathlib syntax
    reqfile  = basedir / req_control

    if req_tarname is None:
        raise ValueError('Undefined TEST_SET_FILE')

    url = req_base + req_tarname

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
