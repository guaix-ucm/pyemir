
import os
import tarfile

import pytest

from numina.tests.testcache import download_cache
from numina.user.cli import main

@pytest.mark.skip
def test_1(tmpdir):

    tarname = "pyemir/pyemir-image-tutorial-v1.tar.gz"
    base = "https://guaix.fis.ucm.es/data/"
    url = base + tarname
    downloaded = download_cache(url)
    tmpdir.chdir()

    # Uncompress
    with tarfile.open(downloaded.name, mode="r:gz") as tar:
        tar.extractall()

    os.remove(downloaded.name)
    print(tmpdir)

    main(['run', 'obs1.yaml', '-r', 'control.yaml'])
    assert False