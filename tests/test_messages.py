"""Details of the errors in their messages"""

import pytest

from emirdrp.processing.wavecal.set_wv_parameters import set_wv_parameters


def test_invalid_filter_and_grism(capsys):
    """A valid filter with a valid grism of another filter"""
    with pytest.raises(ValueError, match="filter_name='J', grism_name='H'"):
        set_wv_parameters("J", "H")
    assert capsys.readouterr().out == ""
