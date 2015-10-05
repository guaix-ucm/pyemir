import yaml
from numina.core import init_drp_system
from numina.core.pipeline import Instrument, Pipeline

def assert_valid_instrument(instrument):
    assert isinstance(instrument, Instrument)

    pipes = instrument.pipelines
    assert 'default' in pipes
    for k, v in pipes.items():
        assert k == v.name
        assert isinstance(v, Pipeline)

def test_recipes_are_defined():

    drp_to_test = """
    id: 1
    mode: bias
    instrument: EMIR
    images:
     - ThAr_LR-U.fits
    """

    loaded_obs = {}
    loaded_ids = []
    for doc in yaml.load_all(drp_to_test):
        loaded_ids.append(doc['id'])
        loaded_obs[doc['id']] = doc

    modes = ['fail', 'IMAGE_BIAS', 'IMAGE_DARK', 'IMAGE_FLAT', 'STARE_IMAGE', 'NODDED_BEAM_SWITCHED_IMAGE', 'DITHERED_IMAGE', 'MICRODITHERED_IMAGE', 'MOSAICED_IMAGE', 'gain_mode1', 'cosmetics', 'dark_current', 'simple_bias', 'TEST0', 'TEST1', 'TEST2', 'TEST3', 'TEST5', 'TEST6', 'IMAGE_SKY', 'TEST7', 'TEST8', 'ARC_CALIBRATION', 'TEST9', 'TEST9', 'FULL_DITHERED_IMAGE']

    m = init_drp_system(loaded_obs)
    for k, v in m.items():
        assert_valid_instrument(v)
        for m in v.modes:
            assert m.key in modes
            modes.remove(m.key)
            print modes

    assert len(modes)== 0

if __name__ == "__main__":
    test_recipes_are_defined()