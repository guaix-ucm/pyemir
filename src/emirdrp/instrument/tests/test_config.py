
import pytest

import emirdrp.loader
import numina.instrument.assembly as asb
import numina.instrument.generic

from .dtuheader import DTU_HEADER_EXAMPLE


@pytest.fixture
def component_store():
    drp = emirdrp.loader.load_drp()
    pkg_paths = [drp.profiles]
    com_store = asb.load_paths_store(pkg_paths, [])
    return com_store


@pytest.fixture
def emir_ins_v2(component_store):
    emir_ins = asb.assembly_instrument(
        component_store,
        "443fc0d1-e09a-48cc-a0fd-02be6f399da2", "2023-10-01T12:00:00",
        by_key="uuid"
    )

    return emir_ins


@pytest.fixture
def emir_ins_v1(component_store):
    emir_ins = asb.assembly_instrument(
        component_store,
        "225fcaf2-7f6f-49cc-972a-70fd0aee8e96", "2022-10-01T12:00:00",
        by_key="uuid")

    return emir_ins


def test_emir2_detector_1(emir_ins_v2):

    emir_det_com = emir_ins_v2.children['detector']

    assert isinstance(emir_det_com, numina.instrument.generic.ComponentGeneric)

    assert emir_det_com.get_property('channels') == 'H2RG_FULL'
    assert emir_det_com.get_property('detector.channels') == 'H2RG_FULL'
    assert str(emir_det_com.uuid) == 'bc90ac77-2fd4-47c0-97ff-88a278a12785'
    assert emir_det_com.name == 'detector'
    assert emir_det_com.origin.name == 'detector'
    assert str(emir_det_com.origin.uuid) == 'bc90ac77-2fd4-47c0-97ff-88a278a12785'
    assert emir_det_com.origin.description == "H2RG detector configuration"
    assert str(emir_det_com.origin.date_start) == '2023-07-01 12:00:00'
    assert emir_det_com.origin.date_end is None

    # Configured only
    assert emir_det_com.is_configured is False


def test_emir2_dtu_1(emir_ins_v2):

    comp = emir_ins_v2.children['dtu']
    # Configured only
    assert comp.is_configured is False

    comp.configure_with_header(DTU_HEADER_EXAMPLE)
    assert comp.is_configured is True

    assert comp.get_property('ydtu_r') == -0.7186057740102463


def test_config_detector2(emir_ins_v1):
    emir_ins = emir_ins_v1

    emir_det_com = emir_ins.children['detector']

    assert isinstance(emir_det_com, numina.instrument.generic.ComponentGeneric)

    assert emir_det_com.get_property('channels') == 'FULL'
    assert emir_det_com.get_property('detector.channels') == 'FULL'

    # Configured only
    assert emir_det_com.is_configured is False
