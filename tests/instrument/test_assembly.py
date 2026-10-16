import numina.instrument.assembly as asb
import numina.instrument.generic


def test_assembly1(component_store):
    emir_ins = asb.assembly_instrument(
        component_store,
        "443fc0d1-e09a-48cc-a0fd-02be6f399da2",
        "2023-10-01T12:00:00",
        by_key="uuid",
    )
    assert isinstance(emir_ins, numina.instrument.generic.InstrumentGeneric)
    assert "dtu" in emir_ins.children


def test_assembly2(component_store):
    emir_ins = asb.assembly_instrument(component_store, "EMIR", "2023-10-01T12:00:00", by_key="name")
    assert isinstance(emir_ins, numina.instrument.generic.InstrumentGeneric)
