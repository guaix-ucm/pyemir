import emirdrp.loader
import numina.instrument.assembly as asb
import numina.instrument.generic


def test_assembly1():
    drp = emirdrp.loader.load_drp()
    pkg_paths = [drp.profiles]
    com_store = asb.load_paths_store(pkg_paths, [])
    emir_ins = asb.assembly_instrument(com_store, "443fc0d1-e09a-48cc-a0fd-02be6f399da2", "2023-10-01T12:00:00",
                                       by_key="uuid")
    assert isinstance(emir_ins, numina.instrument.generic.InstrumentGeneric)


def test_assembly2():
    import emirdrp.loader
    import numina.instrument.assembly as asb

    drp = emirdrp.loader.load_drp()
    pkg_paths = [drp.profiles]
    com_store = asb.load_paths_store(pkg_paths, [])
    emir_ins = asb.assembly_instrument(com_store, "EMIR", "2023-10-01T12:00:00", by_key="name")
    assert isinstance(emir_ins, numina.instrument.generic.InstrumentGeneric)
