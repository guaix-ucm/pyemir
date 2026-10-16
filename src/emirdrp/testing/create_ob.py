import numina.core

from emirdrp.testing.create_base import create_image0
from emirdrp.testing.create_wcs import create_wcs_4829
from emirdrp.testing.create_scenes import (
    create_scene_noise,
    create_scene_4829,
    create_scene_4847,
    create_scene_23,
)


def exposure_keys(idx, exptime=100.0, starttime=0.0):
    """TSUTC1 and TSUTC2 of the image idx of a sequence of exposures without gaps"""
    t1 = starttime + idx * exptime
    return {"TSUTC1": t1, "TSUTC2": t1 + exptime}


def make_frame(scene, hdr=None, keys=None):
    """A frame with the image of scene, with the header hdr and the keywords in keys"""
    return numina.core.DataFrame(frame=create_image0(scene, hdr=hdr, keys=keys))


def make_obresult(frames, mode="UNKNOWN"):
    """An observation result of EMIR with frames"""
    obsresult = numina.core.ObservationResult(instrument="EMIR", mode=mode)
    obsresult.frames = list(frames)
    return obsresult


def create_ob_abba(value, exptime=100.0, starttime=0.0):
    pos = [20, 30, 30, 20]
    size = (100, 100)
    frames = []
    for i in range(4):
        keys = exposure_keys(i, exptime, starttime)
        keys["NUM-SK"] = 1
        scene = create_scene_noise(size, background=value, std=100, pos=pos[i], peak=value)
        frames.append(make_frame(scene, keys=keys))
    return make_obresult(frames, mode="LS_ABBA")


def create_ob_1(value, nimages, crpix, nstare=3, exptime=100.0, starttime=0.0):
    size = (100, 100)
    frames = []
    for i in range(nimages):
        keys = exposure_keys(i, exptime, starttime)
        keys.update({"NUM-NCOM": nstare, "NUM-SK": 1})
        wcs_scene = create_wcs_4829(crpix[i])
        scene = create_scene_4829(size, wcs_scene, background=value, std=100.0)
        frames.append(make_frame(scene, hdr=wcs_scene.to_header(), keys=keys))
    return make_obresult(frames)


def create_ob_2(value, nimages, exptime=100.0, starttime=0.0):
    frames = []
    for i in range(nimages):
        keys = exposure_keys(i, exptime, starttime)
        keys["DATE-OBS"] = "2019-04-12T03:01:05.678"
        frames.append(make_frame(create_scene_4847(val=value), keys=keys))
    return make_obresult(frames)


def create_ob_internal(nimages, exptime, starttime, pos):
    size = (10, 10)
    frames = []
    for j in range(nimages):
        keys = exposure_keys(j, exptime, starttime)
        frames.append(make_frame(create_scene_23(size, peak=2, pos=pos[j]), keys=keys))
    return make_obresult(frames)
