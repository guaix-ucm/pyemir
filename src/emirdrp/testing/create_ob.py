import numina.core

from emirdrp.testing.create_base import create_image0
from emirdrp.testing.create_wcs import create_wcs_4829
from emirdrp.testing.create_scenes import (
    create_scene_noise,
    create_scene_4829,
    create_scene_4847,
    create_scene_23,
)


def create_ob_abba(value, exptime=100.0, starttime=0.0):
    pos = [20, 30, 30, 20]
    size = (100, 100)
    frames = []
    nimages = 4
    for i in range(nimages):
        base = starttime
        off1 = i * exptime
        off2 = off1 + exptime
        t1 = base + off1
        t2 = base + off2
        keys = {"TSUTC1": t1, "TSUTC2": t2, "NUM-SK": 1}
        scene = create_scene_noise(
            size, background=value, std=100, pos=pos[i], peak=value
        )
        frame = numina.core.DataFrame(frame=create_image0(scene, keys=keys))
        frames.append(frame)

    obsresult = numina.core.ObservationResult()
    obsresult.frames = frames
    obsresult.mode = "LS_ABBA"
    return obsresult


def create_ob_1(value, nimages, crpix, nstare=3, exptime=100.0, starttime=0.0):
    size = (100, 100)
    frames = []

    for i in range(nimages):
        base = starttime
        off1 = i * exptime
        off2 = off1 + exptime
        t1 = base + off1
        t2 = base + off2
        keys = {"TSUTC1": t1, "TSUTC2": t2, "NUM-NCOM": nstare, "NUM-SK": 1}
        wcs_scene = create_wcs_4829(crpix[i])
        scene = create_scene_4829(size, wcs_scene, background=value, std=100.0)
        img = create_image0(scene, hdr=wcs_scene.to_header(), keys=keys)
        frame = numina.core.DataFrame(frame=img)
        frames.append(frame)

    obs_result = numina.core.ObservationResult()
    obs_result.frames = frames
    return obs_result


def create_ob_2(value, nimages, exptime=100.0, starttime=0.0):

    frames = []
    for i in range(nimages):
        base = starttime
        off1 = i * exptime
        off2 = off1 + exptime
        t1 = base + off1
        t2 = base + off2
        keys = {"DATE-OBS": "2019-04-12T03:01:05.678", "TSUTC1": t1, "TSUTC2": t2}

        scene = create_scene_4847(val=value)
        frame = numina.core.DataFrame(frame=create_image0(scene=scene, keys=keys))
        frames.append(frame)

    obs_result = numina.core.ObservationResult()
    obs_result.frames = frames
    return obs_result


def create_ob_internal(nimages, exptime, starttime, pos):
    frames = []
    size = (10, 10)
    for j in range(nimages):
        base = starttime
        off1 = j * exptime
        off2 = off1 + exptime
        t1 = base + off1
        t2 = base + off2
        keys = {"TSUTC1": t1, "TSUTC2": t2}
        scene = create_scene_23(size, peak=2, pos=pos[j])
        frame = numina.core.DataFrame(frame=create_image0(scene, keys=keys))
        frames.append(frame)

    obs_result = numina.core.ObservationResult()
    obs_result.frames = frames
    return obs_result
