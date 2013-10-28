#
# Copyright 2008-2013 Universidad Complutense de Madrid
# 
# This file is part of PyEmir
# 
# PyEmir is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
# 
# PyEmir is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# 
# You should have received a copy of the GNU General Public License
# along with PyEmir.  If not, see <http://www.gnu.org/licenses/>.
# 

import itertools as ito

from numina.extraiter import braid


def _channel_gen1(beg, end, step):
    return ito.imap(lambda x: (x, x + step), xrange(beg, end, step))

def _channel_gen2(beg, end, step):
    return ito.imap(lambda x: (x - step, x), xrange(beg, end, -step))

def _ch1():
    return ito.izip(ito.repeat(slice(1024, 2048)), ito.starmap(slice, _channel_gen2(1024, 0, 128)))

def _ch2():  
    return ito.izip(ito.starmap(slice, _channel_gen2(1024, 0, 128)), ito.repeat(slice(0, 1024)))

def _ch3():
    return ito.izip(ito.repeat(slice(0, 1024)), ito.starmap(slice, _channel_gen1(1024, 2048, 128)))

def _ch4():
    return ito.izip(ito.starmap(slice, _channel_gen1(1024, 2048, 128)), ito.repeat(slice(1024, 2048)))

# Channels are listed per quadrant and then in fast readout order
CHANNELS = list(ito.chain(_ch1(), _ch2(), _ch3(), _ch4()))
# Channels as they are populated during reconstruction
CHANNELS_2 = list(ito.chain(_ch3(), _ch4(), _ch1(), _ch2()))
# Channels as listed in Carlos Gonzalez Ph. D. Thesis
CHANNELS_3 = list(ito.chain(reversed(list(_ch2())), 
                            reversed(list(_ch3())),
                            reversed(list(_ch4())),
                            reversed(list(_ch1()))
                            ))

# Channels in read out order
CHANNELS_READOUT = list(braid(_ch1(), _ch2(), _ch3(), _ch4()))
# Channels as they are populated during reconstruction
CHANNELS_READOUT_2 = list(braid(_ch3(), _ch4(), _ch1(), _ch2()))

# Quadrants are listed starting at left-top and counter-clockwise then
QUADRANTS = [(slice(1024, 2048), slice(0, 1024)),
             (slice(0, 1024), slice(0, 1024)),
             (slice(0, 1024), slice(1024, 2048)),
             (slice(1024, 2048), slice(1024, 2048))
             ]

FULL = CHANNELS_2

# FIXME: this is a hack to convert channel name to a structure
def convert_name_to_channels(conf):
    chname = conf.configuration['detector']['channels']
    try:
        allcha = globals()[chname]
        conf.configuration['detector']['channels'] = allcha
        return conf
    except KeyError:
        _logger.warning("incorrect %r name %s", 'channels', chname)
        return None
