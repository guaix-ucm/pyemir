#
# Copyright 2008-2015 Universidad Complutense de Madrid
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

def braid(*iterables):
    """Return the elements of each iterator in turn until some is exhausted.
    This function is similar to the roundrobin example
    in itertools documentation.
    >>> a = iter([1,2,3,4])
    >>> b = iter(['a', 'b'])
    >>> c = iter([1,1,1,1,'a', 'c'])
    >>> d = iter([1,1,1,1,1,1])
    >>> list(braid(a, b, c, d))
    [1, 'a', 1, 1, 2, 'b', 1, 1]
    """

    for itbl in zip(*iterables):
        for it in itbl:
            yield it


_P1 = slice(0, 1024)
_P2 = slice(1024, 2048)
_O1 = [slice(i*128,(i+1)*128) for i in range(8)]
_O2 = [slice(1024 + (7-i)*128, 1024+(8-i)*128) for i in range(8)]


_CH1 = [(_P2, s) for s in _O2]
_CH2 = [(s, _P1) for s in _O2]
_CH3 = [(_P1, s) for s in _O1]
_CH4 = [(s, _P2) for s in _O1]

# Channels are listed per quadrant and then in fast readout order
CHANNELS = [chan for chans in [_CH1, _CH2, _CH3, _CH4] for chan in chans]
CHANNELS_1 = CHANNELS
# Channels as they are populated during reconstruction
CHANNELS_2 = [chan for chans in [_CH3, _CH4, _CH1, _CH2] for chan in chans]
# Channels as listed in Carlos Gonzalez Ph. D. Thesis
CHANNELS_3 = [chan for chans in [_CH2, _CH3, _CH4, _CH1] for chan in reversed(chans)]
# Channels in readout order
CHANNELS_READOUT = list(braid(_CH3,_CH4, _CH1, _CH2))

RCHANNELS_1 = [chan for chans in [_CH3, _CH1, _CH2, _CH4] for chan in chans]
FULL = RCHANNELS_1


# Quadrants are listed starting at left-top and counter-clockwise then
QUADRANTS = [(_P2, _P1), (_P1, _P1), (_P1, _P2), (_P2, _P2)]


# FIXME: this is a hack to convert channel name to a structure
def convert_name_to_channels(conf):
    chname = conf.configuration['detector']['channels']
    allcha = globals()[chname]
    conf.configuration['detector']['channels'] = allcha
    return conf
