import numpy as np
from lib.skymapper import skymapper as skm


def plot_catalog(nside, data):
    ra = list()
    dec = list()
    for d in data:
        ra.append(d[0])
        dec.append(d[1])

    fig, ax, proj = skm.plotDensity(np.asarray(ra),
                                    np.asarray(dec),
                                    nside=nside)
    # fig.savefig(save_as + '_surface.png')


def plot_map(data):
    ra = list()
    dec = list()
    signal = list()
    for d in data:
        ra.append(d[0])
        dec.append(d[1])
        signal.append(d[2])

    fig, ax, proj = skm.plotMap(np.asarray(ra),
                                np.asarray(dec),
                                np.asarray(signal))
    # fig.savefig(save_as + '.png')
