import numpy as np
import matplotlib.pyplot as plt
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


def plot_signal_area(data, NSIDE=4096):
    signals = np.array(data).T[0]
    sig_vals = sorted(list(set(signals)))
    nsig = np.array([signals[signals >= s].size for s in sig_vals])

    print(signals.size)
    area_pix = 4 * np.pi * (180. / np.pi) ** 2 / (12 * 4096 ** 2)

    fig, axes = plt.subplots()
    axes.plot(sig_vals, nsig)
    axes.set_ylabel('$N_{pix}$')
    axes.set_xlabel('signal>')

    ax2 = axes.twinx()
    ax2.plot(sig_vals, 100 * nsig / float(signals.size))
    ylim = ax2.get_ylim()
    for x in range(0, 101, 20):
        ax2.axhline(x, ls='--', color='.5', lw=.8)
    ax2.set_ylim(ylim)
    ax2.set_ylabel('$Area (\%)' + '\; [%.2f\, deg^2]$' % (area_pix * signals.size))
    fig.show()
