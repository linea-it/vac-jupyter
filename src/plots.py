import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import matplotlib.cm as cm
from lib.skymapper import skymapper as skm
import healpy as hp
from astropy.io import fits

from src import db


cmap = mpl.cm.get_cmap("inferno_r")
cmap.set_under('dimgray')
    
    
def plot_map(data, nside=4096, label="titulo"): # Julia e Hillysson: added label
    ra = list()
    dec = list()
    signal = list()
    for d in data:
        ra.append(d[0])
        dec.append(d[1])
        signal.append(d[2])
    ra = np.array(ra)
    ra[ra > 180] = ra[ra > 180] - 360
    
    plot_HPsignal(ra, np.array(dec), np.array(signal), nside, True, label)
    
    
def plot_HPsignal(ra, dec, signal, nside, nest, title, border=1., norm=None):
    plt.clf()
    # from matplotlib import rc
    # rc('font',**{'family':'sans-serif','sans-serif':['DejaVu'], 'size': 8})
    # rc('font',**{'size': 8})
    # rc('text', usetex=True)
    ramin = np.min(ra)
    ramax = np.max(ra)
    decmin = np.min(dec)
    decmax = np.max(dec)
    aspect = (decmax - decmin) / ((ramax - ramin) * np.cos (np.deg2rad((decmin + decmax) / 2.)))
    ipix = hp.ang2pix(nside, (90-dec)/180*np.pi, ra/180*np.pi, nest=True)
    HPX = np.bincount(ipix, minlength=hp.nside2npix(nside), weights=signal)
    ext = [ramax + border, ramin - border, decmin - border, decmax + border]
    mask = np.zeros(hp.nside2npix(nside), dtype=np.bool)
    mask[HPX == 0.] = 1
    m = hp.ma(HPX)
    m.mask = mask
    d_op = hp.visufunc.cartview(map=m, lonra=[ramin - border, ramax + border],
                                latra=[decmin - border, decmax + border], flip='astro', format='%.3g', cmap=cmap,
                                nest=True, norm=norm, return_projected_map=True)
    fig = plt.figure(figsize = [6.4, 6.4 * aspect])
    ax1 = fig.add_subplot(111)
    ax1.set_facecolor('dimgray')
    plt.imshow(d_op, extent=ext, aspect='auto', origin='lower', interpolation=None, vmin=np.min(m),
               vmax=np.max(m), cmap=cmap)
    # plt.title(r'$\mathrm{%s}$' % title.replace('_', '\ '))
    plt.xlabel(r'R.A. (degrees)')
    plt.ylabel(r'Dec. (degrees)')
    plt.ylim([decmin - border, decmax + border])
    plt.xlim([ramax + border, ramin - border])
    plt.grid(c='grey', lw=0.25, alpha=0.5)
    cbaxes = fig.add_axes([0.90, 0.11, 0.015, 0.77])
    cb = plt.colorbar(cax=cbaxes, label=r'Map Value')
    # plt.savefig('%s_signal.jpg' % title, dpi=300, bbox_inches='tight')
    plt.show()
    # plt.close()
    return


def plot_catalog(nside, data, label="titulo"):
    ra = list()
    dec = list()
    for d in data:
        ra.append(d[0])
        dec.append(d[1])
    ra = np.array(ra)
    ra[ra > 180] = ra[ra > 180] - 360
    
    plot_HPcart(ra, np.array(dec), nside, label)
    
    
def plot_HPcart(ra, dec, nside, title, norm=None, border=1.):
    plt.clf()
    # from matplotlib import rc
    # rc('font',**{'family':'sans-serif','sans-serif':['DejaVu'], 'size': 8})
    # rc('font',**{'size': 8})
    # rc('text', usetex=True)
    ramin = np.min(ra)
    ramax = np.max(ra)
    decmin = np.min(dec)
    decmax = np.max(dec)
    aspect = (decmax - decmin) / ((ramax - ramin) * np.cos (np.deg2rad((decmin + decmax) / 2.)))

    ipix = hp.ang2pix(nside, (90-dec)/180*np.pi, ra/180*np.pi, nest=True)
    HPX = np.bincount(ipix, minlength=hp.nside2npix(nside))
    pixarea = hp.pixelfunc.nside2pixarea(nside, degrees=True)
    ext = [ramax + border, ramin - border, decmin - border, decmax + border]
    HPX = HPX / (3600. * pixarea)
    mask = np.zeros(hp.nside2npix(nside), dtype=np.bool)
    mask[HPX == 0.] = 1
    m = hp.ma(HPX)
    m.mask = mask
    d_op = hp.visufunc.cartview(map=m, lonra=[ramin - border, ramax + border],
                                latra=[decmin - border, decmax + border], flip='astro', format='%.3g', cmap=cmap,
                                nest=True, norm=norm, return_projected_map=True)
    fig = plt.figure(figsize = [6.4, 6.4 * aspect])
    ax1 = fig.add_subplot(111)
    ax1.set_facecolor('dimgray')
    
    if norm==None:
        plt.imshow(d_op, extent=ext, aspect='auto', origin='lower', interpolation=None, vmin=np.min(m),
                   vmax=np.max(m), cmap=cmap)
    else:
        plt.imshow(d_op, extent=ext, aspect='auto', origin='lower', interpolation=None, vmin=np.log10(np.min(m)),
                   vmax=np.log10(np.max(m)), cmap=cmap)
    # plt.title('$\mathrm{%s}$' % title.replace('_', '\ '))
    plt.xlabel('ra')
    plt.ylabel('dec')
    plt.ylim([decmin - border, decmax + border])
    plt.xlim([ramax + border, ramin - border])
    plt.grid(c='grey', lw=0.25, alpha=0.5)
    cbaxes = fig.add_axes([0.90, 0.11, 0.015, 0.77])
    if norm==None:
        cb = plt.colorbar(cax=cbaxes, label=r'$\mathrm{density\ (objects\ per\ sq.\ arcmin)}$')
        plt.savefig('linear.png', bbox_inches='tight')
    else:
        cb = plt.colorbar(cax=cbaxes, label=r'$\mathrm{log(density\ [objects\ per\ sq.\ arcmin])}$')
        plt.savefig('%slogscale.png' % title, bbox_inches='tight')
    plt.close()
    return
                  
                  
def plot_catalog(nside, data):
    ra = list()
    dec = list()
    for d in data:
        ra.append(d[0])
        dec.append(d[1])

    fig, ax, proj = skm.plotDensity(np.asarray(ra),
                                    np.asarray(dec),
                                    nside=nside)
    plt.tight_layout()
    # fig.savefig(save_as + '_surface.png')


def plot_signal_area(signals, NSIDE=4096, xlabel='', ylabel=''):   # Julia e Hillysson: added label  
    sig_vals = sorted(list(set(signals)))
    nsig = np.array([signals[signals >= s].size for s in sig_vals])

    print(signals.size)
    area_pix = 4 * np.pi * (180. / np.pi) ** 2 / (12 * 4096 ** 2)

    fig, axes = plt.subplots()
    axes.plot(sig_vals, nsig)
    axes.set_ylabel(ylabel)
    axes.set_xlabel(xlabel)

    ax2 = axes.twinx()
    ax2.plot(sig_vals, 100 * nsig / float(signals.size))
    ylim = ax2.get_ylim()
    for x in range(0, 101, 20):
        ax2.axhline(x, ls='--', color='.5', lw=.8)
    ax2.set_ylim(ylim)
    ax2.set_ylabel('$Area (\%)' + '\; [%.2f\, deg^2]$' % (area_pix * signals.size))
    plt.tight_layout()
    

