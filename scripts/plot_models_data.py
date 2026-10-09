import numpy as np

import camb
from camb import model, initialpower

from astropy.table import Table

from matplotlib import pyplot as plt
from matplotlib import rc

rc('text', usetex=True)
rc('font', family='serif')
rc('font', size=11)

plt.close('all')

# set the list of cosmologies we want to plot
H0_list = np.linspace(50, 90, 5)
H0_list = np.append(H0_list, [67.])

for H0_plot in H0_list:

    # set up a new set of parameters for CAMB
    pars = camb.set_params(H0=H0_plot, ombh2=0.02233, omch2=0.1198, mnu=0.06, omk=0, tau=0.06,  
                           As=2.1e-9, ns=0.965, halofit_version='mead', lmax=6000)

    # perform the camb calculation
    results = camb.get_results(pars)

    # extract dictionary of CAMB power spectra
    powers =results.get_cmb_power_spectra(pars, CMB_unit='muK')
    all_cls = powers['total']

    ells = np.arange(all_cls.shape[0])
    cl_tt = all_cls[:,0]

    if H0_plot == 67.:
        continue
    else:
        plt.plot(ells, cl_tt, label=f'$H_0 = {H0_plot:.1f}$')


planck_data = Table.read('data/COM_PowerSpect_CMB_R2.02.fits', hdu=7)
planck_data_lowell = Table.read('data/COM_PowerSpect_CMB_R2.02.fits', hdu=1)
plt.errorbar(planck_data['ELL'], planck_data['D_ELL'], yerr=planck_data['ERR'], fmt='o', markersize=2, capsize=3, color='k', label='\emph{Planck} PDR2 (2015) data')
plt.errorbar(planck_data_lowell['ELL'], planck_data_lowell['D_ELL'], yerr=[planck_data_lowell['ERRDOWN'], planck_data_lowell['ERRUP']], fmt='o', markersize=2, capsize=3, label=None, color='k')

plt.legend(fontsize='small')
plt.xscale('log')
plt.yscale('log')
plt.xlabel('$l$')
plt.ylabel('$\mathcal{D}_l$')
plt.xlim([2,3000])
plt.ylim([3.e1, 1.e4])

plt.savefig('./plots/cmb_example.png', dpi=300, bbox_inches='tight')

plt.clf()

plt.errorbar(planck_data['ELL'], planck_data['D_ELL'], yerr=planck_data['ERR'], fmt='o', markersize=2, capsize=3, color='k', label='\emph{Planck} PDR2 (2015) data')
plt.errorbar(planck_data_lowell['ELL'], planck_data_lowell['D_ELL'], yerr=[planck_data_lowell['ERRDOWN'], planck_data_lowell['ERRUP']], fmt='o', markersize=2, capsize=3, label=None, color='k')

plt.legend(fontsize='small')
plt.xscale('log')
plt.yscale('log')
plt.xlabel('$l$')
plt.ylabel('$\mathcal{D}_l$')
plt.xlim([2,6000])
plt.ylim([1.e1, 1.e4])

fsky = 0.4
mask = np.logical_and(ells >= 40, ells <= 6000)
TT_ells, TT_Nell = np.loadtxt("data/SO_LAT_Nell_T_atmv1_baseline_fsky0p4_ILC_CMB.txt", usecols = (0, 1), unpack = True)
TT_Nell = TT_Nell[np.logical_and(TT_ells >= 40, TT_ells <= 6000)] / (2.7255e6 ** 2.0)
nltt = np.sqrt(2.0 * (cl_tt[mask] + TT_Nell) ** 2.0 / (fsky * (2.0 * ells[mask] + 1.0)))

plt.errorbar(TT_ells[np.logical_and(TT_ells >= 40, TT_ells <= 6000)], cl_tt[mask], yerr=nltt, fmt='o', markersize=2, capsize=3, color='r', label='SO baseline', zorder=-1000)

plt.legend(fontsize='small')

plt.savefig('./plots/cmb_TT_planck_so.png', dpi=300, bbox_inches='tight')