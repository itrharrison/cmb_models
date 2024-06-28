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

for H0_plot in H0_list:

    # set up a new set of parameters for CAMB
    pars = camb.set_params(H0=H0_plot, ombh2=0.02233, omch2=0.1198, mnu=0.06, omk=0, tau=0.06,  
                           As=2.1e-9, ns=0.965, halofit_version='mead', lmax=3000)

    # perform the camb calculation
    results = camb.get_results(pars)

    # extract dictionary of CAMB power spectra
    powers =results.get_cmb_power_spectra(pars, CMB_unit='muK')
    all_cls = powers['total']

    ells = np.arange(all_cls.shape[0])
    cl_tt = all_cls[:,0]

    plt.plot(ells, cl_tt, label=f'$H_0 = {H0_plot:.1f}$')

cmb_data = Table.read('data/COM_PowerSpect_CMB_R2.02.fits', hdu=7)
cmb_data_lowell = Table.read('data/COM_PowerSpect_CMB_R2.02.fits', hdu=1)
plt.errorbar(cmb_data['ELL'], cmb_data['D_ELL'], yerr=cmb_data['ERR'], fmt='o', markersize=2, capsize=3, color='k', label='\emph{Planck} 2018 data')
plt.errorbar(cmb_data_lowell['ELL'], cmb_data_lowell['D_ELL'], yerr=[cmb_data_lowell['ERRDOWN'], cmb_data_lowell['ERRUP']], fmt='o', markersize=2, capsize=3, label=None, color='k')

plt.legend(fontsize='small')
plt.xscale('log')
plt.xlabel('$l$')
plt.ylabel('$\mathcal{D}_l$')
plt.xlim([2,2000])

plt.savefig('./plots/cmb_example.png', dpi=300, bbox_inches='tight')