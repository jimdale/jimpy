"""
Attempt to take Jim's sims, grid them with yt, and process them with radmc3d.
Doesn't work yet...


1. Need to write temperature, density, velocity, grid
2. Need to copy over dust opacity
3. Compute the line populations first: 
   (need to set lines_mode = 3 before doing this)
    radmc3d calcpop
4. Then compute the "image":
    radmc3d image iline 1 widthkms 10 linenlam 40 linelist nostar


"""
from yt.analysis_modules.radmc3d_export.api import RadMC3DWriter
from yt.utilities.physical_constants import kboltz
import yt

# this loads a specific particle data set
from load_runi import ds

#ds.covering_grid(1, [-1,-1,-1], [256,256,256]).write_out('covering_grid.hdf5')
#ds = yt.load('covering_grid.hdf5')

x_co = 1.0e-4
x_h2co = 1.0e-9
mu_h = yt.YTQuantity(2.34e-24, 'g')

ds.add_deposited_particle_field(('all', 'particle_temperature'), 'cic')
ds.add_field(('gas', 'temperature'), lambda x,y: y[('deposit', 'all_cic_temperature')], units='K')
ds.add_field(('gas', 'density'), lambda x,y: y[('deposit', 'all_cic')], units='g/cm**3')

def _NumberDensityCO(field, data):
    return (x_co/mu_h)*data[('deposit','all_cic')] # data['density']#
ds.add_field(("gas", "number_density_CO"), function=_NumberDensityCO, units="cm**-3")
def _NumberDensityH2(field, data):
    return (1./mu_h)*data[('deposit','all_cic')] # data['density']#
ds.add_field(("gas", "number_density_H2"), function=_NumberDensityH2, units="cm**-3")
def _NumberDensityH2CO(field, data):
    return (1./mu_h)*x_h2co*data[('deposit','all_cic')]*x_h2co # data['density']#
ds.add_field(("gas", "number_density_H2CO"), function=_NumberDensityH2CO, units="cm**-3")
dust_to_gas = 0.01
def _DustDensity(field, data):
    return dust_to_gas * data[('deposit','all_cic')] # data['density']#
ds.add_field(("gas", "dust_density"), function=_DustDensity, units="g/cm**3")

writer = RadMC3DWriter(ds)
# YTFieldNotFound: Could not find field '('all', 'density')' in ParticleData.


writer.write_amr_grid()
writer.write_line_file(("gas", "number_density_CO"), "numberdens_co.inp")
writer.write_line_file(("gas", "number_density_H2"), "numberdens_h2.inp")
writer.write_line_file(("gas", "number_density_H2CO"), "numberdens_h2co.inp")
#writer.write_dust_file(("gas", "temperature"), "gas_temperature.inp")
writer.write_dust_file(("gas", "dust_density"), "dust_density.inp")
velocity_fields = [('deposit',"all_cic_velocity_x"), ('deposit',"all_cic_velocity_y"),
                   ('deposit',"all_cic_velocity_z")]
writer.write_line_file(velocity_fields, "gas_velocity.inp")


import radmc3dPy
import os
import subprocess
import shutil

shutil.copy('/Users/adam/repos/radmc-3d/version_0.39/python/python_examples/datafiles/dustkappa_silicate.inp', '.')
shutil.copy('/Users/adam/repos/radmc-3d/version_0.39/python/python_examples/datafiles/molecule_co.inp', '.')
shutil.copy('/Users/adam/LAMDA/ph2co-h2.dat','molecule_h2co.inp')

params=dict(istar_sphere=0, itempdecoup=1, lines_mode=3, nphot=1000000,
            nphot_scat=30000, nphot_spec=100000, rto_style=3,
            scattering_mode=0, scattering_mode_max=1, tgas_eq_tdust=1,)

params_string = """
istar_sphere = {istar_sphere}
itempdecoup = {itempdecoup}
lines_mode = {lines_mode}
nphot = {nphot}
nphot_scat = {nphot_scat}
nphot_spec = {nphot_spec}
rto_style = {rto_style}
scattering_mode = {scattering_mode}
scattering_mode_max = {scattering_mode_max}
tgas_eq_tdust = {tgas_eq_tdust}
"""

with open('radmc3d.inp','w') as f:
    params['lines_mode'] = 50
    f.write(params_string.format(**params))

assert os.system('radmc3d calcpop writepop') == 0

with open('radmc3d.inp','w') as f:
    params['lines_mode'] = 3
    f.write(params_string.format(**params))

# compute the dust temperature
assert os.system('radmc3d mctherm') == 0

# iline: 1 = CO 1-0, 2 = CO 2-1, etc.
# widthkms = full width of output spectrum, divided by linenlam
# linenlam: number of wavelengtyh bins
# linelist
radmc3dPy.image.makeImage(iline=1, widthkms=5, linenlam=40, nostar=True, npix=500, incl=0, lambdarange=[4115.882695765051*(1-10/3e5),
                                                                                                        4115.882695765051*(1+10/3e5)],
                          nlam=20, sizeau=100000)
shutil.move('image.out', 'h2co_1-0_image.out')
radmc3dPy.image.makeImage(iline=3, widthkms=5, linenlam=40, nostar=True)
shutil.move('image.out', 'h2co_303-202_image.out')
radmc3dPy.image.makeImage(iline=13, widthkms=5, linenlam=40, nostar=True)
shutil.move('image.out', 'h2co_321-220_image.out')
radmc3dPy.image.makeImage(iline=2, widthkms=5, linenlam=40, nostar=True)
shutil.move('image.out', 'co_2-1_image.out')
radmc3dPy.image.makeImage(iline=1, widthkms=5, linenlam=40, nostar=True)
shutil.move('image.out', 'co_1-0_image.out')
#radmc3d image iline 1 widthkms 10 linenlam 40 linelist nostar
