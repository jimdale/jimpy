from fort_dump import fort_dump
import yt
import numpy as np
from yt.units import parsec, Msun, gram, centimeter, second, Kelvin
from yt.fields.particle_fields import add_volume_weighted_smoothed_field

from astropy import units as u
from astropy import constants

def yt_from_jim(fname, n_ref=8):
    """
    Load one of Jim's simulations into yt

    Parameters
    ----------
    n_ref : int
        See http://yt-project.org/docs/3.0/examining/loading_data.html#indexing-criteria
        It is the number of particles included per octree grid before refining that cell.
        Higher numbers -> less noisy, lower numbers -> better resolution
    """
    umassi,utimei,udisti,npart,gt,tkin,tgrav,tterm,escap,rho,poten,x,y,z,m,h,vx,vy,vz,u,iphase = fort_dump(fname)
    # poten = potential energy of cloud [junk]
    # escap = fraction of unbound mass (?)
    # gt = global time
    # tkin = total kinetic energy [code units]
    # tgrav = total gravitational energy [code units]
    # tterm = total t'ermal energy [code units]
    # u = specific internal energy [ergs/g]

    # Only keep the particles that are "turned on", or something like that
    # (this is needed to remove the streaks from accreted particles that track
    # the motion of the stars)
    keep = iphase==0

    ppx=x[keep]
    ppy=y[keep]
    ppz=z[keep]

    ergpergram = udisti**2 / utimei**2
    temperature = 2*ergpergram/constants.R.cgs.value*(2./3.)*u
    # Temperatures above 1000 are ionized, therefore have smaller mean mol weight
    temperature[temperature > 1000] /= 4.

    data = {'particle_position_x': ppx,
            'particle_position_y': ppy,
            'particle_position_z': ppz,
            'particle_velocity_x': vx[keep],
            'particle_velocity_y': vy[keep],
            'particle_velocity_z': vz[keep],
            'particle_mass': m[keep],
            'smoothing_length': h[keep],
            'particle_temperature': temperature[keep]*Kelvin,
           }


    bbox = 1.1*np.array([[min(ppx), max(ppx)],
                         [min(ppy), max(ppy)],
                         [min(ppz), max(ppz)]])

    ds = yt.load_particles(data,
                           length_unit=udisti*centimeter,
                           mass_unit=umassi*gram, n_ref=n_ref,
                           velocity_unit=udisti/utimei*centimeter/second,
                           time_unit=utimei*second,
                           sim_time=gt*utimei*second,
                           periodicity=(False,False,False),
                           bbox=bbox)

    def _temperature(field, data):
        ret = data[('all', "particle_temperature")] * Kelvin
        return ret.in_units('K')

    ds.add_field(('gas', "Temperature"), function=_temperature,
                 particle_type=True, units="K")

    num_neighbors = 64
    #fn = add_volume_weighted_smoothed_field('gas', "particle_position",
    #                                        "particle_mass",
    #                                        "smoothing_length", "density",
    #                                        "Temperature", ds, num_neighbors)

    #ds.alias(("gas", "temperature"), fn[0])
    

    return ds
