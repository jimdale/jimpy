from fort_dump import fort_dump
import yt
import numpy as np
from yt.units import parsec, Msun, gram, centimeter, second

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

    # Only keep the particles that are "turned on", or something like that
    # (this is needed to remove the streaks from accreted particles that track the motion of the stars)
    keep = iphase==0

    ppx=x[keep]
    ppy=y[keep]
    ppz=z[keep]

    data = {'particle_position_x': ppx,
            'particle_position_y': ppy,
            'particle_position_z': ppz,
            'particle_velocity_x': vx[keep],
            'particle_velocity_y': vy[keep],
            'particle_velocity_z': vz[keep],
            'particle_mass': m[keep],
            'smoothing_length': h[keep],
            'particle_temperature': u[keep],}


    bbox = 1.1*np.array([[min(ppx), max(ppx)], [min(ppy), max(ppy)], [min(ppz), max(ppz)]])

    ds = yt.load_particles(data, length_unit=udisti*centimeter,
                           mass_unit=umassi*gram, n_ref=n_ref,
                           velocity_unit=udisti/utimei*centimeter/second,
                           bbox=bbox)

    return ds
