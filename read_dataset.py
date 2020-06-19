import numpy as np
import h5py

snapshot = 'galaxy1'
bigsnapshot = '/media/hunter/0AE21B3CE21B2C07/IC_DATA/Simulations/RefL0012N0188/snapshot_028_z000p000/snap_028_z000p000.0.hdf5'

def read_dataset(itype, att, nfiles=1):
    """ Read a selected dataset, itype is the PartType and att is the attribute name. """
    # Output array.
    data = []
    f = h5py.File(snapshot, 'r')
    tmp = f['PartType%i/%s'%(itype, att)][...]
    data.append(tmp)

    # Get expansion factor and Hubble parameter from the header.
    a       = f['Header'].attrs.get('Time')
    h       = f['Header'].attrs.get('HubbleParam')
    f.close()

    cgs_in, aexp_in, hexp_in = np.loadtxt('conversion_factors_%s'%snapshot, unpack=True, delimiter=',', usecols=(1,2,3))
    if att == 'Coordinates':
        i = 0
    if att == 'Velocity':
        i = 1
    if att == 'Mass':
        i = 2

    if itype == 0:
        cgs = cgs_in[0+i]
        aexp = aexp_in[0+i]
        hexp = hexp_in[0+i]
    if itype == 1:
        cgs = cgs_in[3+i]
        aexp = aexp_in[3+i]
        hexp = hexp_in[3+i]
    if itype == 2:
        cgs = cgs_in[6+i]
        aexp = aexp_in[6+i]
        hexp = hexp_in[6+i]
    if itype == 3:
        cgs = cgs_in[9+i]
        aexp = aexp_in[9+i]
        hexp = hexp_in[9+i]
    if itype == 4:
        cgs = cgs_in[12+i]
        aexp = aexp_in[12+i]
        hexp = hexp_in[12+i]

    print("Extracted CGS: %.2e"%cgs)
    print("Extracted AEXP: %.2e"%aexp)
    print("Extracted HEXP: %.2e"%hexp)

    # # Get conversion factors.
    # f2 = h5py.File(bigsnapshot, 'r')
    # cgs     = f2['PartType%i/%s'%(itype, att)].attrs.get('CGSConversionFactor')
    # aexp    = f2['PartType%i/%s'%(itype, att)].attrs.get('aexp-scale-exponent')
    # hexp    = f2['PartType%i/%s'%(itype, att)].attrs.get('h-scale-exponent')
    # f2.close()
    #
    # print("Header CGS: %.2e"%cgs)
    # print("Header AEXP: %.2e"%aexp)
    # print("Header HEXP: %.2e"%hexp)

    # Combine to a single array.
    if len(tmp.shape) > 1:
        data = np.vstack(data)
    else:
        data = np.concatenate(data)

    # Convert to physical.
    if data.dtype != np.int32 and data.dtype != np.int64:
        data = np.multiply(data, cgs * a**aexp * h**hexp, dtype='f8')

    return data

#def read_dataset_dm_mass():
    #""" Special case for the mass of dark matter particles. """
    #f           = h5py.File(snapshot, 'r')
    #h           = f['Header'].attrs.get('HubbleParam')
    #a           = f['Header'].attrs.get('Time')
    #dm_mass     = f['Header'].attrs.get('MassTable')[1]
    #n_particles = f['Header'].attrs.get('NumPart_Total')[1]

    ## Create an array of length n_particles each set to dm_mass.
    #m = np.ones(n_particles, dtype='f8') * dm_mass

    ## Use the conversion factors from the mass entry in the gas particles.
    #cgs  = f['PartType0/Mass'].attrs.get('CGSConversionFactor')
    #aexp = f['PartType0/Mass'].attrs.get('aexp-scale-exponent')
    #hexp = f['PartType0/Mass'].attrs.get('h-scale-exponent')
    #f.close()

    ## Convert to physical.
    #m = np.multiply(m, cgs * a**aexp * h**hexp, dtype='f8')

    #return m

##halo mass
#dm_mass = read_dataset_dm_mass()

#halo positions
pos = read_dataset(0,'Mass')
x = pos[:,0]
y = pos[:,1]
z = pos[:,2]

# #halo velocities
# vel = read_dataset(1,'Velocity')
# vx = vel[:,0]
# vy = vel[:,1]
# vz = vel[:,2]


#gas mass
#mass = read_dataset(0,'Mass')
