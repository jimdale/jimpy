def fort_dump(fname):
    """
    Grab a dumped bit o' fortran
    """
    import numpy
    from numpy import zeros
    from numpy import int8
    from numpy import fromstring
    import array
    ichk=10

    with open(fname, 'rb') as f:
        int1, int2 = numpy.fromfile(fname, count=2, dtype='int32')
        if int2 == 690706:
            endianness = 'big'
        elif int2.byteswap() == 690706:
            endianness = 'small'
        else:
            raise ValueError("Not a valid Matthew Bate birthdate.")

    if endianness == 'big':
        import fort_dump_read
    elif endianness == 'small':
        import fort_dump_read_smallend as fort_dump_read
    else:
        raise ValueError("How the hell did you get here?")

    #    file=open(fname,'rb')
    #    a=array.array('f',file.read())

    #    swap=fromstring(a,numpy.float32).byteswap()

    #    cobblers=open('SWAP001','wb')
    #    cobblers.write(swap)
    #    cobblers.close()

    #    fswap='SWAP001'

    data=fort_dump_read.fort_dump_read(ichk,fname)
    
    #intialise everything

    umassi=0.0
    udisti=0.0
    utimei=0.0
    npart=0
    tkin=0.0
    tgrav=0.0
    tterm=0.0
    escap=0.0
    rhozero=0.0
    gt=0.0

    idim=1000000
    xyzmh=zeros((5,idim))
    vxyzu=zeros((4,idim))
    rho=zeros((idim))
    poten=zeros((idim))
    iphase=zeros((idim),dtype=int8)

    #move variables from data object into more sensible structures

    ichk=data[0]
    umassi=data[1]
    udisti=data[2]
    utimei=data[3]
    npart=data[4]
    gt=data[5]
    tkin=data[6]
    tgrav=data[7]
    tterm=data[8]
    escap=data[9]
    rhozero=data[10]

    xyzmh=data[11]
    vxyzu=data[12]
    rho=data[13]
    poten=data[14]
    iphase=data[15]

    #unpack and trim arrays

    x=xyzmh[0,:]
    x=x[0:npart]
    y=xyzmh[1,:]
    y=y[0:npart]
    z=xyzmh[2,:]
    z=z[0:npart]
    m=xyzmh[3,:]
    m=m[0:npart]
    h=xyzmh[4,:]
    h=h[0:npart]

    vx=vxyzu[0,:]
    vx=vx[0:npart]
    vy=vxyzu[1,:]
    vy=vy[0:npart]
    vz=vxyzu[2,:]
    vz=vz[0:npart]
    u=vxyzu[3,:]
    u=u[0:npart]
    
    rho=rho[0:npart]
    poten=poten[0:npart]
    iphase=iphase[0:npart]

    

    #return to calling routine

    return umassi,utimei,udisti,npart,gt,tkin,tgrav,tterm,escap,rho,poten,x,y,z,m,h,vx,vy,vz,u,iphase
