c      SUBROUTINE imcol_python(filein,xmin,xmax,ymin,ymax,
c     &inewdump,islice,isink,sinksize,sinkrho,icent,ilinex,
c     &iliney,angle1,angle2,angle3,idt,blah)
      SUBROUTINE imcol_python(filein,xmin,xmax,ymin,ymax,
     &inewdump,islice,isink,sinksize,sinkrho,icent,ilinex,
     &iliney,angle1,angle2,angle3,idt,blah,v1,v2,temp1,temp2)




c***********************************************************
c                                                          *
c  Imaging code modified for wrapping with python          *
c                                                          *
c***********************************************************

cf2py intent(in):: filein,xmin,xmax,ymin,ymax
cf2py intent(in):: inewdump,isink,sinksize,sinkrho
cf2py intent(in):: icent,iline,angle1,angle2,angle3,idt
cf2py intent(in):: temp1,temp2

c FOR SOME STUPID FUCKING REASON, THIS CODE DOES NOT WORK IF THE
C FOLLOWING VARIABLES ARE DECLARED AS INTENT OUT ONLY

cf2py intent(in,out):: blah,v1(2001,2001),v2(2001,2001)

      INCLUDE 'idim'

      REAL*8 udist, umass, utime
      INTEGER*4 ilinex,iliney

      PARAMETER (maxline=2001)

      REAL*4 v1(maxline,maxline),v2(maxline,maxline)

      COMMON /units / udist, umass, utime
      COMMON /part / npart, x(idim), y(idim), z(idim), vx(idim),
     &               vy(idim), vz(idim), u(idim), time
      COMMON /radt/  e(idim),cv(idim)
      COMMON /dens / rho(idim), pmass(idim)
      COMMON /kerne/ h(idim)
      COMMON /aux  / gamma
      COMMON /bodys/ n1, n2
      COMMON /origi/ cmx, cmy, cmz

      COMMON /rhonm/ rhozero
      COMMON /phase/ iphase(idim),isteps(idim),
     &     denshigh, ihigh, nptmass, listpm(10000)


      REAL*8 x,y,z,vx,vy,vz,h,u,rho,pmass

      REAL*4 rhonew(idim),temp1,temp2
      INTEGER*1 iphasenew(idim)
      CHARACTER*100 fileident
      INTEGER*4 first_thing
      INTEGER*8 number8

      CHARACTER*100 filein



      CHARACTER*1 icent, isink, irt, irt2, inewdump, idt,islice
      CHARACTER*1 ihigh



      DATA pi/3.141592654/

      angle1rad = angle1*pi/180.
      angle2rad = angle2*pi/180.
      angle3rad = angle3*pi/180.

      blah=pi

c if islice = 'n', buildcol integrates through the kernel, otherwise
c it is just the usual 3-d kernel

      CALL buildcol(islice)

      izero=0

         WRITE(*,*)'Opening ',filein

         OPEN (UNIT = 11, FILE = filein, FORM = 'unformatted',
     &convert='LITTLE_ENDIAN')

         READ (11, END=20) first_thing

         if (first_thing.eq.690706) then
             rewind(11)
         else
             close(11)
             OPEN (UNIT = 11, FILE = filein, FORM = 'unformatted',
     &convert='BIG_ENDIAN')
             READ (11, END=20) first_thing
             if (first_thing.ne.690706) then
                 goto 20
             endif
         endif

c
c--process file dump
c
         IF (inewdump.EQ.'Y' .OR. inewdump.EQ.'y') THEN
            READ (11, END=20)
            READ (11, END=20) fileident
            write (*,*) fileident
            IF (fileident(1:1).EQ.'F') THEN
               READ (11, END=20) number
               READ (11, END=20) npart,n1,n2,nreassign,naccrete,nkill
               IF (npart.GT.idim) THEN
                  WRITE (*,*) 'npart.GT.idim'
                  STOP
               ENDIF
               DO i = 1, 5
                  READ (11, END=20)
               END DO
               READ (11, END=20) time, dtmaxdp, gamma, rhozero, RK2,
     &         escap, tkin, tgrav, tterm, anglostx, anglosty, anglostz,
     &              specang, ptmassin
               READ (11, END=20)
               READ (11, END=20)
               READ (11, END=20) udist, umass, utime
c     write (*,*)  udist, umass, utime, npart, n1, n2, time
               READ (11, END=20)
               READ (11, END=20)
               READ (11, END=20) number8
               nptmass = number8
               READ (11, END=20) (isteps(i), i=1, npart)
               READ (11, END=20) (iphasenew(i), i=1, npart)
               DO i = 1, npart
                  iphase(i) = iphasenew(i)
               END DO
               READ (11, END=20) (x(i), i=1, npart)
               READ (11, END=20) (y(i), i=1, npart)
               READ (11, END=20) (z(i), i=1, npart)
               READ (11, END=20) (pmass(i), i=1, npart)
               READ (11, END=20) (h(i), i=1, npart)
               READ (11, END=20) (vx(i), i=1, npart)
               READ (11, END=20) (vy(i), i=1, npart)
               READ (11, END=20) (vz(i), i=1, npart)
               READ (11, END=20) (u(i), i=1, npart)
               READ (11, END=20) (rhonew(i), i=1, npart)
               DO i = 1, npart
                  rho(i) = rhonew(i)
c                  write(*,*)pmass(i)
               END DO
               READ (11, END=20)
               READ (11, END=20) (listpm(i), i=1,nptmass)
            ELSE
               READ (11, END=20) number
               READ (11, END=20) npart,n1,n2,nreassign,naccrete,nkill
               IF (npart.GT.idim) THEN
                  WRITE (*,*) 'npart.GT.idim'
                  STOP
               ENDIF
               DO i = 1, 5
                  READ (11, END=20)
               END DO
               READ (11, END=20) time, dtmaxdp, gamma, rhozero, RK2,
     &         escap, tkin, tgrav, tterm, anglostx, anglosty, anglostz,
     &              specang, ptmassin, pmassinitial
               DO i = 1, npart
                  pmass(i) = pmassinitial
               END DO
               READ (11, END=20)
               READ (11, END=20)
               READ (11, END=20) udist, umass, utime
      write (*,*)  udist, umass, utime, npart, n1, n2, time
               READ (11, END=20)
               READ (11, END=20)
               READ (11, END=20) number8
               nptmass = number8
               IF (nptmass.GT.10000) THEN
                  WRITE (*,*) 'nptmass too big'
                  STOP
               ENDIF
               READ (11, END=20) (iphasenew(i), i=1, npart)
               DO i = 1, npart
                  iphase(i) = iphasenew(i)
               END DO
               READ (11, END=20) (x(i), i=1, npart)
               READ (11, END=20) (y(i), i=1, npart)
               READ (11, END=20) (z(i), i=1, npart)
               READ (11, END=20) (rhonew(i), i=1, npart)
               DO i = 1, npart
                  rho(i) = rhonew(i)
               END DO
               READ (11, END=20) (rhonew(i), i=1, npart)
               DO i = 1, npart
                  h(i) = rhonew(i)
               END DO
               READ (11, END=20) (listpm(i), i=1,nptmass)
            ENDIF
            CLOSE(11)
         ELSE

            IF (isink.EQ.'y' .OR. isink.EQ.'Y') THEN
               IF (irt.EQ.'y' .OR. irt.EQ.'Y') THEN
                  IF (irt2.EQ.'y' .OR. irt2.EQ.'Y') THEN
            READ (11, END = 20) udist, umass, utime, npart, n1, n2,
     1        time,gamma,rhozero,RK2,(h(i),i=1,npart),escap,tkin,tgrav,
     2        tterm,trad,(x(i),i=1,npart),(y(i),i=1,npart),
     3        (z(i),i=1,npart),
     3        (vx(i),i=1,npart),(vy(i),i=1,npart),(vz(i),i=1,npart),
     4        (u(i),i=1,npart),(e(i),i=1,npart),(pmass(i),i=1,npart),
     5      (e(i),i=1,npart),(cv(i),i=1,npart),(rho(i),i=1,npart),
     6        (e(i),i=1,npart),(e(i),i=1,npart),
     6        (e(i),i=1,npart),dtmaxdp,
     6        (isteps(i), i=1, npart),(iphase(i), i=1, npart)
            ELSE
            READ (11, END = 20) udist, umass, utime, npart, n1, n2,
     1        time,gamma,rhozero,RK2,(h(i),i=1,npart),escap,tkin,tgrav,
     2        tterm,trad,(x(i),i=1,npart),(y(i),i=1,npart),
     3        (z(i),i=1,npart),
     3        (vx(i),i=1,npart),(vy(i),i=1,npart),(vz(i),i=1,npart),
     4        (u(i),i=1,npart),(e(i),i=1,npart),(pmass(i),i=1,npart),
     5        (e(i),i=1,npart),(rho(i),i=1,npart),
     6        (e(i),i=1,npart),(e(i),i=1,npart),
     6        (e(i),i=1,npart),dtmaxdp,
     6        (isteps(i), i=1, npart),(iphase(i), i=1, npart)
            ENDIF
            ELSE
            READ (11, END = 20) udist, umass, utime, npart, n1, n2,
     1        time,gamma,rhozero,RK2,(h(i),i=1,npart),escap,tkin,tgrav,
     2        tterm,(x(i),i=1,npart),(y(i),i=1,npart),(z(i),i=1,npart),
     3        (vx(i),i=1,npart),(vy(i),i=1,npart),(vz(i),i=1,npart),
     4        (u(i),i=1,npart),(pmass(i),i=1,npart),(rho(i),i=1,
     5        npart),(e(i),i=1,npart),dtmaxdp,
     6        (isteps(i), i=1, npart),(iphase(i), i=1, npart)
            ENDIF
         ELSE
            READ (11, END = 20) udist, umass, utime, npart, n1, n2,
     1        time,gamma,rhozero,RK2,(h(i),i=1,npart),escap,tkin,tgrav,
     2        tterm,(x(i),i=1,npart),(y(i),i=1,npart),(z(i),i=1,npart),
     3        (vx(i),i=1,npart),(vy(i),i=1,npart),(vz(i),i=1,npart),
     4        (u(i),i=1,npart),(pmass(i),i=1,npart),(rho(i),i=1,
     5        npart),(e(i),i=1,npart)
         ENDIF

         ENDIF

         WRITE(*,*)'DUMP READ SUCCESSFULLY'

c turn off particles outside temp range
            ikill=0
            WRITE(*,*)'Temps: ',temp1,temp2
            DO i=1,npart
               IF(u(i).LT.temp1 .OR. u(i).GT.temp2)THEN
                 iphase(i)=-1
                 ikill=ikill+1
               END IF
            END DO

            WRITE(*,*)'killed ',ikill,' particles.'

c
c--reset coordinates to center of mass
c
         frac = 1.0
         IF (icent.EQ.'m') THEN
            CALL origin (n1, n2, frac)
         ELSEIF (icent.EQ.'d') THEN
            irhomax = 1
            DO i = 1, npart
               IF (rho(i).GT.rho(irhomax)) irhomax = i
            END DO
            WRITE (*,*) 'irhomax ',irhomax, rho(irhomax),
     &           x(irhomax), y(irhomax), z(irhomax)
            ncen = 0
            xcen = 0.
            ycen = 0.
            zcen = 0.
            DO i = 1, npart
               dx = x(i) - x(irhomax)
               dy = y(i) - y(irhomax)
               dz = z(i) - z(irhomax)
               r2 = dx*dx + dy*dy + dz*dz
               IF (r2.LT.16.0*h(irhomax)*h(irhomax)) THEN
                  ncen = ncen + 1
                  xcen = xcen + x(i)
                  ycen = ycen + y(i)
                  zcen = zcen + z(i)
               ENDIF
            END DO
            xcen = xcen/ncen
            ycen = ycen/ncen
            zcen = zcen/ncen
            WRITE (*,*) 'centre: ',ncen, xcen, ycen, zcen
            DO i = 1, npart
               x(i) = x(i) - xcen
               y(i) = y(i) - ycen
               z(i) = z(i) - zcen
            END DO
         ENDIF
c
c--Rotate by angles
c
         DO i = 1, npart
            r = SQRT(x(i)*x(i) + y(i)*y(i))

            th1 = ATAN2(y(i),x(i))

            th1 = th1 + angle1rad
            x(i) = r*COS(th1)
            y(i) = r*SIN(th1)

            r = SQRT(vx(i)*vx(i) + vy(i)*vy(i))

            th1 = ATAN2(vy(i),vx(i))

            th1 = th1 + angle1rad
            vx(i) = r*COS(th1)
            vy(i) = r*SIN(th1)

            r = SQRT(y(i)*y(i) + z(i)*z(i))

            th2 = ATAN2(z(i),y(i))

            th2 = th2 + angle2rad
            y(i) = r*COS(th2)
            z(i) = r*SIN(th2)

            r = SQRT(vy(i)*vy(i) + vz(i)*vz(i))

            th2 = ATAN2(vz(i),vy(i))

            th2 = th2 + angle2rad
            vy(i) = r*COS(th2)
            vz(i) = r*SIN(th2)

            r = SQRT(x(i)*x(i) + z(i)*z(i))

            th3 = ATAN2(z(i),x(i))

            th3 = th3 + angle3rad
            x(i) = r*COS(th3)
            z(i) = r*SIN(th3)

            r = SQRT(vx(i)*vx(i) + vz(i)*vz(i))

            th3 = ATAN2(vz(i),vx(i))

            th3 = th3 + angle3rad
            vx(i) = r*COS(th3)
            vz(i) = r*SIN(th3)
         END DO
c
c--Centre on frame
c
         xcen = (xmax+xmin)/2.0
         ycen = (ymax+ymin)/2.0
         xmaxp = xmax - xcen
         xminp = xmin - xcen
         ymaxp = ymax - ycen
         yminp = ymin - ycen

         DO i = 1, npart
            x(i) = x(i) - xcen
            y(i) = y(i) - ycen
         END DO

c
c--create images
c
         CALL grid (xminp, xmaxp, yminp, ymaxp, ilinex, iliney, idt,
     &sinksize, sinkrho,v1,v2,islice)

      RETURN

 20   WRITE(*,*)'ERROR READING DUMP'
      STOP

      END

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      SUBROUTINE buildcol(islice)

      COMMON /column/ coltable(1000)

      CHARACTER*1 islice

      DATA pi/3.14159265/


      DO i=1,1000
         r=(i-1)/500.
         dist=SQRT(4.0-r*r)
         step=dist/4000.0
         ypos=0.

         IF(islice.EQ.'y')THEN
            v=SQRT(r*r+ypos*ypos)
            IF (v.LT.1.0) THEN
               v2=v*v
               val=1.0-1.5*v2+0.75*v2*v
            ELSE
               v2m=2.0-v
               val=0.25*v2m*v2m*v2m
            ENDIF

            coltable(i)=val/pi

         ELSE
            coldens=0.0
            DO j=1,4000
               v=SQRT(r*r+ypos*ypos)
               IF (v.LT.1.0) THEN
                  v2=v*v
                  val=1.0-1.5*v2+0.75*v2*v
                  coldens=coldens+val*step
               ELSE
                  v2m=2.0-v
                  val=0.25*v2m*v2m*v2m
                  coldens=coldens+val*step
               ENDIF
               ypos=ypos+step
            END DO

            coltable(i)=2.0*coldens/pi
         END IF
      END DO

c      check=0.0
c      DO i=1,1000
c         r=(i-1)/500.
c         WRITE(85,*) r, coltable(i)
c         check=check+2*pi*r*coltable(i)*(1.0/500.)
c      END DO
c      WRITE(*,*)check

      RETURN
      END

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      SUBROUTINE grid (xmin, xmax, ymin, ymax, ilinex, iliney, idt,
     &sinksize, sinkrho,v1,v2,islice)
c
c*****************************************************************
c                                                                *
c  This subroutines interpolates quantitites on a given grid.    *
c  The quantity plotted is column density through the simulation *
c     calculated by integrating density through the SPH          *
c     smoothing kernel for all particles contributing to the     *
c     line of sight density.                                     *
c                                                                *
c*****************************************************************
c
      INCLUDE 'idim'
      PARAMETER (maxline=2001)
      PARAMETER (mline2=maxline*maxline)
      PARAMETER (npert=30)
      PARAMETER (npert1=npert+1)

      REAL*8 udist, umass, utime, udens, ucolumndens, uvel
      INTEGER*4 ilinex,iliney
      REAL*4 v1(maxline,maxline),v2(maxline,maxline)
      REAL*4 timeout, xmindensout, xmaxdensout
      CHARACTER*1 idt

      COMMON /units / udist, umass, utime
      COMMON /part / npart, x(idim), y(idim), z(idim), vx(idim),
     &               vy(idim), vz(idim), u(idim), time
      COMMON /radt/  e(idim),cv(idim)
      COMMON /partu/ npp
      COMMON /dens / rho(idim), pmass(idim)
      COMMON /tempp/ temperaturep(idim)
      COMMON /kerne/ h(idim)
      COMMON /aux  / gamma

      COMMON /bodys/ n1, n2
      COMMON /origi/ cmx, cmy, cmz
      COMMON /scale/ vmin, vmax
      COMMON /rhonm/ rhozero
      COMMON /phase/ iphase(idim),isteps(idim),
     &     denshigh, ihigh, nptmass, listpm(10000)
      COMMON /grids/ xmaxdens, xmindens,
     & xmaxdens2, xmindens2


      REAL*8 x,y,z,vx,vy,vz,h,u,rho,pmass

c
      CHARACTER*1 ihigh,islice
      REAL*4 xmaxdens, xmindens
      REAL*4 xmaxdens2, xmindens2
c
      DATA pi/3.141592654/
c
      IF (ilinex.GT.maxline.OR.iliney.GT.maxline) THEN
         WRITE(*,*) 'ERROR - Dimensions not big enough'
         STOP
      ENDIF
      udens = umass/(udist**3)
      ucolumndens = udens*udist
      uvel = udist/utime/1e5
      
c      write(*,*)'islice: ',islice
c
c--compute steps and narrow the number of particles
c
      x0 = xmin
      y0 = ymin
      stepx = (xmax - xmin) / (ilinex - 1)
      stepy = (ymax - ymin) / (iliney - 1)
      npp = 0
      nlistpm = nptmass
      nlistpmadd = 0
      DO 12 i = 1, npart
         IF (iphase(i).GE.0) THEN
            IF (x(i)+2.0*h(i).LT.xmin.OR.
     &           x(i)-2.0*h(i).GT.xmax) GOTO 12
            IF (y(i)+2.0*h(i).LT.ymin.OR.
     &           y(i)-2.0*h(i).GT.ymax) GOTO 12
      IF((islice.EQ.'y').AND.
     &(z(i).GT.12.*h(i).OR.z(i)
     &.LT.-12.*h(i)))GOTO 12
            npp = npp + 1
            x(npp) = x(i)
            y(npp) = y(i)
            z(npp) = z(i)
            vx(npp) = vx(i)
            vy(npp) = vy(i)
            vz(npp) = vz(i)
            u(npp)=u(i)
            pmass(npp) = pmass(i)
c            write(*,*)pmass(i)
            rho(npp) = rho(i)
            temperaturep(npp) = u(i)/cv(i)
            IF (iphase(i).GE.1) THEN
               DO j = 1, nptmass
                  IF (listpm(j).EQ.i) THEN
                     listpm(j)=npp
                     GOTO 9
                  ENDIF
               END DO
            ENDIF
 9          CONTINUE

c            write (*,*) i,iphase(i),rho(i),udens,denshigh

            IF (iphase(i).EQ.0)THEN
               h(npp) = h(i)
            ELSE
 11            h(npp)=sinksize
               pmass(npp) = sinkrho*3.142*(2.0*sinksize)**2
c               pmass(npp)=0.
            ENDIF
         ENDIF
 12   CONTINUE

      write(*,*)'Number of particles kept   :  ',npp
      write(*,*)'Number of sinks   :  ',nptmass

      xmaxdens = -1.0
      xmindens = 1.0E+20
      xmaxdens2 = -1.0
      xmindens2 = 1.0E+20
c
c--initialize grid
c
C$OMP PARALLEL default(none)
C$OMP& shared(ilinex,iliney,v1,v2,npp,x,y,z,pmass,rho,h)
C$OMP& shared(stepx,stepy,x0,y0,ucolumndens,temperaturep,idt)
C$OMP& private(i,j,xi,yi,zi,pmassi,rhoi,hi,hproj,ix,iy,icx,icy)
C$OMP& private(idepx,idepy,ifinx,ifiny,kx,ky,gx,gy,val,valuei)
C$OMP& reduction(MIN:xmindens)
C$OMP& reduction(MIN:xmindens2)
C$OMP& reduction(MAX:xmaxdens)
C$OMP& reduction(MAX:xmaxdens2)
C$OMP DO SCHEDULE(dynamic, 10)

      DO  i = 1, ilinex
         DO  j = 1, iliney
            v1(i,j) = 0.0
            v2(i,j) = 0.0
         END DO
      END DO

C$OMP END DO
c
c--compute values on grid
c
C$OMP DO SCHEDULE(dynamic, 10)
      DO 30 i = 1, npp
c         IF (MOD(i,1).EQ.0) WRITE (*,*) i,npp

         xi = x(i)
         yi = y(i)
         zi = z(i)
         valuei=pmass(i)
         vzi=vz(i)
c         valuei = u(i)
c         write(*,*)pmass(i)
         rhoi = rho(i)
         hi = h(i)
         hproj = 2.0*hi

c         WRITE(*,*)xi,yi,zi,hi,stepx,stepy,x0,y0
c
c--find index of closest grid point
c
         ix = nint( (xi - x0) / stepx) + 1
         iy = nint( (yi - y0) / stepy) + 1
         icx = hproj / stepx + 2
         icy = hproj / stepy + 2
         idepx = max0(1, ix - icx)
         idepy = max0(1, iy - icy)
         ifinx = min0(ilinex, ix + icx)
         ifiny = min0(iliney, iy + icy)

c         WRITE(*,*)ix,iy,icx,icy
c         WRITE(*,*)idepx,idepy,ifinx,ifiny
c
c--compute particle's contribution to all grid points
c
         DO 20 ky = idepy, ifiny
            gy = y0 + (ky - 1) * stepy
            DO 20 kx = idepx, ifinx
               gx = x0 + (kx - 1) * stepx
               gz=zi

         CALL value (xi,yi,zi,hi,gx,gy,gz,islice,val,val2,valuei,vzi)

c               WRITE(*,*)val2
C$OMP ATOMIC
               v1(kx,ky) = v1(kx,ky) + val
               v2(kx,ky) = v2(kx,ky) + val2

C$OMP ATOMIC
 20      CONTINUE
 30   CONTINUE
C$OMP END DO
c
c--create output file
c
C$OMP DO SCHEDULE(dynamic, 10)
      DO j=1, iliney
         DO i=1, ilinex
            IF (v1(i,j).GT.0.) v2(i,j)=v2(i,j)/v1(i,j)
            IF(islice.EQ.'y')THEN
               v1(i,j)=v1(i,j)*udens
            ELSE
               v1(i,j)=v1(i,j)*ucolumndens
               v2(i,j)=v2(i,j)*uvel
            END IF

            IF(v1(i,j).LT.0.)v1(i,j)=0.

            xmaxdens = MAX(xmaxdens, v1(i,j))
            xmindens = MIN(xmindens, v1(i,j))
            xmaxdens2 = MAX(xmaxdens2, v2(i,j))
            xmindens2 = MIN(xmindens2, v2(i,j))
         END DO
      END DO
C$OMP END DO
C$OMP END PARALLEL
      WRITE(*,*) 'Max col density (g cm^-2)=',xmaxdens
      WRITE(*,*) 'Min col density (g cm^-2)=',xmindens
      WRITE(*,*) 'Max temperature (K)      =',xmaxdens2
      WRITE(*,*) 'Min temperature (K)      =',xmindens2

c      v1(1,1)=10000.0

      IF (idt.EQ.'t') THEN
         timeout = time
         xmaxdensout = xmaxdens2
         xmindensout = xmindens2
         write (*,*) ilinex,iliney,timeout,xmaxdensout,xmindensout
c      WRITE(16,IOSTAT=io) iline,iline,timeout,xmaxdensout,xmindensout
c      WRITE(16,IOSTAT=io) ((v2(i,j), i=1,iline), j=1,iline)
      ELSE
         timeout = time
         xmaxdensout = xmaxdens
         xmindensout = xmindens
         write (*,*) ilinex,iliney,timeout,xmaxdensout,xmindensout
c      WRITE(16,IOSTAT=io) iline,iline
c      WRITE(16,IOSTAT=io) ((v1(i,j), i=1,iline), j=1,iline)
      ENDIF

      RETURN
      END

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      SUBROUTINE value (xi,yi,zi,hi,gx,gy,gz,islice,val,val2,
     &pmassi,vzi)
c***********************************************************
c                                                          *
c  this subroutine computes the value of the quantities at *
c  one given point.                                        *
c                                                          *
c***********************************************************
c
      INCLUDE 'idim'

      COMMON /column/ coltable(1000)
      CHARACTER*1 islice

      INTEGER pos1,pos2,pos
c
      pi=3.14159265
c
      val = 0.0
      val2 = 0.0
c
c      IF(islice.EQ.'y')THEN
c         vr2 = ((xi-gx)**2 + (yi-gy)**2+gz**2.)/hi**2
c      ELSE
         vr2 = ((xi-gx)**2 + (yi-gy)**2)/hi**2
c      END IF


      IF (vr2.GE.4.0) RETURN

      vr = SQRT(vr2)
c
c--find nearest table points
c
      pos=INT(vr*500.0)+1
      IF (pos.EQ.1000) THEN
         pos1=999
         pos2=1000
      ELSE
         pos1=pos
         pos2=pos1+1
      ENDIF
      r1=(pos1-1)/500.
      r2=(pos2-1)/500.
      valuemean = (vr-r1)/(r2-r1)*(coltable(pos2)-coltable(pos1))
      valuemean = valuemean + coltable(pos1)
c
c      write(*,*)valuemean,pmassi,hi

      IF(islice.EQ.'y')THEN
         val=valuemean*pmassi/hi/hi/hi
         val2=0.
      ELSE
         val = valuemean*pmassi/hi/hi
         val2 = valuemean*pmassi*vzi/hi/hi
      END IF
c
      RETURN
      END

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      SUBROUTINE origin (n1, n2, frac)
c************************************************************
c                                                           *
c  this subroutine computes the origin of the coordinate    *
c  system.                                                  *
c                                                           *
c************************************************************
c
      INCLUDE 'idim'
      PARAMETER (iline=1001)
c
      COMMON /part / npart, x(idim), y(idim), z(idim), vx(idim),
     1               vy(idim), vz(idim), u(idim), time
      COMMON /dens / rho(idim), pmass(idim)
      COMMON /kerne/ h(idim)
      COMMON /aux  / gamma
      COMMON /origi/ cmx, cmy, cmz
c
c--center of mass
c
      cmx = 0.
      cmy = 0.
      cmz = 0.
      tmass = 0.
c
      DO 10 i=1,npart
         tmass = tmass + pmass(i)
         cmx = cmx + pmass(i)*x(i)
         cmy = cmy + pmass(i)*y(i)
   10    cmz = cmz + pmass(i)*z(i)

      cmx = cmx/tmass
      cmy = cmy/tmass
      cmz = cmz/tmass
c
      DO 20 i=1,npart
         x(i) = x(i) - cmx*frac
         y(i) = y(i) - cmy*frac
         z(i) = z(i) - cmz*frac
   20 CONTINUE

      WRITE(*,*) cmx,cmy,cmz
c
      RETURN
      END
