      SUBROUTINE fort_dump_read(ichk,fname,umassi,udisti,
     &utimei,npart,gt,tkin,tgrav,tterm,escap,rhozero,xyzmh,
     &vxyzu,rho,poten,iphase,listpm)

cf2py intent(in,out):: ichk
cf2py intent(in):: fname

cf2py intent(out):: umassi, udisti, utimei
cf2py intent(out):: npart
cf2py intent(out):: gt,tkin,tgrav,tterm,escap,rhozero
cf2py intent(out):: xyzmh(5,idim)
cf2py intent(out):: vxyzu(5,idim)
cf2py intent(out):: rho(idim)
cf2py intent(out):: poten(idim)
cf2py intent(out):: iphase(idim)
cf2py intent(out):: listpm(iptdim)

      INCLUDE 'idim'
      INCLUDE 'COMMONS/typef'
      INCLUDE 'COMMONS/cgas'
      INCLUDE 'COMMONS/kerne'
      INCLUDE 'COMMONS/bodys'
c      INCLUDE 'COMMONS/ptmass'
      INCLUDE 'COMMONS/binary'
      INCLUDE 'COMMONS/timei'
      INCLUDE 'COMMONS/stepopt'
      INCLUDE 'COMMONS/curlist'
      INCLUDE 'COMMONS/numpa'
      INCLUDE 'COMMONS/treecom_P'
      INCLUDE 'COMMONS/radtrans'
      INCLUDE 'COMMONS/mhd'
      INCLUDE 'COMMONS/units'


      REAL*8 umassi, udisti, utimei
      REAL*8 xyzmh(5,idim), vxyzu(4,idim)
      REAL*4 rho(idim)
      REAL*8 gt,dtmaxdp
      REAL*4 poten(idim), dq(idim)
      REAL*8 trotz, trotx, tkin, tgrav, tterm, tmag, trad
      REAL*4 dgrav(idim)
      REAL*8 escap
      REAL*8 rhozero, RK2, rhocrit, rhocrit2, rhocrit3
      INTEGER*1 iphase(idim)
      CHARACTER*100 fileident
      CHARACTER*100 fname
      INTEGER*4 int1, int2, int1i, int2i, int3i
      INTEGER*8 number8
      DIMENSION nums(8)
      INTEGER nptmass, listpm(iptdim),iptmass
      REAL*8 spinx(iptdim), spiny(iptdim), spinz(iptdim), hacc, haccall,
     &ptmcrit, radcrit,angaddx(iptdim), angaddy(iptdim),
     &angaddz(iptdim), spinadx(iptdim),spinady(iptdim), spinadz(iptdim)

      REAL*8 r1i

      ichk=1
      idisk1=10


      OPEN(idisk1,file=fname,form='unformatted',convert='big_endian')
c
c--Dump file
c-------------
c
c--Standard numbers
c
      int1 = 690706
      int2 = 780806
c
c--Read ouput file
c


      READ (idisk1, END=100) int1i,r1i,int2i,i1i,int3i
c      READ (idisk1, END=100) int1i
c      write(*,*)'hello ',int1i
      IF (int1i.NE.int1) THEN
         WRITE (*,*) 'ERROR 1 in rdump: ENDIANNESS wrong?'
         GOTO 100
      ENDIF
      IF (int2i.NE.int2) THEN
         WRITE (*,*) 'ERROR 2 in rdump: default integer size wrong'
         GOTO 100
      ENDIF
      IF (int3i.NE.int1) THEN
         WRITE (*,*) 'ERROR 3 in rdump: default real size wrong'
         GOTO 100
      ENDIF
      READ (idisk1, END=100) fileident
c
c--Single values
c
c--Default int
      READ (idisk1, END=100) number
      IF (number.NE.6) THEN
         WRITE (*,*) 'ERROR 4 in rdump'
         GOTO 100
      ENDIF
      READ (idisk1, END=100) npart,n1,n2,nreassign,naccrete,nkill
c--int*1, int*2, int*4, int*8
      DO i = 1, 4
         READ (idisk1, END=100) number
      END DO
c--Default real
      READ (idisk1, END=100) number
      IF (number.NE.14) THEN
         WRITE (*,*) 'ERROR 5 in rdump'
         GOTO 100
      ENDIF
      READ (idisk1, END=100) gt, dtmaxdp, gamma, rhozero, RK2,
     &     escap, tkin, tgrav, tterm, anglostx, anglosty, anglostz,
     &     specang, ptmassin

c      WRITE(*,*)'FORTRAN: ',gamma,rhozero,RK2,escap,tkin

c--real*4
      READ (idisk1, END=100) number
c--real*8
      READ (idisk1, END=100) number
      IF (number.NE.3) THEN
         WRITE (*,*) 'ERROR 6 in rdump'
         GOTO 100
      ENDIF
      READ (idisk1, END=100) udisti, umassi, utimei
c
c--Arrays
c
c--Number of array lengths
c
      READ (idisk1, END=100) number
      IF (number.LT.2 .OR. number.GT.4) THEN
         WRITE (*,*) 'ERROR 7 in rdump'
         GOTO 100
      ENDIF
c
c--Read array type 1 header
c
      READ (idisk1, END=100) number8, (nums(i), i=1,8)
      IF (number8.NE.npart) THEN
         WRITE (*,*) 'ERROR 8 in rdump: npart wrong'
         GOTO 100
      ENDIF
      npart = number8
c
c--Read array type 2 header
c
      READ (idisk1, END=100) number8, (nums(i), i=1,8)
      nptmass = number8
c
c--Read array type 3 header
c
      IF (number.GE.3) THEN
         READ (idisk1, END=100) number8, (nums(i), i=1,8)
         IF (number8.GT.iradtrans .OR. number8.NE.1 .AND.
     &        number8.NE.npart) THEN
            WRITE (*,*) 'ERROR 9 in rdump: iradtrans wrong ',number8,
     &           iradtrans,npart
            GOTO 100
         ENDIF
         nradtrans = number8
      ENDIF
c
c--Read array type 4 header
c
      IF (number.GE.4) THEN
         READ (idisk1, END=100) number8, (nums(i), i=1,8)
         IF (number8.GT.imhd .OR. number8.NE.1 .AND.
     &        number8.NE.npart) THEN
            WRITE (*,*) 'ERROR 10 in rdump: imhd wrong ',number8,
     &           imhd,npart
            GOTO 100
         ENDIF
         nmhd = number8
      ENDIF
c
c--Read array type 1 arrays
c
c--Default int
      READ (idisk1, END=100) (isteps(i), i=1, npart)
c--int*1
      READ (idisk1, END=100) (iphase(i), i=1, npart)
c--int*2

c--int*4

c--int*8

c--Default real
      DO j = 1, 5
         READ (idisk1, END=100) (xyzmh(j,i), i=1, npart)
      END DO
      DO j = 1, 4
         READ (idisk1, END=100) (vxyzu(j,i), i=1, npart)
      END DO
c--real*4
      READ (idisk1, END=100) (rho(i), i=1, npart)
      READ (idisk1, END=100) (poten(i), i=1, npart)
c      READ (idisk1, END=100) (dgrav(i), i=1, npart)
c     READ (idisk1, END=100) (alphaMM(i), i=1, npart)
c--real*8

c
c--Read array type 2 arrays
c
c--Default int
      READ (idisk1, END=100) (listpm(i), i=1,nptmass)
c--int*1

c--int*2

c--int*4

c--int*8

c--Default real
c      READ (idisk1, END=100) (spinx(i),i=1,nptmass)
c      READ (idisk1, END=100) (spiny(i),i=1,nptmass)
c      READ (idisk1, END=100) (spinz(i),i=1,nptmass)
c      READ (idisk1, END=100) (angaddx(i),i=1,nptmass)
c      READ (idisk1, END=100) (angaddy(i),i=1,nptmass)
c      READ (idisk1, END=100) (angaddz(i),i=1,nptmass)
c      READ (idisk1, END=100) (spinadx(i),i=1,nptmass)
c      READ (idisk1, END=100) (spinady(i),i=1,nptmass)
c      READ (idisk1, END=100) (spinadz(i),i=1,nptmass)
c--real*4

c--real*8

      IF (number.GE.3 .AND. nradtrans.GT.1) THEN
c
c--Array length 3 arrays
c
c--Default int

c--int*1

c--int*2

c--int*4

c--int*8

c--Default real
         DO j = 1, 5
            READ (idisk1, END=100) (ekcle(j,i), i=1, npart)
         END DO
c--real*4

c--real*8

      ENDIF
      IF (number.GE.4 .AND. nmhd.GT.1) THEN
c
c--Array length 4 arrays
c
c--Default int

c--int*1

c--int*2

c--int*4

c--int*8

c--Default real
         DO j = 1, 3
            READ (idisk1, END=100) (Bevolxyz(j,i), i=1, npart)
         END DO
c--real*4

c--real*8

      ENDIF
c
c--End reading of dump file
c--------------------------
c

      ichk=0

      WRITE(*,*)'DUMP READ OKAY'

      CLOSE(10)

      RETURN


 100  ichk = 1

      CLOSE(10)

      WRITE(*,*)'ERROR READING DUMP - STOP'

      RETURN

      END
