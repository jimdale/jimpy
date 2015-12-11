def imager(fname,xmin,xmax,ymin,ymax,zmin,zmax,logmax,logrange,logmin,ang1,ang2,ang3,iline,coldtemp,hottemp):

  import numpy
  from numpy import zeros
  from numpy import float64
  from numpy import where
  from numpy import log10

  #import imcol_python_ions
  import imcol_python_ions_vels


  import colorsys
  from colorsys import hsv_to_rgb

#define other necessary options

  islice='n'
  icolmap='roy'
  inewdump='y'
  isink='y'
  sinksize=0.0
  sinkrho=0.0
  icent='n'
  maxline=2001
  angle1=ang1
  angle2=ang2
  angle3=ang3
  idt='d'
  blah=3.14159265359
  v1=zeros((maxline,maxline))
  v2=zeros((maxline,maxline))


  coldens=zeros((maxline,maxline),float64)
  velocities=zeros((maxline,maxline),float64)

  data=imcol_python_ions_vels.imcol_python(fname,xmin,xmax,ymin,ymax,inewdump,
                               islice,isink,sinksize,sinkrho,icent,iline,
                               iline,angle1,angle2,angle3,idt,blah,v1,v2,coldtemp,hottemp)

  gmw=2.

  surfdensfac=1.991e33/3.086e18**2.

  coldens=data[1]

  coldens=coldens[0:iline,0:iline]
  
  coldens=numpy.where(coldens<=0.0,10.**(logmin-10.),coldens)
  
  vels=data[2]

  vels=vels[0:iline,0:iline]
  
  vels=numpy.where(coldens<=10.**logmin,0.0,vels)
  
  coldens=log10(coldens)

  return coldens,vels
