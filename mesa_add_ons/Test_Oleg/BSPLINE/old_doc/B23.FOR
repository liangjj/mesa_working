C----------------------------------------------------------------------
C     generation of B2 and B3 splines:
C----------------------------------------------------------------------

      Open(1,file='b2.dat')
      Do i=1,31
       x=-1.6+0.1*i
       y=B2(x)
       write(1,'(2f10.5)') x,y
      End do
      Close(1)

      Open(1,file='b3.dat')
      Do i=1,41
       x=-2.1+0.1*i
       y=B3(x)
       write(1,'(2f10.5)') x,y
      End do
      Close(1)

      End

      Function B3(x)
      B3=0.0
      if(x.ge.2.0.or.x.le.-2.0) Return
      if(x.ge.-2.0.and.x.lt.-1.0) B3=(2.0+x)**3/6.0
      if(x.ge.-1.0.and.x.lt.+0.0) B3=(2.0/3.0-x**2-x**3/2)
      if(x.ge.+0.0.and.x.lt.+1.0) B3=(2.0/3.0-x**2+x**3/2)
      if(x.ge.+1.0.and.x.lt.+2.0) B3=(2.0-x)**3/6.0
      Return
      End

      Function B2(x)
      B2=0.0
      if(x.ge.1.5.or.x.le.-1.5) Return
      if(x.ge.-1.5.and.x.lt.-0.5) B2=0.5*(1.5+x)**2
      if(x.ge.-0.5.and.x.lt.+0.5) B2=(0.75-x**2)
      if(x.ge.+0.5.and.x.lt.+1.5) B2=0.5*(1.5-x)**2
      Return
      End


