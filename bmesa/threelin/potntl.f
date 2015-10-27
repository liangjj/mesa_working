c $Header: potntl.f,v 1.2 92/12/12 09:35:06 bis Exp $
*deck potntl.f
c***begin prologue     potntl
c***date written       910128   (yymmdd)
c***revision date               (yymmdd)
c***keywords           potential
c***author             schneider, barry(lanl)
c***source             @(#)m6020
c***purpose            calculate potential for on grid
c***
c***
c***references
c
c***routines called    util
c***end prologue
      subroutine potntl(v,rmin,rmax,stp,n,type)
      implicit integer (a-z)
      dimension v(n,n,n)
      real *8 v, rmin, rmax, stp, xstp, ystp, zstp, xx, yy, zz
      character *(*) type
      common /io/ inp, iout
      if (type.eq.'exponential') then
          xstp=rmin          
          do 10 x=1,n
             xx=xstp*xstp
             ystp=rmin   
             do 20 y=1,n
                yy=ystp*ysp
                zstp=rmin
                do 30 z=1,n
                   zz=zstp*zstp
                   rval=sqrt( xx + yy + zz )
                   v(x,y,z) = -exp(-rval)
                   zstp=zstp+stp
   30           continue
                ystp=ystp+stp
   20        continue
             xstp=xstp+stp
   10     continue
      elseif (type.eq.'none') then
          call rzero(v,n*n*n)
      elseif (type.eq.'one') then
          do 40 x=1,n
             do 50 y=1,n
                do 60 z=1,n
                   v(x,y,z)=-1.d0
   60           continue
   50        continue
   40     continue
      elseif (type.eq.'half') then
          do 70 x=1,n
             do 80 y=1,n
                do 90 z=1,n
                   v(x,y,z)=.5d0
   90           continue
   80        continue
   70     continue
      elseif (type.eq.'yukawa') then
          xstp=rmin          
          do 100 x=1,n
             xx=xstp*xstp
             ystp=rmin   
             do 200 y=1,n
                yy=ystp*ysp
                zstp=rmin
                do 300 z=1,n
                   zz=zstp*zstp
                   rval=sqrt( xx + yy + zz )
                   v(x,y,z) = -exp(-rval)/rval
                   zstp=zstp+stp
  300           continue
                ystp=ystp+stp
  200        continue
             xstp=xstp+stp
  100     continue
      else        
          call lnkerr('error in call to potential')
      endif
      write(iout,*) '     potential type   = ',type 
      return
      end
