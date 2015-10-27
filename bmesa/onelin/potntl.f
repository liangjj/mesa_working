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
      subroutine potntl(v,x,n,last,type)
      implicit integer (a-z)
      dimension v(0:n), x(0:n)
      real *8 v, x
      character *(*) type
      common /io/ inp, iout
      call rzero(v(0),n+1)
      if (type.eq.'exponential') then
          do 10 i=0,last
             v(i) = -exp(-x(i))
   10     continue
      elseif (type.eq.'none') then
          return
      elseif (type.eq.'one') then
          do 20 i=0,last
             v(i)=-1.d0
   20     continue
      elseif (type.eq.'half') then
          do 30 i=0,last
             v(i)=.5d0
   30     continue
      else
          call lnkerr('error in call to potential')
      endif
      write(iout,*) '     potential type   = ',type 
      return
      end
