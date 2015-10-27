c $Header: mkfx.f,v 1.2 92/12/12 09:34:39 bis Exp $
*deck @(#)mkfx.f	1.1 9/8/91
c***begin prologue     mkfx
c***date written       921206 (yymmdd)
c***revision date             (yymmdd)
c***keywords           
c***author             schneider, barry(lanl)
c***source             @(#)util
c***purpose            set up effective potential for numerov
c***                   integration of one-dimensional
c***                   schroedinger equation.
c***                                         ''
c***                   the equation here is y  + f(x) y = g
c***                   and we are computing f.
c***routines called  
c***end prologue
      subroutine mkfx(v,f,energy,refe,n,dir,urefe)
      implicit integer (a-z)
      real *8  v, f, energy, refe
      character*(*) dir
      logical urefe
      dimension v(0:n), f(0:n)
      common /io/ inp, iout
      if (dir.eq.'with v') then
          do 10 i=0,n
             f(i) =  energy - 2.d0*v(i)
   10     continue
      elseif (dir.eq.'without v'.and.urefe.eq..false.) then
          do 20 i=0,n
             f(i) =  energy
   20     continue
      elseif (dir.eq.'without v'.and.urefe.eq..true.) then
      write(iout,*) ' using reference energy'
          do 30 i=0,n
             f(i) =  refe
   30     continue
      endif
      return
      end



