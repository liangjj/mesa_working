*deck potntl.f
c***begin prologue     potntl
c***date written       910128   (yymmdd)
c***revision date               (yymmdd)
c***keywords           potential
c***author             schneider, barry(lanl)
c***source             @(#)m6020
c***purpose            potential on grid
c***
c***
c***references
c
c***routines called    util
c***end prologue
      subroutine potntl(v,r,n,type,l,prnt)
      implicit integer (a-z)
      dimension v(0:n-1), r(0:n-1)
      real *8 v, r
      character *(*) type
      character*80 title
      logical prnt
      common /io/ inp, iout
      call rzero(v(0),n-1)
      if (type.eq.'exponential') then
          do 10 i=0,n-1
             v(i) = - exp(-r(i))
   10     continue
      elseif (type.eq.'none'.or.type.eq.'zero') then
          return
      elseif (type.eq.'one') then
          do 20 i=0,n-1
             v(i)=-1.d0
   20     continue
      elseif (type.eq.'half') then
          do 30 i=0,n-1
             v(i)=.5d0
   30     continue
      elseif (type.eq.'yukawa') then
          do 40 i=0,n-1
             v(i) = - exp(-r(i))/r(i)
   40     continue
      elseif (type.eq.'coulomb') then
          do 50 i=0,n-1
             v(i) = - 1.d0/r(i)
   50     continue
      else
          call lnkerr('error in call to potential')
      endif
      do 60 i=0,n-1
         v(i) = v(i) +l*(l+1)/(2.d0*r(i)*r(i))
   60 continue     
      write(iout,*)
      write(iout,*) '     potential type   = ',type
      if (prnt) then
          title='potential on grid'
          call prntrm(title,v(0),n,1,n,1,iout)
      endif         
      return
      end







