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
      subroutine potntl(v,r,n,offset,type,prnt)
      implicit integer (a-z)
      dimension v(*), r(*)
      real *8 v, r
      character *(*) type
      character*80 title
      logical prnt
      common /io/ inp, iout
      call rzero(v,n+offset)
      if (type.eq.'exponential') then
          do 10 i=1,n
             v(i) = - exp(-r(i))
   10     continue
      elseif (type.eq.'none') then
          return
      elseif (type.eq.'one') then
          do 20 i=1,n
             v(i)=-1.d0
   20     continue
      elseif (type.eq.'half') then
          do 30 i=1,n
             v(i)=.5d0
   30     continue
      else
          call lnkerr('error in call to potential')
      endif
      write(iout,*) '     potential type   = ',type
      if (prnt) then
          title='potential on grid'
          call prntrm(title,v,n+offset,1,n+offset,1,iout)
      endif         
      return
      end







