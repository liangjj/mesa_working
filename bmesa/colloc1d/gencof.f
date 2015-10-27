*deck gencof.f
c***begin prologue     gencof
c***date written       950704   (yymmdd)
c***revision date               (yymmdd)
c***keywords           functions, fit
c***author             schneider, barry(nsf)
c***source             
c***purpose            calculate collocation vector
c***
c***
c***references
c
c***routines called
c***end prologue
      subroutine gencof(r,a,c,work,npt,type,prnt)
      implicit integer (a-z)
      dimension r(npt), work(npt), a(npt,npt), c(npt)
      real*8 r, a, c, work
      character*(*) type
      character*80 title
      logical prnt
      common /io/ inp, iout
      len=length(type)
      if (type.eq.'sine') then
          do 10 i=1,npt
             work(i)=sin(r(i))
 10       continue
      elseif (type.eq.'cosine') then
          do 20 i=1,npt
             work(i)=cos(r(i))
 20       continue
      elseif (type.eq.'exponential') then
          do 30 i=1,npt
             work(i)=exp(-r(i))
 30       continue
      else
          call lnkerr('error in function call')
      endif
      call ebc(c,a,work,npt,npt,1)
      if (prnt) then
          title='fitting vector for '//type(1:len)//' function'
          call prntrm(title,c,npt,1,npt,1,iout)
      endif
      return
      end







































