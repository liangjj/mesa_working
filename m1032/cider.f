*deck @(#)cider.f	5.1  11/6/94
      subroutine cider(num,nder,nat,lag,u,intd1e,d1e,ops,maxshl,nshell)
c
c***begin prologue     cider
c***date written       871202   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c
c***keywords           ci gradients, ci derivatives
c***author             saxe, paul (lanl)
c***source             @(#)cider.f	5.1   11/6/94
c
c***purpose            to add up the last terms for ci gradients.
c
c***description
c
c    adds in the lagrangian term to the ci gradients:
c
c       all  occ
c       sum  sum l(r,i) * u(r,i,a)
c        r    i
c
c***references
c
c***routines called    (none)
c
c***end prologue       cider
c
      implicit integer (a-z)
c
      integer maxshl(nshell)
      real*8 lag(num,num)
      real*8 u(num,num,nder)
      real*8 intd1e(nder)
      real*8 d1e(nder)
      character*(*) ops
      logical logkey
c
      common /io/ inp,iout
c
      call iosys('read integer "last of shell" from rwf',nshell,
     #            maxshl,0,' ')
      call iosys('read real "ci lagrangian" from rwf',num**2,lag,0,' ')
      call iosys('read real "cphf solutions" from rwf',
     $     num**2*nder,u,0,' ')
      call iosys('read real "ci integral first derivatives" from rwf',
     $     nder,intd1e,0,' ')
c
      call rzero(d1e,nder)
      do 3 der=1,nder
         do 2 i=1,maxshl(nshell)
            do 1 r=1,num
               d1e(der)=d1e(der)+lag(r,i)*u(r,i,der)
 1          continue
 2       continue
 3    continue
c
      if (logkey(ops,'print=gradient=lagrangian-term',.false.,' '))
     $     then
         write (iout,4)
 4       format (/,10x,'lagrangian contribution to ci derivatives')
         call matout(d1e,3,nat,3,nat,iout)
      end if
c
      do 5 der=1,nder
         d1e(der)=d1e(der)+intd1e(der)
 5    continue
c
      call iosys('write real "ci first derivatives" to rwf',
     $     nder,d1e,0,' ')
      call iosys('write real "cartesian first derivatives" to rwf',
     $     nder,d1e,0,' ')
c
      if (logkey(ops,'print=gradient=total',.false.,' '))
     $     then
         write (iout,6)
 6       format (/,10x,'total ci first derivatives')
         call matout(d1e,3,nat,3,nat,iout)
      end if
c
c
      return
      end
