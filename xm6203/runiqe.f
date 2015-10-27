*deck @(#)runiqe.f	1.2  10/27/94
c***begin prologue     runiqe
c***date written       930922   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c***keywords           runiqe, link m6203
c***author             schneider, barry (nsf)
c***source             m6203
c***purpose            unique radial points
c***references         none
c
c***routines called
c***end prologue       runiqe
      subroutine runiqe (r,ru,rtyp,n)
      implicit integer (a-z)
      real*8 r, ru
      character*(*) rtyp
      common /io/inp, iout
      dimension r(n), ru(*)
      if (rtyp.eq.'newton-cotes') then
          do 10 i=2,n
             ru(i-1)=r(i)
   10     continue
      else 
          do 20 i=1,n
             ru(i)=r(i)
   20     continue
      endif          
      return
      end

