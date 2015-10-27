*deck @(#)matm.f	1.1 9/8/91
c***begin prologue     matm
c***date written       880423   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c***keywords           matm, link 1106, kohn variational
c***author             schneider, barry (lanl), rescigno, tom(llnl)  
c***source             m1106
c***purpose            negative of a matrix
c***description        a = -a
c***references         schneider and rescigno, physical review
c***routines called    iosys, util and mdutil
c***end prologue       matm
      subroutine matm(a,n)
      implicit integer(a-z)
      complex *16 a
      dimension a(n)
      do 10 i=1,n
         a(i)=-a(i)
   10 continue
      return
      end
