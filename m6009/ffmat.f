*deck @(#)ffmat.f	1.1 9/8/91
c***begin prologue     ffmat
c***date written       880423   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c***keywords           ffmat, link 1106, kohn variational
c***author             schneider, barry (lanl), rescigno, tom(llnl)  
c***source             m1106
c***purpose            effective free-free matrix elements
c***description        calculate effective free-free kohn matrix elements
c***references         schneider and rescigno, physical review
c***routines called    iosys, util and mdutil
c***end prologue       ffmat
      subroutine ffmat(a,b,c,d,n,m,dir,prnt)
      implicit integer(a-z)
      real *8 d
      complex *16 a,b,c
      character *(*) dir
      character *80 title
      logical prnt
      dimension a(n,n), b(n,m), c(m,n)
      dimension d(m,n)
      common /io/ inp, iout
      if (dir.eq.'complex') then
          call cambc(a,b,c,n,m,n)
      else
          call amcbc(a,b,d,n,m,n)
      endif
      if (prnt) then
          title='new-right-hand-side'
          call prntcm(title,a,n,n,n,n,iout)
      endif
      return
      end
