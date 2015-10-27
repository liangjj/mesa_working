*deck @(#)fildf1.f	1.1 9/8/91
c***begin prologue     fildf1
c***date written       880423   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c***keywords           utility, matrix
c***author             schneider, barry (lanl)
c***source             m6005
c***purpose            fill matrix in packed form
c***description
c***
c***
c***references         none
c
c***routines called
c***end prologue       fildf1
      subroutine fildf1(d,a,b,nd,l,m,maxl,npt,lmax,dim,bcond)
      implicit integer (a-z)
      complex *16  b, d
      real *8 a
      character *(*) bcond
      dimension d(nd,npt), a(npt,0:lmax,0:2*lmax), b(npt,maxl)
      dimension l(dim), m(dim)
c----------------------------------------------------------------------c
c                 test boundary conditions                             c
c----------------------------------------------------------------------c
      if (bcond.eq.'s-matrix') then
          do 10 nolm1=1,nd
             l1=l(nolm1)
             m1=m(nolm1)
             do 20 grpt=1,npt
                d(nolm1,grpt)=a(grpt,l1,m1)*conjg(b(grpt,nolm1))
   20        continue
   10     continue
      elseif (bcond.eq.'t-matrix') then
           do 30 nolm1=1,nd
              l1=l(nolm1)
              m1=m(nolm1)
              do 40 grpt=1,npt
                 d(nolm1,grpt)=a(grpt,l1,m1)*imag(b(grpt,nolm1))
   40        continue
   30     continue
      else
          call lnkerr('boundary condition error in ffints')
      endif
      return
      end
