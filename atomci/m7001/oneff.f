*deck @(#)oneff.f	1.1 9/8/91
c***begin prologue     m7001
c***date written       920525   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c***keywords           m7001, link 7001, spline
c***author             schneider, b. (nsf)
c***source             m7001
c***purpose            free-free one electron integrals
c***
c***description        compute all standard one electron integrals
c***                   between free orbitals using numerical quadrature.
c***                   integrals are computed which are never used
c***                   due to vectorization of algorithm.
c***                 
c***                 
c***                 
c***references       
c
c***routines called    cvmul(math), iosys(io), cebtc(mylib), sscal(clams)
c***end prologue       m7001
      subroutine oneff(frefn0,frefn1,ddfre0,ddfre1,s,t,v,scr,mask,
     1                 pt,z,n,sze,prnt)
      implicit integer (a-z)
      real *8 pt, z, mask
      complex *16 frefn0, ddfre0, frefn1, ddfre1, t, v, scr
      character *80 title
      logical prnt
      dimension frefn0(n,sze), ddfre0(n,sze), pt(n)
      dimension s(sze,sze), t(sze,sze), v(sze,sze), mask(sze,sze)
      dimension scr(n,sze), frefn1(n,sze), ddfre1(n,sze)
      common /io/ inp, iout
c**********************************************************************c
c                       H0 - E                                         c
c**********************************************************************c
      call cebtc(t,frefn0,ddfre0,sze,n,sze)
      call cvmul(t,t,mask,mask,sze*sze,'real')
      call iosys ('write real "fcfc kinetic energy '//
     1            'integrals" to atomci',2*sze*sze,t,0,' ')
      if (prnt) then
          title='free-free H0 - E integrals of complex-complex type'
          call prntcm(title,t,sze,sze,sze,sze,iout)
      endif
      call cebtc(t,frefn0,ddfre1,sze,n,sze)
      call cvmul(t,t,mask,mask,sze*sze,'real')
      call iosys ('write real "fcfr kinetic energy '//
     1            'integrals" to atomci',2*sze*sze,t,0,' ')
      if (prnt) then
          title='free-free H0 - E integrals of complex real type'
          call prntcm(title,t,sze,sze,sze,sze,iout)
      endif
c**********************************************************************c
c                       potential energy                               c
c**********************************************************************c
      do 10 i=1,sze
         do 20 j=1,n
            scr(j,i)=frefn0(j,i)/pt(j)
   20    continue
   10 continue     
      call cebtc(v,frefn0,scr,sze,n,sze)
      call sscal(2*sze*sze,-z,v,1)
      call cvmul(v,v,mask,mask,sze*sze,'real')
      call iosys ('write real "fcfc potential energy '//
     1            'integrals" to atomci',2*sze*sze,v,0,' ')
      if (prnt) then
          title='free-free potential energy integrals of complex '//
     1          'complex type'
          call prntcm(title,v,sze,sze,sze,sze,iout)
      endif
      call cebtc(v,scr,frefn1,sze,n,sze)
      call sscal(2*sze*sze,-z,v,1)
      call cvmul(v,v,mask,mask,sze*sze,'real')
      call iosys ('write real "fcfr potential energy '//
     1            'integrals" to atomci',2*sze*sze,v,0,' ')
      if (prnt) then
          title='free-free potential energy integrals of complex '//
     1           'real type'
          call prntcm(title,v,sze,sze,sze,sze,iout)
      endif
      return
      end






