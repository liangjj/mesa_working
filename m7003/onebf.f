*deck @(#)onebf.f	1.1 9/8/91
c***begin prologue     m7001
c***date written       920525   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c***keywords           m7001, link 7001, spline
c***author             schneider, b. (nsf)
c***source             m7001
c***purpose            bound-free one electron integrals
c***
c***description        compute all standard one electron integrals
***                    between real and free orbitals by numerical
c***                   quadrature. vectorized for efficiency.
c***                 
c***                 
c***references       
c
c***routines called    ebtcc(mylib),cvmul(math),iosys(io),vmul(math)
c***                   sscal(clams)
c***end prologue       m7001
      subroutine onebf(fns,frefn0,frefn1,ddfre0,ddfre1,s,t,v,scr,mask,
     1                 pt,z,n,nbfb,nbff,prnt)
      implicit integer (a-z)
      real *8 fns, scr, pt, z, mask
      complex *16 frefn0, ddfre0, frefn1, ddfre1, s, t, v
      logical prnt
      character *80 title
      common /io/ inp, iout
      dimension fns(n,nbfb), frefn0(n,nbff), ddfre0(n,nbff), pt(n)
      dimension s(nbfb,nbff), t(nbfb,nbff), v(nbfb,nbff)
      dimension mask(nbfb,nbff), scr(n,nbfb), frefn1(n,nbff)
      dimension ddfre1(n,nbff)
c**********************************************************************c
c                       overlap                                        c
c**********************************************************************c
      call ebtcc(s,fns,frefn0,nbfb,n,nbff)
      call cvmul(s,s,mask,mask,nbfb*nbff,'real')
      call iosys ('write real "bfc overlap '//
     1            'integrals" to atomci',2*nbfb*nbff,s,0,' ')
      if (prnt) then
          title='bound-free overlap integrals of complex type'
          call prntcm(title,s,nbfb,nbff,nbfb,nbff,iout)
      endif
      call ebtcc(s,fns,frefn1,nbfb,n,nbff)
      call cvmul(s,s,mask,mask,nbfb*nbff,'real')
      call iosys ('write real "bfr overlap '//
     1            'integrals" to atomci',2*nbfb*nbff,s,0,' ')
      if (prnt) then
          title='bound-free overlap integrals of real type'
          call prntcm(title,s,nbfb,nbff,nbfb,nbff,iout)
      endif
c**********************************************************************c
c                       ( H0 - E )                                     c
c**********************************************************************c
      call ebtcc(t,fns,ddfre0,nbfb,n,nbff)
      call cvmul(t,t,mask,mask,nbfb*nbfc,'real')
      call iosys ('write real "bfc kinetic energy '//
     1            'integrals" to atomci',2*nbfb*nbff,t,0,' ')
      if (prnt) then
          title='bound-free H0 -E integrals of complex type'
          call prntcm(title,t,nbfb,nbff,nbfb,nbff,iout)
      endif
      call ebtcc(t,fns,ddfre1,nbfb,n,nbff)
      call cvmul(t,t,mask,mask,nbfb*nbfc,'real')
      call iosys ('write real "bfr kinetic energy '//
     1            'integrals" to atomci',2*nbfb*nbff,t,0,' ')
      if (prnt) then
          title='bound-free H0 -E integrals of real type'
          call prntcm(title,t,nbfb,nbff,nbfb,nbff,iout)
      endif
c**********************************************************************c
c                       potential energy                               c
c**********************************************************************c
      do 10 i=1,nbfb
         do 20 j=1,n
            scr(j,i)=fns(j,i)/pt(j)
   20    continue
   10 continue     
      call ebtcc(v,scr,frefn0,nbfb,n,nbff)
      call sscal(2*nbfb*nbff,-z,v,1)
      call cvmul(v,v,mask,mask,nbfb*nbff,'real')
      call iosys ('write real "bfc potential energy '//
     1            'integrals" to atomci',2*nbfb*nbff,v,0,' ')
      if (prnt) then
          title='bound-free potential energy integrals of complex type'
          call prntcm(title,v,nbfb,nbff,nbfb,nbff,iout)
      endif
      call ebtcc(v,scr,frefn1,nbfb,n,nbff)
      call sscal(2*nbfb*nbff,-z,v,1)
      call cvmul(v,v,mask,mask,nbfb*nbff,'real')
      call iosys ('write real "bfr potential energy '//
     1            'integrals" to atomci',2*nbfb*nbff,v,0,' ')
      if (prnt) then
          title='bound-free potential energy integrals of real type'
          call prntcm(title,v,nbfb,nbff,nbfb,nbff,iout)
      endif
      return
      end


