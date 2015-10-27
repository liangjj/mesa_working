*deck @(#)onecf.f	1.1 9/8/91
c***begin prologue     m7001
c***date written       920525   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c***keywords           m7001, link 7001, spline
c***author             schneider, b. (nsf)
c***source             m7001
c***purpose            complex-free one electron integrals
c***
c***description        compute all standard one electron integrals
c***                   between complex and free orbitals using numerical
c***                   quadrature. vectorized for efficiency.
c***                 
c***                 
c***references       
c
c***routines called    cebtc(mylib),cvmul(math),iosys(io),sscal(clams)
c***                   ecbtc(mylib)
c***end prologue       m7001
      subroutine onecf(fnsc,frefn0,frefn1,ddfre0,ddfre1,s,t,v,scr,
     1                  mask,pt,z,n,nbfc,nbff,prnt)
      implicit integer (a-z)
      real *8 pt, z, mask
      complex *16 fnsc, frefn0, ddfre0, frefn1, ddfre1, scr, s, t, v
      character *80 title
      logical prnt
      common /io/ inp, iout
      dimension fnsc(n,nbfc), frefn0(n,nbff), ddfre0(n,nbff), pt(n)
      dimension s(nbfc,nbff), t(nbfc,nbff), v(nbfc,nbff), frefn1(n,nbff)
      dimension ddfre1(n,nbff), mask(nbfc,nbff), scr(n,nbfc)
c**********************************************************************c
c                       overlap                                        c
c**********************************************************************c
      call cebtc(s,fnsc,frefn0,nbfc,n,nbff)
      call cvmul(s,s,mask,mask,nbfc*nbff,'real')
      call iosys ('write real "cfc overlap '//
     1            'integrals" to atomci',2*nbfc*nbff,s,0,' ')
      if (prnt) then
          title='complex-free overlap integrals of complex type'
          call prntcm(title,s,nbfc,nbff,nbfc,nbff,iout)
      endif
      call cebtc(s,fnsc,frefn1,nbfc,n,nbff)
      call cvmul(s,s,mask,mask,nbfc*nbff,'real')
      call iosys ('write real "cfr overlap '//
     1            'integrals" to atomci',2*nbfc*nbff,s,0,' ')
      if (prnt) then
          title='complex-free overlap integrals of real type'
          call prntcm(title,s,nbfc,nbff,nbfc,nbff,iout)
      endif
c**********************************************************************c
c                       ( H0 - E )                                     c
c**********************************************************************c
      call cebtc(t,fnsc,ddfre0,nbfc,n,nbff)
      call cvmul(t,t,mask,mask,nbfc*nbff,'real')
      call iosys ('write real "cfc kinetic energy '//
     1            'integrals" to atomci',2*nbfc*nbff,t,0,' ')
      if (prnt) then
          title='complex free H0 - E integrals of complex type'
          call prntcm(title,t,nbfc,nbff,nbfc,nbff,iout)
      endif
      call cebtc(t,fnsc,ddfre1,nbfc,n,nbff)
      call cvmul(t,t,mask,mask,nbfc*nbff,'real')
      call iosys ('write real "cfr kinetic energy '//
     1            'integrals" to atomci',2*nbfc*nbff,t,0,' ')
      if (prnt) then
          title='complex-free H0 - E integrals of real type'
          call prntcm(title,t,nbfc,nbff,nbfc,nbff,iout)
      endif
c**********************************************************************c
c                       potential energy                               c
c**********************************************************************c
      do 10 i=1,nbfc
         do 20 j=1,n
            scr(j,i)=fnsc(j,i)/pt(j)
   20    continue
   10 continue     
      call cebtc(v,scr,frefn0,nbfc,n,nbff)
      call sscal(2*nbfc*nbff,-z,v,1)
      call cvmul(v,v,mask,mask,nbfc*nbff,'real')
      call iosys ('write real "cfc potential energy '//
     1            'integrals" to atomci',2*nbfc*nbff,v,0,' ')
      if (prnt) then
          title='complex-free potential energy integrals of '//
     1           'complex type'
          call prntcm(title,v,nbfc,nbff,nbfc,nbff,iout)
      endif
      call cebtc(v,scr,frefn1,nbfc,n,nbff)
      call sscal(2*nbfc*nbff,-z,v,1)
      call cvmul(v,v,mask,mask,nbfc*nbff,'real')
      call iosys ('write real "cfr potential energy '//
     1            'integrals" to atomci',2*nbfc*nbff,v,0,' ')
      if (prnt) then
          title='complex-free potential energy integrals of real type'
          call prntcm(title,v,nbfc,nbff,nbfc,nbff,iout)
      endif
      return
      end


