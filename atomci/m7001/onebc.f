*deck @(#)onebc.f	1.1 9/8/91
c***begin prologue     m7001
c***date written       920525   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c***keywords           m7001, link 7001, spline
c***author             schneider, b. (nsf)
c***source             m7001
c***purpose            bound-complex one electron integrals
c***
c***description        compute all standard one electron integrals
c***                   between real and complex orbitals by numerical
c***                   quadrature. routine is vectorized for efficiency.
c***                 
c***                 
c***references       
c
c***routines called    ebtcc(mylib),cvmul(math),iosys(io),,sscal(clams)
c***end prologue       m7001
      subroutine onebc(fns,fnsc,delfnc,s,t,v,scr,mask,pt,z,n,
     1                 nbfb,nbfc,prnt)
      implicit integer (a-z)
      real *8 fns, scr, pt, z, mask
      complex *16 s, t, v, fnsc, delfnc
      character *80 title
      logical prnt
      common /io/ inp, iout
      dimension fns(n,nbfb), fnsc(n,nbfc), delfnc(n,nbfc), pt(n)
      dimension s(nbfb,nbfc), t(nbfb,nbfc), v(nbfb,nbfc)
      dimension mask(nbfb,nbfc), scr(n,nbfb)
c**********************************************************************c
c                       overlap                                        c
c**********************************************************************c
      call ebtcc(s,fns,fnsc,nbfb,n,nbfc)
      call cvmul(s,s,mask,mask,nbfb*nbfc,'real')
      call iosys ('write real "bc overlap '//
     1            'integrals" to atomci',2*nbfb*nbfc,s,0,' ')
      if (prnt) then
          title='bound-complex overlap integrals'
          call prntcm(title,s,nbfb,nbfc,nbfb,nbfc,iout)
      endif
c**********************************************************************c
c                       kinetic energy                                 c
c**********************************************************************c
      call ebtcc(t,fns,delfnc,nbfb,n,nbfc)
      call cvmul(t,t,mask,mask,nbfb*nbfc,'real')
      call iosys ('write real "bc kinetic energy '//
     1            'integrals" to atomci',2*nbfb*nbfc,t,0,' ')
      if (prnt) then
          title='bound-complex kinetic energy integrals'
          call prntcm(title,t,nbfb,nbfc,nbfb,nbfc,iout)
      endif
c**********************************************************************c
c                       potential energy                               c
c**********************************************************************c
      do 10 i=1,nbfb
         do 20 j=1,n
            scr(j,i)=fns(j,i)/pt(j)
   20    continue
   10 continue     
      call ebtcc(v,scr,fnsc,nbfb,n,nbfc)
      call sscal(2*nbfb*nbfc,-z,v,1)
      call cvmul(v,v,mask,mask,nbfb*nbfc,'real')
      call iosys ('write real "bc potential energy '//
     1            'integrals" to atomci',2*nbfb*nbfc,v,0,' ')
      if (prnt) then
          title='bound-complex potential energy integrals'
          call prntcm(title,v,nbfb,nbfc,nbfb,nbfc,iout)
      endif
      return
      end















