*deck @(#)onecc.f	1.1 9/8/91
c***begin prologue     m7001
c***date written       920525   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c***keywords           m7001, link 7001, spline
c***author             schneider, b. (nsf)
c***source             m7001
c***purpose            complex-complex one electron integrals. 
c***
c***description        compute all standard one electron integrals
c***                   between complex orbitals using numerical quadrature.
c***                   vectorized for efficiency.
c***                 
c***                 
c***references       
c
c***routines called    cebtc(mylib),cvmul(math),iosys(io),sscal(clams)
c***end prologue       m7001
      subroutine onecc(fnsc,ddfnsc,s,t,v,scr,mask,pt,z,n,nbfc,prnt)
      implicit integer (a-z)
      real *8 pt, z, mask
      complex *16 fnsc, ddfnsc, s, t, v, scr
      character *80 title
      logical prnt
      common /io/ inp, iout
      dimension  fnsc(n,nbfc), ddfnsc(n,nbfc), pt(n), scr(n,nbfc)
      dimension s(nbfc,nbfc), t(nbfc,nbfc), v(nbfc,nbfc)
      dimension mask(nbfc,nbfc)
c**********************************************************************c
c                       overlap                                        c
c**********************************************************************c
      call cebtc(s,fnsc,fnsc,nbfc,n,nbfc)
      call cvmul(s,s,mask,mask,nbfc*nbfc,'real')
      if (prnt) then
          title='complex-complex overlap integrals'
          call prntcm(title,s,nbfc,nbfc,nbfc,nbfc,iout)
      endif
c**********************************************************************c
c                       kinetic energy                                 c
c**********************************************************************c
      call cebtc(t,fnsc,ddfnsc,nbfc,n,nbfc)
      call cvmul(t,t,mask,mask,nbfc*nbfc,'real')
      call iosys ('write real "cc kinetic energy '//
     1            'integrals" to atomci',2*nbfc*nbfc,t,0,' ')
      if (prnt) then
          title='complex-complex kinetic energy integrals'
          call prntcm(title,t,nbfc,nbfc,nbfc,nbfc,iout)
      endif
c**********************************************************************c
c                       potential energy                               c
c**********************************************************************c
      do 10 i=1,nbfc
         do 20 j=1,n
            scr(j,i)=fnsc(j,i)/pt(j)
   20    continue
   10 continue     
      call cebtc(v,fnsc,scr,nbfc,n,nbfc)
      call sscal(2*nbfc*nbfc,-z,v,1)
      call cvmul(v,v,mask,mask,nbfc*nbfc,'real')
      call iosys ('write real "cc potential energy '//
     1            'integrals" to atomci',2*nbfc*nbfc,v,0,' ')
      if (prnt) then
          title='complex-complex potential energy integrals'
          call prntcm(title,v,nbfc,nbfc,nbfc,nbfc,iout)
      endif
      return
      end






