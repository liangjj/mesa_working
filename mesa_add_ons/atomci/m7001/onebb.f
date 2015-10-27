*deck @(#)onebb.f	1.1 9/8/91
c***begin prologue     m7001
c***date written       920525   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c***keywords           m7001, link 7001, spline
c***author             schneider, b. (nsf)
c***source             m7001
c***purpose            bound-bound one electron integrals
c***
c***description        compute all standard one electron integrals
c***                   between real orbitals by numerical quadrature.
c***                   the routine is vectorized for efficiency.
c***                 
c***                 
c***references       
c
c***routines called    ebtc(mylib),vmul(math),iosys(io),sscal(clams)
c***end prologue       m7001
      subroutine onebb(fns,ddfns,lb,mb,s,t,v,scr,mask,eig,work,pt,z,
     1                 tol,nsym,n,nbf,nout,dimsym,prnt,orth)
      implicit integer (a-z)
      real *8 fns, ddfns, scr, pt, z, s, t, v, mask, eig, work, tol
      character *80 title
      logical prnt, orth
      dimension fns(n,nbf), ddfns(n,nbf), pt(n), s(nbf,nbf)
      dimension t(nbf,nbf), v(nbf,nbf), mask(nbf,nbf), scr(n,nbf)
      dimension eig(nbf), work(nbf), nsym(dimsym), lb(nbf), mb(nbf)
      common /io/ inp, iout
c**********************************************************************c
c                       overlap                                        c
c**********************************************************************c
      call ebtc(s,fns,fns,nbf,n,nbf)
      call vmul(s,s,mask,nbf*nbf)
c**********************************************************************c
c                  go to orthonormal basis                             c
c                  fns and ddfns now in new basis                      c
c**********************************************************************c
      call tomobb(fns,ddfns,mask,lb,mb,scr,s,eig,work,tol,nsym,n,nbf,
     1            nout,dimsym,prnt,orth)
      call mskone(mask,lb,mb,nout,lb,mb,nout)
c**********************************************************************c
c                      overlap                                         c
c                should be orthonormal                                 c
c**********************************************************************c
      call ebtc(s,fns,fns,nout,n,nout)
      call vmul(s,s,mask,nout*nout)
      if (prnt) then
          title='real-real overlap integrals'
          call prntrm(title,s,nout,nout,nout,nout,iout)
      endif
c**********************************************************************c
c                       kinetic energy                                 c
c**********************************************************************c
      call ebtc(t,fns,ddfns,nout,n,nout)
      call vmul(t,t,mask,nout*nout)
      call iosys ('write real "real bb kinetic energy '//
     1            'integrals" to atomci',nout*nout,t,0,' ')
      if (prnt) then
          title='real-real kinetic energy integrals'
          call prntrm(title,t,nout,nout,nout,nout,iout)
      endif
c**********************************************************************c
c                       potential energy                               c
c**********************************************************************c
      do 10 i=1,nout
         do 20 j=1,n
            scr(j,i)=fns(j,i)/pt(j)
   20    continue
   10 continue     
      call ebtc(v,fns,scr,nout,n,nout)
      call vmul(v,v,mask,nout*nout)
      call sscal(nout*nout,-z,v,1)
      call iosys ('write real "real bb potential energy '//
     1            'integrals" to atomci',nout*nout,v,0,' ')
      if (prnt) then
          title='real-real potential energy integrals'
          call prntrm(title,v,nout,nout,nout,nout,iout)
      endif
      return
      end







