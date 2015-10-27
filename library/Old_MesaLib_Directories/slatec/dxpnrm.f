*deck dxpnrm
      subroutine dxpnrm (nu1, nu2, mu1, mu2, pqa, ipqa, ierror)
c***begin prologue  dxpnrm
c***subsidiary
c***purpose  to compute the values of legendre functions for dxlegf.
c            this subroutine transforms an array of legendre functions
c            of the first kind of negative order stored in array pqa
c            into normalized legendre polynomials stored in array pqa.
c            the original array is destroyed.
c***library   slatec
c***category  c3a2, c9
c***type      double precision (xpnrm-s, dxpnrm-d)
c***keywords  legendre functions
c***author  smith, john m., (nbs and george mason university)
c***routines called  dxadj
c***revision history  (yymmdd)
c   820728  date written
c   890126  revised to meet slatec cml recommendations.  (dwl and jms)
c   901019  revisions to prologue.  (dwl and wrb)
c   901106  changed all specific intrinsics to generic.  (wrb)
c           corrected order of sections in prologue and added type
c           section.  (wrb)
c   920127  revised purpose section of prologue.  (dwl)
c***end prologue  dxpnrm
      double precision c1,dmu,nu,nu1,nu2,pqa,prod
      dimension pqa(*),ipqa(*)
c***first executable statement  dxpnrm
      ierror=0
      l=(mu2-mu1)+(nu2-nu1+1.5d0)
      mu=mu1
      dmu=mu1
      nu=nu1
c
c         if mu .gt.nu, norm p =0.
c
      j=1
  500 if(dmu.le.nu) go to 505
      pqa(j)=0.d0
      ipqa(j)=0
      j=j+1
      if(j.gt.l) return
c
c        increment either mu or nu as appropriate.
c
      if(mu2.gt.mu1) dmu=dmu+1.d0
      if(nu2-nu1.gt..5d0) nu=nu+1.d0
      go to 500
c
c         transform p(-mu,nu,x) into normalized p(mu,nu,x) using
c              norm p(mu,nu,x)=
c                 sqrt((nu+.5)*factorial(nu+mu)/factorial(nu-mu))
c                              *p(-mu,nu,x)
c
  505 prod=1.d0
      iprod=0
      k=2*mu
      if(k.le.0) go to 520
      do 510 i=1,k
      prod=prod*sqrt(nu+dmu+1.d0-i)
  510 call dxadj(prod,iprod,ierror)
      if (ierror.ne.0) return
  520 do 540 i=j,l
      c1=prod*sqrt(nu+.5d0)
      pqa(i)=pqa(i)*c1
      ipqa(i)=ipqa(i)+iprod
      call dxadj(pqa(i),ipqa(i),ierror)
      if (ierror.ne.0) return
      if(nu2-nu1.gt..5d0) go to 530
      if(dmu.ge.nu) go to 525
      prod=sqrt(nu+dmu+1.d0)*prod
      if(nu.gt.dmu) prod=prod*sqrt(nu-dmu)
      call dxadj(prod,iprod,ierror)
      if (ierror.ne.0) return
      mu=mu+1
      dmu=dmu+1.d0
      go to 540
  525 prod=0.d0
      iprod=0
      mu=mu+1
      dmu=dmu+1.d0
      go to 540
  530 prod=sqrt(nu+dmu+1.d0)*prod
      if(nu.ne.dmu-1.d0) prod=prod/sqrt(nu-dmu+1.d0)
      call dxadj(prod,iprod,ierror)
      if (ierror.ne.0) return
      nu=nu+1.d0
  540 continue
      return
      end
