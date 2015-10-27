*deck dxpmup
      subroutine dxpmup (nu1, nu2, mu1, mu2, pqa, ipqa, ierror)
c***begin prologue  dxpmup
c***subsidiary
c***purpose  to compute the values of legendre functions for dxlegf.
c            this subroutine transforms an array of legendre functions
c            of the first kind of negative order stored in array pqa
c            into legendre functions of the first kind of positive
c            order stored in array pqa. the original array is destroyed.
c***library   slatec
c***category  c3a2, c9
c***type      double precision (xpmup-s, dxpmup-d)
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
c***end prologue  dxpmup
      double precision dmu,nu,nu1,nu2,pqa,prod
      dimension pqa(*),ipqa(*)
c***first executable statement  dxpmup
      ierror=0
      nu=nu1
      mu=mu1
      dmu=mu
      n=int(nu2-nu1+.1d0)+(mu2-mu1)+1
      j=1
      if(mod(real(nu),1.).ne.0.) go to 210
  200 if(dmu.lt.nu+1.d0) go to 210
      pqa(j)=0.d0
      ipqa(j)=0
      j=j+1
      if(j.gt.n) return
c        increment either mu or nu as appropriate.
      if(nu2-nu1.gt..5d0) nu=nu+1.d0
      if(mu2.gt.mu1) mu=mu+1
      go to 200
c
c        transform p(-mu,nu,x) to p(mu,nu,x) using
c        p(mu,nu,x)=(nu-mu+1)*(nu-mu+2)*...*(nu+mu)*p(-mu,nu,x)*(-1)**mu
c
  210 prod=1.d0
      iprod=0
      k=2*mu
      if(k.eq.0) go to 222
      do 220 l=1,k
      prod=prod*(dmu-nu-l)
  220 call dxadj(prod,iprod,ierror)
      if (ierror.ne.0) return
  222 continue
      do 240 i=j,n
      if(mu.eq.0) go to 225
      pqa(i)=pqa(i)*prod*(-1)**mu
      ipqa(i)=ipqa(i)+iprod
      call dxadj(pqa(i),ipqa(i),ierror)
      if (ierror.ne.0) return
  225 if(nu2-nu1.gt..5d0) go to 230
      prod=(dmu-nu)*prod*(-dmu-nu-1.d0)
      call dxadj(prod,iprod,ierror)
      if (ierror.ne.0) return
      mu=mu+1
      dmu=dmu+1.d0
      go to 240
  230 prod=prod*(-dmu-nu-1.d0)/(dmu-nu-1.d0)
      call dxadj(prod,iprod,ierror)
      if (ierror.ne.0) return
      nu=nu+1.d0
  240 continue
      return
      end
