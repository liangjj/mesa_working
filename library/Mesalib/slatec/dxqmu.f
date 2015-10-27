*deck dxqmu
      subroutine dxqmu (nu1, nu2, mu1, mu2, theta, x, sx, id, pqa, ipqa,
     1   ierror)
c***begin prologue  dxqmu
c***subsidiary
c***purpose  to compute the values of legendre functions for dxlegf.
c            method: forward mu-wise recurrence for q(mu,nu,x) for fixed
c            nu to obtain q(mu1,nu,x), q(mu1+1,nu,x), ..., q(mu2,nu,x).
c***library   slatec
c***category  c3a2, c9
c***type      double precision (xqmu-s, dxqmu-d)
c***keywords  legendre functions
c***author  smith, john m., (nbs and george mason university)
c***routines called  dxadd, dxadj, dxpqnu
c***revision history  (yymmdd)
c   820728  date written
c   890126  revised to meet slatec cml recommendations.  (dwl and jms)
c   901019  revisions to prologue.  (dwl and wrb)
c   901106  corrected order of sections in prologue and added type
c           section.  (wrb)
c   920127  revised purpose section of prologue.  (dwl)
c***end prologue  dxqmu
      dimension pqa(*),ipqa(*)
      double precision dmu,nu,nu1,nu2,pq,pqa,pq1,pq2,sx,x,x1,x2
      double precision theta
c***first executable statement  dxqmu
      ierror=0
      mu=0
c
c        call dxpqnu to obtain q(0.,nu1,x)
c
      call dxpqnu(nu1,nu2,mu,theta,id,pqa,ipqa,ierror)
      if (ierror.ne.0) return
      pq2=pqa(1)
      ipq2=ipqa(1)
      mu=1
c
c        call dxpqnu to obtain q(1.,nu1,x)
c
      call dxpqnu(nu1,nu2,mu,theta,id,pqa,ipqa,ierror)
      if (ierror.ne.0) return
      nu=nu1
      k=0
      mu=1
      dmu=1.d0
      pq1=pqa(1)
      ipq1=ipqa(1)
      if(mu1.gt.0) go to 310
      k=k+1
      pqa(k)=pq2
      ipqa(k)=ipq2
      if(mu2.lt.1) go to 330
  310 if(mu1.gt.1) go to 320
      k=k+1
      pqa(k)=pq1
      ipqa(k)=ipq1
      if(mu2.le.1) go to 330
  320 continue
c
c        forward recurrence in mu to obtain
c                  q(mu1,nu,x),q(mu1+1,nu,x),....,q(mu2,nu,x) using
c             q(mu+1,nu,x)=-2.*mu*x*sqrt(1./(1.-x**2))*q(mu,nu,x)
c                               -(nu+mu)*(nu-mu+1.)*q(mu-1,nu,x)
c
      x1=-2.d0*dmu*x*sx*pq1
      x2=(nu+dmu)*(nu-dmu+1.d0)*pq2
      call dxadd(x1,ipq1,-x2,ipq2,pq,ipq,ierror)
      if (ierror.ne.0) return
      call dxadj(pq,ipq,ierror)
      if (ierror.ne.0) return
      pq2=pq1
      ipq2=ipq1
      pq1=pq
      ipq1=ipq
      mu=mu+1
      dmu=dmu+1.d0
      if(mu.lt.mu1) go to 320
      k=k+1
      pqa(k)=pq
      ipqa(k)=ipq
      if(mu2.gt.mu) go to 320
  330 return
      end
