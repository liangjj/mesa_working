*deck xqnu
      subroutine xqnu (nu1, nu2, mu1, theta, x, sx, id, pqa, ipqa,
     1   ierror)
c***begin prologue  xqnu
c***subsidiary
c***purpose  to compute the values of legendre functions for xlegf.
c            method: backward nu-wise recurrence for q(mu,nu,x) for
c            fixed mu to obtain q(mu1,nu1,x), q(mu1,nu1+1,x), ...,
c            q(mu1,nu2,x).
c***library   slatec
c***category  c3a2, c9
c***type      single precision (xqnu-s, dxqnu-d)
c***keywords  legendre functions
c***author  smith, john m., (nbs and george mason university)
c***routines called  xadd, xadj, xpqnu
c***revision history  (yymmdd)
c   820728  date written
c   890126  revised to meet slatec cml recommendations.  (dwl and jms)
c   901019  revisions to prologue.  (dwl and wrb)
c   901106  corrected order of sections in prologue and added type
c           section.  (wrb)
c   920127  revised purpose section of prologue.  (dwl)
c***end prologue  xqnu
      dimension pqa(*),ipqa(*)
      real dmu,nu,nu1,nu2,pq,pqa,pq1,pq2,sx,x,x1,x2
      real theta,pql1,pql2
c***first executable statement  xqnu
      ierror=0
      k=0
      pq2=0.0
      ipq2=0
      pql2=0.0
      ipql2=0
      if(mu1.eq.1) go to 290
      mu=0
c
c        call xpqnu to obtain q(0.,nu2,x) and q(0.,nu2-1,x)
c
      call xpqnu(nu1,nu2,mu,theta,id,pqa,ipqa,ierror)
      if (ierror.ne.0) return
      if(mu1.eq.0) return
      k=(nu2-nu1+1.5)
      pq2=pqa(k)
      ipq2=ipqa(k)
      pql2=pqa(k-1)
      ipql2=ipqa(k-1)
  290 mu=1
c
c        call xpqnu to obtain q(1.,nu2,x) and q(1.,nu2-1,x)
c
      call xpqnu(nu1,nu2,mu,theta,id,pqa,ipqa,ierror)
      if (ierror.ne.0) return
      if(mu1.eq.1) return
      nu=nu2
      pq1=pqa(k)
      ipq1=ipqa(k)
      pql1=pqa(k-1)
      ipql1=ipqa(k-1)
  300 mu=1
      dmu=1.
  320 continue
c
c        forward recurrence in mu to obtain q(mu1,nu2,x) and
c              q(mu1,nu2-1,x) using
c              q(mu+1,nu,x)=-2.*mu*x*sqrt(1./(1.-x**2))*q(mu,nu,x)
c                   -(nu+mu)*(nu-mu+1.)*q(mu-1,nu,x)
c
c              first for nu=nu2
c
      x1=-2.*dmu*x*sx*pq1
      x2=(nu+dmu)*(nu-dmu+1.)*pq2
      call xadd(x1,ipq1,-x2,ipq2,pq,ipq,ierror)
      if (ierror.ne.0) return
      call xadj(pq,ipq,ierror)
      if (ierror.ne.0) return
      pq2=pq1
      ipq2=ipq1
      pq1=pq
      ipq1=ipq
      mu=mu+1
      dmu=dmu+1.
      if(mu.lt.mu1) go to 320
      pqa(k)=pq
      ipqa(k)=ipq
      if(k.eq.1) return
      if(nu.lt.nu2) go to 340
c
c              then for nu=nu2-1
c
      nu=nu-1.
      pq2=pql2
      ipq2=ipql2
      pq1=pql1
      ipq1=ipql1
      k=k-1
      go to 300
c
c         backward recurrence in nu to obtain
c              q(mu1,nu1,x),q(mu1,nu1+1,x),....,q(mu1,nu2,x)
c              using
c              (nu-mu+1.)*q(mu,nu+1,x)=
c                       (2.*nu+1.)*x*q(mu,nu,x)-(nu+mu)*q(mu,nu-1,x)
c
  340 pq1=pqa(k)
      ipq1=ipqa(k)
      pq2=pqa(k+1)
      ipq2=ipqa(k+1)
  350 if(nu.le.nu1) return
      k=k-1
      x1=(2.*nu+1.)*x*pq1/(nu+dmu)
      x2=-(nu-dmu+1.)*pq2/(nu+dmu)
      call xadd(x1,ipq1,x2,ipq2,pq,ipq,ierror)
      if (ierror.ne.0) return
      call xadj(pq,ipq,ierror)
      if (ierror.ne.0) return
      pq2=pq1
      ipq2=ipq1
      pq1=pq
      ipq1=ipq
      pqa(k)=pq
      ipqa(k)=ipq
      nu=nu-1.
      go to 350
      end
