*deck xpmu
      subroutine xpmu (nu1, nu2, mu1, mu2, theta, x, sx, id, pqa, ipqa,
     1   ierror)
c***begin prologue  xpmu
c***subsidiary
c***purpose  to compute the values of legendre functions for xlegf.
c            method: backward mu-wise recurrence for p(-mu,nu,x) for
c            fixed nu to obtain p(-mu2,nu1,x), p(-(mu2-1),nu1,x), ...,
c            p(-mu1,nu1,x) and store in ascending mu order.
c***library   slatec
c***category  c3a2, c9
c***type      single precision (xpmu-s, dxpmu-d)
c***keywords  legendre functions
c***author  smith, john m., (nbs and george mason university)
c***routines called  xadd, xadj, xpqnu
c***revision history  (yymmdd)
c   820728  date written
c   890126  revised to meet slatec cml recommendations.  (dwl and jms)
c   901019  revisions to prologue.  (dwl and wrb)
c   901106  changed all specific intrinsics to generic.  (wrb)
c           corrected order of sections in prologue and added type
c           section.  (wrb)
c   920127  revised purpose section of prologue.  (dwl)
c***end prologue  xpmu
      real pqa,nu1,nu2,p0,x,sx,theta,x1,x2
      dimension pqa(*),ipqa(*)
c
c        call xpqnu to obtain p(-mu2,nu,x)
c
c***first executable statement  xpmu
      ierror=0
      call xpqnu(nu1,nu2,mu2,theta,id,pqa,ipqa,ierror)
      if (ierror.ne.0) return
      p0=pqa(1)
      ip0=ipqa(1)
      mu=mu2-1
c
c        call xpqnu to obtain p(-mu2-1,nu,x)
c
      call xpqnu(nu1,nu2,mu,theta,id,pqa,ipqa,ierror)
      if (ierror.ne.0) return
      n=mu2-mu1+1
      pqa(n)=p0
      ipqa(n)=ip0
      if(n.eq.1) go to 300
      pqa(n-1)=pqa(1)
      ipqa(n-1)=ipqa(1)
      if(n.eq.2) go to 300
      j=n-2
  290 continue
c
c        backward recurrence in mu to obtain
c              p(-mu2,nu1,x),p(-(mu2-1),nu1,x),....p(-mu1,nu1,x)
c              using
c              (nu-mu)*(nu+mu+1.)*p(-(mu+1),nu,x)=
c                2.*mu*x*sqrt((1./(1.-x**2))*p(-mu,nu,x)-p(-(mu-1),nu,x)
c
      x1=2.*mu*x*sx*pqa(j+1)
      x2=-(nu1-mu)*(nu1+mu+1.)*pqa(j+2)
      call xadd(x1,ipqa(j+1),x2,ipqa(j+2),pqa(j),ipqa(j),ierror)
      if (ierror.ne.0) return
      call xadj(pqa(j),ipqa(j),ierror)
      if (ierror.ne.0) return
      if(j.eq.1) go to 300
      j=j-1
      mu=mu-1
      go to 290
  300 return
      end
