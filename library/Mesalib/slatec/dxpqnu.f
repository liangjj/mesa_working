*deck dxpqnu
      subroutine dxpqnu (nu1, nu2, mu, theta, id, pqa, ipqa, ierror)
c***begin prologue  dxpqnu
c***subsidiary
c***purpose  to compute the values of legendre functions for dxlegf.
c            this subroutine calculates initial values of p or q using
c            power series, then performs forward nu-wise recurrence to
c            obtain p(-mu,nu,x), q(0,nu,x), or q(1,nu,x). the nu-wise
c            recurrence is stable for p for all mu and for q for mu=0,1.
c***library   slatec
c***category  c3a2, c9
c***type      double precision (xpqnu-s, dxpqnu-d)
c***keywords  legendre functions
c***author  smith, john m., (nbs and george mason university)
c***routines called  dxadd, dxadj, dxpsi
c***common blocks    dxblk1
c***revision history  (yymmdd)
c   820728  date written
c   890126  revised to meet slatec cml recommendations.  (dwl and jms)
c   901019  revisions to prologue.  (dwl and wrb)
c   901106  changed all specific intrinsics to generic.  (wrb)
c           corrected order of sections in prologue and added type
c           section.  (wrb)
c   920127  revised purpose section of prologue.  (dwl)
c***end prologue  dxpqnu
      double precision a,nu,nu1,nu2,pq,pqa,dxpsi,r,theta,w,x,x1,x2,xs,
     1 y,z
      double precision di,dmu,pq1,pq2,factmu,flok
      dimension pqa(*),ipqa(*)
      common /dxblk1/ nbitsf
      save /dxblk1/
c
c        j0, ipsik, and ipsix are initialized in this subroutine.
c        j0 is the number of terms used in series expansion
c        in subroutine dxpqnu.
c        ipsik, ipsix are values of k and x respectively
c        used in the calculation of the dxpsi function.
c
c***first executable statement  dxpqnu
      ierror=0
      j0=nbitsf
      ipsik=1+(nbitsf/10)
      ipsix=5*ipsik
      ipq=0
c        find nu in interval [-.5,.5) if id=2  ( calculation of q )
      nu=mod(nu1,1.d0)
      if(nu.ge..5d0) nu=nu-1.d0
c        find nu in interval (-1.5,-.5] if id=1,3, or 4  ( calc. of p )
      if(id.ne.2.and.nu.gt.-.5d0) nu=nu-1.d0
c        calculate mu factorial
      k=mu
      dmu=mu
      if(mu.le.0) go to 60
      factmu=1.d0
      if=0
      do 50 i=1,k
      factmu=factmu*i
   50 call dxadj(factmu,if,ierror)
      if (ierror.ne.0) return
   60 if(k.eq.0) factmu=1.d0
      if(k.eq.0) if=0
c
c        x=cos(theta)
c        y=sin(theta/2)**2=(1-x)/2=.5-.5*x
c        r=tan(theta/2)=sqrt((1-x)/(1+x)
c
      x=cos(theta)
      y=sin(theta/2.d0)**2
      r=tan(theta/2.d0)
c
c        use ascending series to calculate two values of p or q
c        for use as starting values in recurrence relation.
c
      pq2=0.0d0
      do 100 j=1,2
      ipq1=0
      if(id.eq.2) go to 80
c
c        series for p ( id = 1, 3, or 4 )
c        p(-mu,nu,x)=1./factorial(mu)*sqrt(((1.-x)/(1.+x))**mu)
c                *sum(from 0 to j0-1)a(j)*(.5-.5*x)**j
c
      ipq=0
      pq=1.d0
      a=1.d0
      ia=0
      do 65 i=2,j0
      di=i
      a=a*y*(di-2.d0-nu)*(di-1.d0+nu)/((di-1.d0+dmu)*(di-1.d0))
      call dxadj(a,ia,ierror)
      if (ierror.ne.0) return
      if(a.eq.0.d0) go to 66
      call dxadd(pq,ipq,a,ia,pq,ipq,ierror)
      if (ierror.ne.0) return
   65 continue
   66 continue
      if(mu.le.0) go to 90
      x2=r
      x1=pq
      k=mu
      do 77 i=1,k
      x1=x1*x2
   77 call dxadj(x1,ipq,ierror)
      if (ierror.ne.0) return
      pq=x1/factmu
      ipq=ipq-if
      call dxadj(pq,ipq,ierror)
      if (ierror.ne.0) return
      go to 90
c
c        z=-ln(r)=.5*ln((1+x)/(1-x))
c
   80 z=-log(r)
      w=dxpsi(nu+1.d0,ipsik,ipsix)
      xs=1.d0/sin(theta)
c
c        series summation for q ( id = 2 )
c        q(0,nu,x)=sum(from 0 to j0-1)((.5*ln((1+x)/(1-x))
c    +dxpsi(j+1,ipsik,ipsix)-dxpsi(nu+1,ipsik,ipsix)))*a(j)*(.5-.5*x)**j
c
c        q(1,nu,x)=-sqrt(1./(1.-x**2))+sqrt((1-x)/(1+x))
c             *sum(from 0 t0 j0-1)(-nu*(nu+1)/2*ln((1+x)/(1-x))
c                 +(j-nu)*(j+nu+1)/(2*(j+1))+nu*(nu+1)*
c     (dxpsi(nu+1,ipsik,ipsix)-dxpsi(j+1,ipsik,ipsix))*a(j)*(.5-.5*x)**j
c
c        note, in this loop k=j+1
c
      pq=0.d0
      ipq=0
      ia=0
      a=1.d0
      do 85 k=1,j0
      flok=k
      if(k.eq.1) go to 81
      a=a*y*(flok-2.d0-nu)*(flok-1.d0+nu)/((flok-1.d0+dmu)*(flok-1.d0))
      call dxadj(a,ia,ierror)
      if (ierror.ne.0) return
   81 continue
      if(mu.ge.1) go to 83
      x1=(dxpsi(flok,ipsik,ipsix)-w+z)*a
      ix1=ia
      call dxadd(pq,ipq,x1,ix1,pq,ipq,ierror)
      if (ierror.ne.0) return
      go to 85
   83 x1=(nu*(nu+1.d0)*(z-w+dxpsi(flok,ipsik,ipsix))+(nu-flok+1.d0)
     1  *(nu+flok)/(2.d0*flok))*a
      ix1=ia
      call dxadd(pq,ipq,x1,ix1,pq,ipq,ierror)
      if (ierror.ne.0) return
   85 continue
      if(mu.ge.1) pq=-r*pq
      ixs=0
      if(mu.ge.1) call dxadd(pq,ipq,-xs,ixs,pq,ipq,ierror)
      if (ierror.ne.0) return
      if(j.eq.2) mu=-mu
      if(j.eq.2) dmu=-dmu
   90 if(j.eq.1) pq2=pq
      if(j.eq.1) ipq2=ipq
      nu=nu+1.d0
  100 continue
      k=0
      if(nu-1.5d0.lt.nu1) go to 120
      k=k+1
      pqa(k)=pq2
      ipqa(k)=ipq2
      if(nu.gt.nu2+.5d0) return
  120 pq1=pq
      ipq1=ipq
      if(nu.lt.nu1+.5d0) go to 130
      k=k+1
      pqa(k)=pq
      ipqa(k)=ipq
      if(nu.gt.nu2+.5d0) return
c
c        forward nu-wise recurrence for f(mu,nu,x) for fixed mu
c        using
c        (nu+mu+1)*f(mu,nu,x)=(2.*nu+1)*f(mu,nu,x)-(nu-mu)*f(mu,nu-1,x)
c        where f(mu,nu,x) may be p(-mu,nu,x) or if mu is replaced
c        by -mu then f(mu,nu,x) may be q(mu,nu,x).
c        note, in this loop, nu=nu+1
c
  130 x1=(2.d0*nu-1.d0)/(nu+dmu)*x*pq1
      x2=(nu-1.d0-dmu)/(nu+dmu)*pq2
      call dxadd(x1,ipq1,-x2,ipq2,pq,ipq,ierror)
      if (ierror.ne.0) return
      call dxadj(pq,ipq,ierror)
      if (ierror.ne.0) return
      nu=nu+1.d0
      pq2=pq1
      ipq2=ipq1
      go to 120
c
      end
