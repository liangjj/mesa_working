c----------------------------------------------------------------------c
c                 coulomb function code of barnett from                c
c                             cpc library                              c
c----------------------------------------------------------------------c
*deck grncal
      subroutine grncal (rn,ln,en,charge,gr1,gr1p,gr2,gr2p)
c
      parameter (nc8=12, np8=61, ns8=4, nr8=10, no8=10, nx8=78, mx8=10)
c
      dimension fc(np8), gc(np8), fcp(np8), gcp(np8), sig(np8)
c
      complex cark, xx, eta, zlmin, fc, gc, fcp, gcp, sig, civv, sc, sc1
      complex fc1, gc1, civ, cn1, cn11, ws, cexp2, fc1p, gc1p
      data stw /1.4142135e+0/
      data zero, one, pi2 /0.e0,1.e0,1.5707963e0/
c
c this subroutine constructs the regular and irregular modulating
c functions(gr1 and gr2) and their respective derivatives(gr1p,gr2p)
c as a function of r, the radial variable, and the channel
c orbital angular momentum, l, from the products of the general
c coulomb package , coulcc , of thompson and barnett (cpc 36,
c 363(1985).
c
c these modulating functions solve the following second-order,
c ordinary differential equation:
c
c ( d2/d2x + ( 1. - (2.*eta)/x - (l*(l+1))/x*x)) gi(x,eta,l) = 0
c
c where x and eta can be complex.
c
c description of parameters :
c
c* calling variables
c
c rn = radius in a0
c ln = orbital angular momentum
c en = energy in rydbergs(k2)
c charge = residual tartget charge = z
c gr1,gr1p - regular function and derivative
c gr2,gr2p - irregular function and derivative
c
c** conventions
c
c* 1) z = 0 , k2 .ge. 0
c
c     x = k*r     eta = 0
c
c     g1 = ricatti-bessel function of order l = x*jl(x)
c     g2 = ricatti-neumann function of order l = x*nl(x)
c
c     jl(x) = spherical bessel function of order l
c     nl(x) = spherical neumann function of order l
c
c coulcc: jl(x) , nl(x) , jl(x)' , nl(x)'
c              where ' indicates a derivative wrt x (not r)
c
c asymptotic forms
c
c      g1 = sin(x - l*pi/2)
c      g2 = cos(x - l*pi/2)
c
c* 2) z = 0 , k2 .lt. 0
c
c    x = i*k*r    eta = 0    i = sqrt(-1)
c
c    g1 = sqrt(2)*((-i)**l)*x*jl(x)
c    g2 = (-(i)**l)*x*hl+(x)/sqrt(2)
c
c    jl(x) = spherical bessel function of order l (x imaginary)
c    hl+(x) = jl(x) + i*nl(x)
c    nl(x) = spherical neumann function of order l
c
c       jl and hl+ are either purely real or purely imaginary quantities
c      the scaling factors proportional to (i)**l are used to
c      obtain a real quantitiy
c
c coulcc: jl(x) , hl+(x) , jl(x)' , hl+(x)'
c
c asymptotic form
c
c     g1 = exp(+k*r)/sqrt(2)
c     g2 = exp(-k*r)/sqrt(2)
c
c* 3) z .ne. 0 , k2 .ge. 0
c
c     x = k*r    eta = -z/k
c
c     g1 = regular coulomb function of order l = fl(x)
c     g2 = irregular coulomb function of order l = gl(x)
c
c coulc c: fl(x) , gl(x) , fl(x)' , gl(x)' , sig(l)
c
c asymptotic form
c    g1 = sin(thetal)
c    g2 = cos(thetal)
c
c    thetal = x - eta*alog(2*x) - pi*l/2 + sig(l)
c    sig(l) = coulomb phase shift
c
c* 4) z .ne. 0 , k2 .lt. 0
c
c     x = i*k*r    eta = i*z/k
c
c     g1 = regular whitaker function = (-i*sqrt(2)/ws)*fl(x)
c     g2 = irregualr whitaker function = ws*hl+(x)/sqrt(2)
c
c       ws = exp(i(pi/2(l + i*eta) - sig(l)))
c       fl(x) = regaular coulomb function of imaginary argument
c       gl(x) = irregular coulomb function
c       hl+(x) = fl(x) + i*gl(x)
c
c coulcc: fl(x) , hl+(x) , fl(x)' , hl+(x)' , sig(l)
c
c asymptotic forms
c
c     g1 = exp(+kr)*((2*kr)**etaa)/sqrt((2)
c     g2 = -exp(-kr)*((2*kr)**-etaa)/sqrt((2)
c
c          etaa = -z/k
c
c***********************************************************************
c
      civ=cmplx(zero,one)
      civv=cmplx(zero,-one)
      if (charge.ne.0.) go to 10
      eta=cmplx(zero,zero)
      l=ln
      l1=l+1
      xl=float(l)
      zlmin=cmplx(xl,0.)
      e=en
      xk=sqrt(abs(e))
      if (e.lt.zero) then
      sc=cmplx(one,zero)
      if (l.ne.0) sc=civv**l
      sc1=cmplx(-one,zero)
      if (l.ne.0) sc1=-(civ**l)
      endif
      r=rn
      ark=xk*r
c
      if (e.ge.zero) then
      xx=cmplx(ark,zero)
      call coulcc (xx,eta,zlmin,1,fc,gc,fcp,gcp,sig,1,1,ifail)
      gr1=ark*real(fc(1))
      gr2=-ark*real(gc(1))
      gr1p=xk*(ark*real(fcp(1))+real(fc(1)))
      gr2p=-xk*(ark*real(gcp(1))+real(gc(1)))
      else
      xx=cmplx(zero,ark)
      call coulcc (xx,eta,zlmin,1,fc,gc,fcp,gcp,sig,11,1,ifail)
      fc1=sc*fc(1)
      gc1=sc1*gc(1)
      gr1=ark*real(fc1)*stw
      gr2=ark*real(gc1)/stw
      fc1p=fc(1)+xx*fcp(1)
      fc1p=sc*fc1p
      gc1p=gc(1)+xx*gcp(1)
      gc1p=sc1*gc1p
      gr1p=xk*real(fc1p)*stw
      gr2p=xk*real(gc1p)/stw
      endif
      go to 20
   10 l=ln
      zlmin=l
      e=en
      xk=sqrt(abs(e))
      etaa=charge/xk
      if (e.ge.zero) then
      eta=cmplx(-etaa,zero)
      else
      eta=cmplx(zero,etaa)
      cn11=civ*eta
      cn11=pi2*(zlmin+cn11)
      endif
      r=rn
      ark=xk*r
      if (e.ge.zero) then
      xx=cmplx(ark,zero)
      call coulcc (xx,eta,zlmin,1,fc,gc,fcp,gcp,sig,1,0,ifail)
      gr1=real(fc(1))+aimag(fc(1))
      gr2=real(gc(1))+aimag(gc(1))
      gr1p=xk*(real(fcp(1))+aimag(fcp(1)))
      gr2p=xk*(real(gcp(1))+aimag(gcp(1)))
      else
      xx=cmplx(zero,ark)
      call coulcc (xx,eta,zlmin,1,fc,gc,fcp,gcp,sig,11,0,ifail)
      cn1=cn11-sig(1)
      cn1=civ*cn1
      ws=cexp(cn1)
      gc1=ws*gc(1)
      gc1p=ws*civ*gcp(1)
      gr2=real(gc1)/stw
      gr2p=xk*real(gc1p)/stw
      cexp2=-civ/ws
      fc1=cexp2*fc(1)
      fc1p=civ*cexp2*fcp(1)
      gr1=real(fc1)*stw
      gr1p=xk*real(fc1p)*stw
      endif
c
   20 return
      end
