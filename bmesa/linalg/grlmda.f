      subroutine grlmda (f1,f2,gam,pt,wt,lam,npts)
c***begin prologue     grlmda
c***date written       861117   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c***keywords           greens function integral, r12 moment integral
c***                   exchange operator integral
c***author             schneider, barry (lanl)
c***source             grlmda
c***purpose            calculate radial integral for exchange operator
c***description        radial integrals needed in scattering calculation
c***                   between two functions and the electron repulsion
c***                   computed for particular value of lamda
c***                   using numerical quadrature.
c
c***references
c
c***routines called    none
c
      implicit integer(a-z)
      real *8f1, f2, pt, wt, sumf, sumb, gam
      dimension f1(npts), f2(npts), pt(npts), wt(npts), gam(npts)
c
c     ----- do forward recursion -----
      lampls=lam+1
      sumf=0.d+00
      do 10 i=1,npts
      sumf=sumf+f1(i)*f2(i)*wt(i)*(pt(i)**lam)
   10 gam(i)=sumf/(pt(i)**lampls)
c
c     ----- do backward recursion -----
c
      nn=npts-1
      sumb=0.d+00
      do 20 i=nn,1,-1
      i1=i+1
      sumb=sumb+f1(i1)*f2(i1)*wt(i1)/(pt(i1)**lampls)
      gam(i)=gam(i)+sumb*(pt(i)**lam)
   20 continue
      return
      end
