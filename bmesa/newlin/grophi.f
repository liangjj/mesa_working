      subroutine grophi (reg,ireg,phi,gam,pt,wt,npts)
c***begin prologue     grophi
c***date written       861117   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c***keywords           greens function integral
c***                   integral of greens function with a function
c***author             schneider, barry (lanl)
c***source             grophi
c***purpose            radial integral of greens function with
c***                   arbitrary function
c***description        radial integrals needed in scattering calculation
c***                   when greens function acts on various functions
c***                   computed for particular value of l
c***                   using numerical quadrature.
c
c***references
c
c***routines called    none
      implicit integer(a-z)
      real *8reg, ireg, phi, pt, wt, sumf, sumb, gam
      dimension reg(npts), ireg(npts), pt(npts), wt(npts), gam(npts)
      dimension phi(npts)
c
c     ----- do forward recursion -----
      sumf=0.d+00
      do 10 i=1,npts
      sumf=sumf+reg(i)*phi(i)*wt(i)
   10 gam(i)=ireg(i)*sumf
c
c     ----- do backward recursion -----
c
      nn=npts-1
      sumb=0.d+00
      do 20 i=nn,1,-1
      i1=i+1
      sumb=sumb+ireg(i1)*phi(i1)*wt(i1)
      gam(i)=gam(i)+sumb*reg(i)
   20 continue
      return
      end
