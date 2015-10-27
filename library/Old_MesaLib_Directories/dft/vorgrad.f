*deck @(#)vorgrad.f	5.2 4/17/95
      subroutine vorgrad(natoms,c,xyzgrid,wts,mxgrd,ngrid,atom,
     $                   vwts,rnuc,amu,pwtx,rr,adjust,radii,akl,
     $                   dograd,gradwt)
c***begin prologue     vorgrad.f
c***date written       930601  
c***revision date      4/17/95
c
c***keywords           
c***author             P. J. Hay            
c***source             @(#)vorgrad.f	5.2 4/17/95
c***purpose            calculates "fuzzy" voronoi weights and gradients
c***description
c   natoms     ...     number of atoms. 
c   c          ...     atomic coordinates.
c   xyzgrid    ...     grid points at which to evaluate voronoi weights.
c   wts        ...     atomic quadrature weights.
c   mxgrd      ...     maximum number of grid points.
c   ngrid      ...     number of grid points.
c   atom       ...     the sequence number of the atom whose grid we are
c                      evaluating
c   vwts       ...     voronoi weights.
c   rnuc       ...     scratch(natoms,natoms).
c   amu        ...     scratch(natoms,natoms).
c   pwtx       ...     scratch(natoms).
c   rr         ...     scratch(natoms).
c   adjust     ...     flag which determines if the cell functions are
c                      to be scaled according to atomic radii.
c   radii      ...     atomic radii used for cell adjustment.
c
c***references
c
c                      A.D. Becke, J. Chem. Phys. 88, 2547 (1988)
c***routines called
c
c***end prologue       vorgrad.f
      implicit none
c     --- input variables -----
      integer natoms,mxgrd,ngrid,atom
      logical adjust
      logical dograd
c     --- input arrays (unmodified) ---
      real*8 c(3,natoms)
      real*8 xyzgrid(mxgrd,3)
      real*8 radii(natoms)
c     --- input arrays (scratch) ---
c     --- output arrays ---
      real*8 vwts(mxgrd),wts(mxgrd)
      real*8 gradwt(mxgrd,3,natoms)
c     --- output variables ---
c     --- scratch arrays ---
      real*8 rnuc(natoms,natoms),amu(natoms,natoms)
      real*8 akl(natoms,natoms)
      real*8 pwtx(natoms),rr(natoms)
      real*8 term1(3),term2(3),term3(3),x(3)
c     --- local variables ---
      integer igrid,k,l,kparm,kk,iatom,jatom
      integer ix
      real*8 xx,yy,zz,xam,yam,denom,dist
      real*8 zero,one,half,onehf
      real*8 chi,u
      real*8 tmu,gmu,t,gradmu
      integer iout,inp
      common /io/inp,iout
      parameter(zero=0.0d0,one=1.0d0,half=0.5d0,onehf=1.5d0)
c
c     --- check for atomic case, note early return
      if(natoms.eq.1) then
         call vfill(vwts,one,ngrid)
         return
      endif
c
c     --- calculate internuclear distances 
      do 100 iatom=1,natoms
         do 100 jatom=1,natoms
            dist=sqrt((c(1,iatom)-c(1,jatom))**2
     $     	+ (c(2,iatom)-c(2,jatom))**2  
     $    	+ (c(3,iatom)-c(3,jatom))**2)
            rnuc(iatom,jatom)=dist
 100  continue 
c
c     --- calculate scaling parameters if adjusting cell radii.
      if(adjust) then
         do 105 k=1,natoms
            do 104 l=1,natoms
               chi=radii(k)/radii(l)
               u=(chi-one)/(chi+one)
               akl(k,l)=u/(u*u-one)
  104      continue
  105   continue
      endif
c
c     --- kparm is value of repeat in recursion function
      kparm=3
c
c     --- part for energy only
      if(.not.dograd)then
c
c        --- determine the hyperboloidal coordinates for this grid 
         do 110 igrid=1,ngrid
            xx=xyzgrid(igrid,1)
            yy=xyzgrid(igrid,2)
            zz=xyzgrid(igrid,3)
            do 111 k=1,natoms
               dist=(xx-c(1,k))**2 + (yy-c(2,k))**2
     $             + (zz-c(3,k))**2
               rr(k)=sqrt(dist)
 111        continue 
            do 112 k=1,natoms
               do 112 l=1,natoms
                  if (k.ne.l) then
                     amu(k,l)=(rr(k)-rr(l))/rnuc(k,l)
                     if(adjust) then
                        amu(k,l)=amu(k,l)
     $                          +akl(k,l)*(one-amu(k,l)*amu(k,l))
                     endif
                  end if
 112        continue
c
c           --- generate the cell functions
            do 120 k=1,natoms
               pwtx(k)=one
               do 130 l=1,natoms
                  if(k.ne.l) then
                     xam=amu(k,l)
                     do 131 kk=1,kparm
                        yam=onehf*xam-half*xam**3
                        xam=yam
 131                 continue 
                     yam=half*(one-yam)
                     pwtx(k)=pwtx(k)*yam
                  end if
 130           continue
 120        continue
c
c           --- normalize cell functions
            denom=zero
            do 113 k=1,natoms
               denom=denom+pwtx(k)
 113        continue 
            vwts(igrid)=pwtx(atom)/denom
 110     continue
c
c        --- multiply voronoi weights and numerical weights
         do 140 igrid=1,ngrid
            wts(igrid)=wts(igrid)*vwts(igrid)
  140    continue
      else
c
c        --- part for energy and gradient
c
c        --- determine the hyperboloidal coordinates for this grid 
         do 210 igrid=1,ngrid
            x(1)=xyzgrid(igrid,1)
            x(2)=xyzgrid(igrid,2)
            x(3)=xyzgrid(igrid,3)
            do 211 k=1,natoms
               dist=(x(1)-c(1,k))**2 + (x(2)-c(2,k))**2
     $             + (x(3)-c(3,k))**2
               rr(k)=sqrt(dist)
 211        continue 
            do 212 k=1,natoms
               do 212 l=1,natoms
                  if (k.ne.l) then
                     amu(k,l)=(rr(k)-rr(l))/rnuc(k,l)
                     if(adjust) then
                        amu(k,l)=amu(k,l)
     $                          +akl(k,l)*(one-amu(k,l)*amu(k,l))
                     endif
                  end if
 212        continue
c
c           --- generate the cell functions
            do 220 k=1,natoms
               pwtx(k)=one
               do 230 l=1,natoms
                  if(k.ne.l) then
                     xam=amu(k,l)
                     do 233 kk=1,kparm
                        yam=onehf*xam-half*xam**3
                        xam=yam
 233                 continue 
                     yam=half*(one-yam)
                     pwtx(k)=pwtx(k)*yam
                  end if
 230           continue
 220        continue
c
c           --- calculate denominator
c               wait until end to normalize cell functions
            denom=zero
            do 213 k=1,natoms
               denom=denom+pwtx(k)
 213        continue 
c
c           start gradient part
            do 300 k=1,natoms
               do 301 ix=1,3
                  term1(ix)=zero
                  term2(ix)=zero
                  term3(ix)=zero
                  gradwt(igrid,ix,k)=zero
  301          continue
               if(k.ne.atom) then
c                 term1 -gradb(pa)/denom
c                 grad(B) P(A) = -P(A) t(mu(AB)) gradB(mu(BA))
                  tmu=t(amu(atom,k))
                  do 303 ix=1,3
c
c                    need grad w.r.t. atom k of mu(k,atom)
c
                     gmu=gradmu(x(ix),c(ix,k),c(ix,atom),rr(k),       
     $                          rr(atom),rnuc(k,atom),
     $                          adjust,akl(k,atom))
                     term2(ix) = tmu*gmu*pwtx(atom)
                     term1(ix) = -term2(ix)/denom
  303             continue
c
c                 term2  sum(c.ne.b) gradb(pc) * pa/denom**2
c                 grad(B)P(c) = -Pc*t(mu(CB))*grad(B)(mu(BC))
c                 but since we want to subtract this term off later,
c                 drop the minus and just add later.
c                 Oh, and we already stuck the gradB pA term upstairs
c                 so don't bother recalculating it.
                  do 310 l=1,natoms
                     if(l.ne.atom.and.l.ne.k)then
                        tmu=t(amu(l,k))
                        do 311 ix=1,3
                           gmu=gradmu(x(ix),c(ix,k),c(ix,l),
     $                                rr(k),rr(l),rnuc(k,l),
     $                                adjust,akl(k,l))
                           term2(ix)=term2(ix)+tmu*gmu*pwtx(l)
  311                   continue
                     endif
  310             continue
                  do 312 ix=1,3
                     term2(ix) = term2(ix)*pwtx(atom)/denom**2
  312             continue
               
c                 term3  -gradb(pb)*pa/denom**2
c                 gradB P(B) = P(B) sum(c<>B) t(mu(BC)) gradB mu(BC)
c
                  do 320 l=1,natoms
                     if(l.ne.k)then
                        tmu=t(amu(k,l))
                        do 321 ix=1,3
                           gmu=gradmu(x(ix),c(ix,k),c(ix,l),
     $                                rr(k),rr(l),rnuc(k,l),
     $                                adjust,akl(k,l))
                           term3(ix)=term3(ix)+tmu*gmu
  321                   continue
                     endif
  320             continue
                  do 322 ix=1,3
                     term3(ix) = -term3(ix)*pwtx(atom)*pwtx(k)/denom**2
  322             continue
c
c                 add all parts together
c
                  do 330 ix=1,3
                     gradwt(igrid,ix,k)=term1(ix)+term2(ix)+term3(ix)
  330             continue
               endif
  300       continue
c
c           now have grad of weight except self-term
c           use translational invariance
c
            do 340 k=1,natoms
               if (k .eq. atom) goto 340
               do 341 ix=1,3
                  gradwt(igrid,ix,atom)=gradwt(igrid,ix,atom)
     $                                 - gradwt(igrid,ix,k)
  341          continue
  340       continue
c
c           normalize weights
c
            vwts(igrid)=pwtx(atom)/denom
 210     continue
c
c        --- multiply voronoi weights and numerical weights
         do 240 igrid=1,ngrid
            do 242 k=1,natoms
               do 241 ix=1,3
                  gradwt(igrid,ix,k)=gradwt(igrid,ix,k)*wts(igrid)
 241           continue
 242        continue 
            wts(igrid)=wts(igrid)*vwts(igrid)
 240     continue
      endif
c
c
      return
      end
      function t(amu)
      implicit none
      real*8 t,amu,s,p(3),xam,one,onehf,half,zero
      integer kparm, k
      parameter (kparm=3,zero=0.0d0,one=1.0d0,onehf=1.5d0,half=0.5d0)
      xam=amu
      do 100 k=1,kparm
         p(k) = onehf*xam - half*xam**3
         xam=p(k)
  100 continue
      s=half*(one-p(3))
      t=-(27.d0/16.d0)*(one-p(2)**2)*(one-p(1)**2)*(one-amu**2)
c
c     --- ASSume that if s=0 then t will be damned near zero too!      
      if (s.eq.zero) then
         t=zero
      else
         t=t/s
      endif
c
c
      return
      end 
      function gradmu(x,xa,xb,ra,rb,rab,adjust,akl)
c     calculates grad(a) of mu(a,b)
c     x = coord of grid point
c     xa = coord of atom a
c     xb = coord of atom b
c     ra = distance from grid point to atom a
c     rb = distance from grid point to atom b
c     rab distance from a to b
c     grada(mu(ab)) = 1/rab*U(a) - (ra-rb)/rab**2 * U(ab)
c     where u(a) is unit vector from grid point to atom a
c     u(ab) is unit vector from a->b
c
c     if adjusted cell radii are used, nu(ab), then adjust=.true.
c     amu = mu(ab)
c     akl = scale factor -- see Becke paper.
c     grada(nu(ab))=grada(mu(ab)) -2*akl*mu(ab)*grada(mu(ab))
c
      real*8 gradmu
      real*8 x,xa,xb,ra,rb,rab,akl
      logical adjust
c     --- local
      real*8 one,two,amu
c
      parameter (one=1.0d0,two=2.0d0)
c
      gradmu=(xa-x)/(ra*rab) - (ra-rb)*(xa-xb)/rab**3
      if(adjust) then
         amu=(ra-rb)/rab
         gradmu=(one-two*akl*amu)*gradmu
      endif
c
c
      return
      end
