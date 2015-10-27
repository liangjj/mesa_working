*deck @(#)rsolver.f	5.1 11/6/94
      subroutine rsolver(nrings,nrtot,nr,lmax,nlm,ptlm,
     $                   rpts,wt,nwtot,flm,j,y,scr,
     $                   tmp,psilm,int0f,int0b)
c***begin prologue     rsolver.f
c***date written       940304    (yymmdd)  
c***revision date      11/6/94      
c
c***keywords           
c***author             martin, richard(lanl) 
c***source             @(#)rsolver.f	5.1   11/6/94
c***purpose            solves the radial equation for the Poisson problem
c                                       [rho(r2)]
c                         U(r1)=r1* Int ---------
c                                        |r1-r2|
c        
c                      which results from the Laplace expansion of 1/r12.
c
c                      Int(0,Inf) [4*pi*r*f(r)* (r<)**l+1 ]
c                                                ---------
c                                              (2l+1)*(r>)**l
c***description       
c                      
c
c***references         
c                      Newton-Cotes quadrature for the radial equation:
c                      D.H. Oza and J. Callaway, J. Comp.Phys 68, 89(1987).
c                      B. Basden and R.R. Lucchese, J. Comp. Phys. 77,524(1988).
c
c***routines called
c
c***end prologue       rsolver.f
      implicit none
c     --- input variables -----
      integer nrings,nrtot,nwtot,lmax,nlm
c     --- input arrays (unmodified) ---
      integer nr(nrings),ptlm(0:lmax,-lmax:lmax)
      real*8 rpts(nrtot),wt(nwtot)
      real*8 flm(nrtot,nlm)
c     --- input arrays (scratch) ---
      real*8 j(nrtot,0:lmax),y(nrtot,0:lmax)
      real*8 scr(nrtot),tmp(nrtot)
c     --- output arrays ---
      real*8 int0f(nlm),int0b(nlm)
      real*8 psilm(nrtot,nlm)
c     --- output variables ---
c     --- scratch arrays ---
c     --- local variables ---
      integer l,m,lm,ptr,ptwt,ns,i
      integer inp,iout
      real*8 zero,one,two,four,pi,atan
      real*8 fourpi,fact
      logical debug
c
      parameter (debug=.false.)
c
      data zero,one,two,four/0.0d+00,1.0d+00,2.0d+00,4.0d+00/
      save zero,one,four
c
      common/io/inp,iout
c
      pi=four*atan(one)
      fourpi=four*pi
c
c     --- get the homogeneous solutions which define the greens's function
      call laplace(rpts,j,y,nrtot,lmax)
c
c     --- for each lm component, solve the radial equation by splitting
c         the singularity and finding the forward and backward integrals.
c         the density, decomposed into ylm components, enters in flm,
c         the solution is returned in psilm.
      lm=0
      do 100 l=0,lmax
         do 90 m=-l,l
            lm=lm+1
c           --- scale density by 4*pi*r/2l+1
            call vmul(flm(1,ptlm(l,m)),flm(1,ptlm(l,m)),rpts,nrtot)
            fact=fourpi/float((2*l+1))
            call sscal(nrtot,fact,flm(1,ptlm(l,m)),1)
c
            ptr=1
            ptwt=1
c           --- initialize forward integral and u at the origin.
            int0f(lm)=zero
            psilm(1,ptlm(l,m))=int0f(lm)
            do 60 ns=1,nrings
               call forwrd(flm(ptr,ptlm(l,m)),j(ptr,l),y(ptr,l),
     $                     wt(ptwt),int0f(lm),scr(ptr),nr(ns))
c              --- see routine newton for an explanation of the funny
c                  business regarding these pointers.
               ptr=ptr+nr(ns)-1
               ptwt=ptwt+nr(ns)*(nr(ns)-1)
   60       continue
c           --- accumulate the forward piece of u. 
            call vmove(psilm(2,ptlm(l,m)),scr,nrtot-1)
c
c           --- now do the backward integration. note the pointers
c               go backward over the rings
            int0b(lm)=zero
            scr(nrtot)=int0b(lm)
            do 70 ns=nrings,1,-1
               ptr=ptr-nr(ns)+1
               ptwt=ptwt-nr(ns)*(nr(ns)-1)
               call bakwrd(flm(ptr,ptlm(l,m)),j(ptr,l),y(ptr,l),
     $                    wt(ptwt),int0b(lm),scr(ptr),nr(ns))
   70       continue
c
c           --- combine the forward and backward integrals.
            call vadd(psilm(1,ptlm(l,m)),psilm(1,ptlm(l,m)),
     $                scr,nrtot)
            if(debug) then
               write(iout,*) 'l,m',l,m,nrtot
               write(iout,*) '   forward sum',
     $                       int0f(lm)/sqrt(four*pi*(l+l+1))
               write(iout,*) '   backward sum',sqrt(four*pi)*int0b(lm)
               write(iout,*) 'ulm',(psilm(i,ptlm(l,m)),i=1,nrtot)
            endif
   90    continue
  100 continue
c
c
      return
      end
