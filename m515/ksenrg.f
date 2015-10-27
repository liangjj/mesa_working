*deck @(#)ksenrg.f	5.1  11/28/95
      subroutine ksenrg(jmat,exc,kmat,dmat,nnp,nbf,ncoul,nexch,ndmat,h,
     $                  enuc,energy,calc,eonel,etwoel,
     $                  nshell,fcoef,alpha,beta,hf,ehf,khf)
c***begin prologue     ksenrg.f
c***date written       930501 
c***revision date      11/6/94      
c
c***keywords           dft, energy
c***author             martin, richard (lanl) 
c***source             @(#)ksenrg.f	5.1   11/28/95
c***purpose            computes the kohn-sham total energy 
c***description
c     
c    
c
c***references
c
c***routines called    sdot
c
c***end prologue       ksenrg.f
      implicit none
c     --- input variables -----
      integer nbf,nnp,ncoul,nexch,ndmat
      integer nshell
      logical hf
      real*8 exc,ehf
      character*(*) calc
c     --- input arrays (unmodified) ---
      real*8 jmat(nnp,ncoul),kmat(nnp,nexch),dmat(nnp,ndmat),
     $     khf(nnp,nexch)
      real*8 h(nnp)
      real*8 fcoef(nshell),alpha(nshell,nshell),beta(nshell,nshell)
c     --- input arrays (scratch) ---
c     --- output arrays ---
c     --- output variables ---
      real*8 enuc,energy,eonel,etwoel
c     --- scratch arrays ---
c     --- local variables ---
      integer i,j,jj
      real*8 two,half,four,one,three
      real*8 sdot
c
      parameter (two=2.0d+00, half=1.0d+00/two,four=4.0d+00)
      parameter (one=1.0d+00 ,three=3.0d+00)
c
c     --- double the off-diagonal part of the density matrices
c           to account for symmetry
      do 30 i=1,ndmat
         do 10 j=1,nnp
            dmat(j,i)=dmat(j,i)*two
   10    continue
         do 20 j=1,nbf
            jj=j*(j+1)/2
            dmat(jj,i)=dmat(jj,i)*half
   20    continue
   30 continue
c    
c     these sections evaluate the contribution from the one-electron
c     and coulomb operators to the total energy.  the exchange contribution
c     (exc) is evaluated from the total and spin-densities elsewhere.
      if (calc.eq.'closed') then
c
c        --- closed shell scf ---
c            e= 2*h(ii) + 2*j(ii,jj) - k(ii,jj)
c
         eonel=two*sdot(nnp,h,1,dmat,1)
         etwoel=two*sdot(nnp,jmat,1,dmat,1)
         if (hf) ehf=-sdot(nnp,khf,1,dmat,1)
c
      else if (calc.eq.'open') then
c
c        --- high-spin open-shell scf ---
c
c        e = 2 * h(ii) +  2  * j(ii,jj) -       k(ii,jj)
c              + h(aa) + 1/2 * j(aa,bb) - 1/2 * k(aa,bb)
c              +          2  * j(ii,aa) -       k(ii,aa)
c
c          i,j = d.o. ; a,b = h.o
c
         eonel=two*sdot(nnp,h,1,dmat(1,1),1)+
     $             sdot(nnp,h,1,dmat(1,2),1)
         etwoel= two*sdot(nnp,jmat(1,1),1,dmat(1,1),1)
     $         +half*sdot(nnp,jmat(1,2),1,dmat(1,2),1)
     $          +two*sdot(nnp,jmat(1,1),1,dmat(1,2),1)
c
      else if (calc.eq.'general') then
         eonel=0.d0
         etwoel=0.d0
         do 40 i=1,nshell-1
            eonel=eonel+  two*fcoef(i)*sdot(nnp,h,1,dmat(1,i),1)
            etwoel=etwoel+alpha(i,i)*sdot(nnp,jmat(1,i),1,dmat(1,i),1)
  40     continue
         if(nshell.gt.2) then
            do 60 i=2,nshell-1
               do 50 j=1,i-1
                  etwoel= etwoel
     $                 +alpha(i,j)*sdot(nnp,jmat(1,i),1,dmat(1,j),1)
     $                 +alpha(j,i)*sdot(nnp,jmat(1,j),1,dmat(1,i),1)
  50           continue
  60        continue
         end if
      endif
      etwoel=etwoel+exc
      energy=enuc+eonel+etwoel
      if (hf) energy=energy+ehf
c
c     --- halve the off-diagonal part of the density matrices
c           to return them to normal
c
      do 6 i=1,ndmat
         do 4 j=1,nnp
            dmat(j,i)=dmat(j,i)*half
    4    continue
         do 5 j=1,nbf
            jj=j*(j+1)/2
            dmat(jj,i)=dmat(jj,i)*two
    5    continue
    6 continue
c
c
      return
      end
