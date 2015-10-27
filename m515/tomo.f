*deck @(#)tomo.f	5.1  11/28/95
      subroutine tomo(nbf,nnp,nshell,ncoul,nexch,hao,h,jmat,kmat,
     $                c,t1,t2)
c***begin prologue     tomo.f
c***date written       930610  
c***revision date      11/6/94      
c
c***keywords           
c***author             martin,richard and saxe,paul (lanl) 
c***source             @(#)tomo.f	5.1   11/28/95
c***purpose            to transform the one-electron,j,and k matrices
c                      to a molecular orbital representation. 
c***description
c     this routine transforms the one-electron, coulomb and
c     exchange matrices to a molecular orbital representation.
c
c     
c***references
c
c***routines called
c
c***end prologue       tomo.f
      implicit none
c     --- input variables -----
      integer nbf,nnp,nshell,ncoul,nexch
c     --- input arrays (unmodified) ---
      real*8 hao(nnp)
      real*8 c(nbf,nbf)
c     --- input arrays (scratch) ---
c     --- output arrays ---
      real*8 h(nnp),jmat(nnp,nshell),kmat(nnp,nshell)
c     --- output variables ---
c     --- scratch arrays ---
      real*8 t1(nbf,nbf),t2(nbf,nbf)
c     --- local variables ---
      integer i
      integer inp,iout
c
      common /io/ inp,iout
c
c     --- transform the one-electron integrals to the mo basis ---
      call trtosq(t2,hao,nbf,nnp)
      call ebc(t1,t2,c,nbf,nbf,nbf)
      call ebtc(t2,c,t1,nbf,nbf,nbf)
      call sqtotr(h,t2,nbf,nnp)
c
c     --- transform the j and k matrices to the mo basis ---
      do 10 i=1,ncoul
         call trtosq(t2,jmat(1,i),nbf,nnp)
         call ebc(t1,t2,c,nbf,nbf,nbf)
         call ebtc(t2,c,t1,nbf,nbf,nbf)
         call sqtotr(jmat(1,i),t2,nbf,nnp)
 10   continue
      do 20 i=1,nexch
         call trtosq(t2,kmat(1,i),nbf,nnp)
         call ebc(t1,t2,c,nbf,nbf,nbf)
         call ebtc(t2,c,t1,nbf,nbf,nbf)
         call sqtotr(kmat(1,i),t2,nbf,nnp)
 20   continue
c
c
      return
      end
