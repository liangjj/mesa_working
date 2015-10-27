*deck @(#)tr1dm.f	1.1  4/25/95
      subroutine tr1dm(t1,t2,dao,prtoao,dpr,nbf,nnp,npf,nnprim)
c***begin prologue     tr1dm.f
c***date written       860820  (yymmdd)  
c***revision date      11/6/94      
c
c***keywords           
c***author             saxe, paul  
c***source             @(#)tr1dm.f	1.1   4/25/95
c***purpose            to transform the hf density matrices from the ao
c                      to the primitive basis
c***description
c       prtoao   ...   the primitive to contracted basis transformation matrix
c       dao      ...   the contracted basis density matrix.
c       t1       ...   scratch (npf*npf).
c       t2       ...   scratch (nbf*npf).
c       nbf      ...   number of contracted basis functions.
c       npf      ...   number of primitive basis functions.
c       nnp      ...   nbf*(nbf+1)/2
c       nnprim   ...   npf*(npf+1)/2  
c
c       dpr      ...   the primitive basis density matrix.
c
c***references
c
c***routines called
c
c***end prologue       tr1dm.f
      implicit none
c     --- input variables -----
      integer nbf,nnp,npf,nnprim
c     --- input arrays (unmodified) ---
      real*8 prtoao(npf,nbf)
      real*8 dao(nnp)
c     --- input arrays (scratch) ---
      real*8 t1(npf*npf)
c     t1 must hold (npf*npf)
      real*8 t2(nbf*npf)
c     --- output arrays ---
      real*8 dpr(nnprim)
c     --- output variables ---
c     --- scratch arrays ---
c     --- local variables ---
      logical debug
      integer inp,iout
c
      parameter (debug=.false.)
c
 1000 format(5x,'the primitive-to-ao transformation matrix')
 1010 format (5x,'the ao hf density matrix ')
 1020 format (5x,'the primitive hf density matrix:',i2)
c
      common /io/     inp,iout
c
      if (debug) then
         write (iout,1000)
         call matout(prtoao,npf,nbf,npf,nbf,iout)
         write (iout,1010) 
         call print(dao,nnp,nbf,iout)
      end if
c
c     --- transform all the one-particle density matrices to the
c         primitive basis
      call trtosq(t1,dao,nbf,nnp)
      call ebct(t2,t1,prtoao,nbf,nbf,npf)
      call ebc(t1,prtoao,t2,npf,nbf,npf)
      call sqtotr(dpr,t1,npf,nnprim)
c
      if (debug) then
         write (iout,1020) 
         call print(dpr,nnprim,npf,iout)
      end if
c
c
      return
      end
