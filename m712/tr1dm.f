*deck @(#)tr1dm.f	5.1  11/6/94
      subroutine tr1dm(t1,t2,dao,prtoao,dpr,nbf,nnp,npf,nnprim,ndmat)
c
c***begin prologue     tr1dm
c***date written       860820  (yymmdd)
c***revision date      yymmdd  (yymmdd)
c***keywords
c***author             saxe, paul (lanl)
c***source             @(#)tr1dm.f	5.1   11/6/94
c***purpose
c***description
c
c***references
c***routines called
c***end prologue       tr1dm
c
      implicit integer (a-z)
c
      real*8 t1(npf,npf),t2(npf,npf),dao(nnp,ndmat),prtoao(npf,nbf)
      real*8 dpr(nnprim,ndmat)
c
c     ----- transform all the one-particle density matrices to the primitive
c           basis
c
      do 1 dmat=1,ndmat
         call trtosq(t1,dao(1,dmat),nbf,nnp)
         call ebct(t2,t1,prtoao,nbf,nbf,npf)
         call ebc(t1,prtoao,t2,npf,nbf,npf)
         call sqtotr(dpr(1,dmat),t1,npf,nnprim)
    1 continue
c
c
      return
      end
