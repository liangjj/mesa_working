*deck @(#)trn1e.f	5.1  11/6/94
      subroutine trn1e(h,c,nbf,norbs,nnp,numij,t1,t2,hmo)
c
c***purpose: transformation of one-electron integrals.
c
c paul saxe             21 august 1984                   lanl
c
      implicit integer (a-z)
c
      real*8 h(nnp),c(nbf,norbs),t1(norbs,nbf),t2(norbs,norbs)
      real*8 hmo(numij)
c
      call trtosq(t2,h,nbf,nnp)
      call ebtc(t1,c,t2,norbs,nbf,nbf)
      call ebc(t2,t1,c,norbs,nbf,norbs)
      call sqtotr(hmo,t2,norbs,numij)
c
c
      return
      end
