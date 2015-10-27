*deck @(#)trn1dm.f	5.1  11/6/94
      subroutine trn1dm(ao1dm,c,nbf,nnp,t1,t2,mo1dm,nactiv,nnpact,
     $     ncore)
c
c***purpose: transformation of one-electron integrals.
c
c paul saxe             21 august 1984                   lanl
c
      implicit integer (a-z)
c
      real*8 ao1dm(nnp)
      real*8 c(nbf,nbf)
      real*8 t1(nbf,nbf)
      real*8 t2(nbf,nbf)
      real*8 mo1dm(nnpact)
c
      call trtosq(t2,mo1dm,nactiv,nnpact)
      call ebc(t1,c(1,ncore+1),t2,nbf,nactiv,nactiv)
      call ebct(t2,t1,c(1,ncore+1),nbf,nactiv,nbf)
      call sqtotr(ao1dm,t2,nbf,nnp)
c
c
      return
      end
