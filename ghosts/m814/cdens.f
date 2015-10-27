*deck @(#)cdens.f	1.1  11/30/90
      subroutine cdens(corden,nnp,c,nbf,ncore,t1)
c
      implicit integer (a-z)
c
      real*8 corden(nnp),c(nbf,ncore),t1(nbf,nbf)
c
      call ebct(t1,c,c,nbf,ncore,nbf)
      call sqtotr(corden,t1,nbf,nnp)
      do 1 i=1,nnp
         corden(i)=corden(i)*2.0d+00
    1 continue
c
      return
      end
