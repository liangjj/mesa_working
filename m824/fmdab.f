*deck @(#)fmdab.f	5.1  11/6/94
      subroutine fmdab(d,c,t1,t2,nactiv,nnpact,nbf,nnp,ncore)
c
      implicit integer (a-z)
c
      real*8 d(nnp),c(nbf,nbf),t1(nbf,nbf),t2(nbf,nbf)
      real*8 trace
c
      call iosys('read real mcscf_mo_1pdm from rwf',nnpact,d,
     #            0,' ')
c
c     ----- half off-diagonals to return to true density matrix
c
      ij=0
      do 2 i=1,nactiv
         do 1 j=1,i-1
            ij=ij+1
            d(ij)=d(ij)*0.5d+00
    1    continue
         ij=ij+1
    2 continue
      do 33 ij=1,nnpact
         d(ij)=d(ij)*2.0d+00
   33 continue
c
cps      call print(d,nnpact,nactiv,6)
c
      trace=0.0d+00
      ii=0
      do 3 i=1,nactiv
         ii=ii+i
         trace=trace+d(ii)
    3 continue
c
cps      write (6,4) trace
cps    4 format(' trace=',f12.6)
c
c
      call trtosq(t1,d,nactiv,nnpact)
      call ebc (t2,c(1,ncore+1),t1,nbf,nactiv,nactiv)
      call ebct(t1,t2,c(1,ncore+1),nbf,nactiv,nbf)
      call sqtotr(d,t1,nbf,nnp)
c
c
      return
      end
