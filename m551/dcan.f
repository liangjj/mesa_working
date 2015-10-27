*deck @(#)dcan.f	1.1  11/30/90
      subroutine dcan(h,g,gcan,norbs,nnp,nnp2)
c
      implicit integer (a-z)
c
      real*8 h(nnp),g(nnp,nnp),gcan(nnp2)
c
c
      ijkl=0
      do 9 i=1,norbs
         do 8 j=1,i
            ij=i*(i-1)/2+j
            do 7 k=1,i
               if (k.eq.i) then
                   lmax=j
               else
                   lmax=k
               end if
               do 6 l=1,lmax
                  kl=k*(k-1)/2+l
                  ijkl=ijkl+1
                  gcan(ijkl)=g(ij,kl)*8.0d+00
                  if (i.eq.j) gcan(ijkl)=gcan(ijkl)*0.5d+00
                  if (k.eq.l) gcan(ijkl)=gcan(ijkl)*0.5d+00
                  if (ij.eq.kl) gcan(ijkl)=gcan(ijkl)*0.5d+00
    6          continue
    7       continue
    8    continue
    9 continue
c
      return
      end
