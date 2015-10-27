*deck @(#)maketa.f	1.1  11/30/90
      subroutine maketa(tanco,tanao,ta,nco,nao,nob)
      implicit integer(a-z)
      real*8 tanco(nob,nco),tanao(nob,nao),ta(nob,nob)
c
      if(nco.eq.0) go to 100
c
      do 1 i=1,nco
         do 2 j=1,nob
            tanco(j,i)=0.0d+00
  2      continue
 1    continue
c
      do 3 j=1,nco
         do 4 i=1,j
            tanco(i,j)=-ta(i,j)
  4      continue
 3    continue
c
      do 5 i=1,nco
         tanco(i,i)=(0.5d+00)*tanco(i,i)
  5   continue
c
  100 continue
c
      if(nao.eq.0) go to 200
c
      do 11 i=1,nao
         do 12 j=1,nob
            tanao(j,i)=0.0d+00
  12     continue
 11   continue
c
      jend=nco
      do 13 i=1,nao
         jend=jend+1
         do 14 j=1,jend
            tanao(j,i)=-ta(j,jend)
  14     continue
 13   continue
c
      do 15 i=1,nao
       tanao(i+nco,i)=tanao(i+nco,i)*.5d+00
  15  continue
c
 200  continue
c
      return
      end
