*deck %W%  %G%
      subroutine soaot(vec,tr,scr,nao,nso)
      implicit real*8(a-h,o-z)
      dimension vec(nso,nao),tr(nao,nso),scr(nao)
      if(nso.ne.nao)stop ' quit in soaot '
      do 100 i=1,nso
      do 110 j=1,nao
      x=0.0d0
      do 120 k=1,nso
  120 x=x+vec(k,i)*tr(j,k)
  110 scr(j)=x
      do 111 j=1,nao
  111 vec(j,i)=scr(j)
  100 continue
      return
      end
