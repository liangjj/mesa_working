*deck %W%  %G%
      subroutine mcbcku(temp,u,vec,nmo,nao,kk,id2h)
c
c***begin prologue
c***date written       871022   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c
c***keywords
c***author             lengsfield, byron (brl)
c***source             %W%   %G%
c
c***purpose
c
c***description
c
c***references
c
c***routines called    (none)
c
c***end prologue
c
      implicit real*8 (a-h,o-z)
cc
c     extended dummy u,vec,temp
cc
      dimension u(nmo,nmo),vec(nao,nmo),temp(2)
c
      do 110 i=1,nao
         do 100 j=1,nmo
            xx=0.0d0
            do 90 k=1,nmo
               xx=xx+vec(i,k)*u(j,k)
 90         continue
            temp(j)=xx
 100     continue
         do 105 j=1,nmo
            vec(i,j)=temp(j)
 105     continue
 110  continue
c
      if(id2h.eq.0)go to 160
c
      nao2=nao/2
      do 130 j=1,nmo
         do 120 i=1,nao2
            xx=(abs(vec(i,j))+abs(vec(i+nao2,j)))*.5d0
            vec(i,j)=sign(xx,vec(i,j))
            vec(i+nao2,j)=sign(xx,vec(i+nao2,j))
 120     continue
 130  continue
c
 160  continue
c
      return
      end
