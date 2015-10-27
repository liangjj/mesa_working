*deck %W%  %G%
      subroutine facab (na,nb,cax,cbx,xab)
      implicit integer(a-z)
c
c     ----- arguments unchanged -----
      integer na,nb
      real*8 cax,cbx
c     ------ arguments returned -----
      real*8 xab(0:na+nb)
c
c     ----- local variables -----
      real*8 zero,one
      real*8 caxp,cbxp,afac
c     ----- external function -----
      real*8 bico
c
      parameter (zero=0.0d+00, one=1.0d+00)
c     
      do 10 i=0,na+nb
         xab(i)=zero
   10 continue
c
c****** rewrite to avoid exponential
      do 60 ia=0,na
         if ((na-ia).eq.0) then
            caxp=one
         else 
            caxp=cax**(na-ia)
         endif
         afac=bico(na,ia)*caxp
         do 55 ib=0,nb
            if ((nb-ib).eq.0) then
               cbxp=one
            else 
               cbxp=cbx**(nb-ib)
            endif
            xab(ia+ib)=xab(ia+ib)+afac*bico(nb,ib)*cbxp
   55    continue
   60 continue
c
c
      return
      end
