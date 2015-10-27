*deck @(#)facab.f	5.1  11/6/94
      subroutine facab (na,nb,cax,cbx,xab)
      implicit real*8(a-h,o-z)
c
      integer ind(6)
      real*8 xab(1)
      real*8 binom(21)
      data binom /1., 1.,1., 1.,2.,1., 1.,3.,3.,1., 1.,4.,6.,4.,1.,
     $            1.,5.,10.,10.,5.,1./
      data ind /1,2,4,7,11,16/
c
      save binom,ind
c
      parameter (zero=0.0d+00, one=1.0d+00)
c     
      iend=na+nb+1
      do 10 i=1,iend
         xab(i)=zero
   10 continue
      na1=na+1
      nb1=nb+1
      naind=ind(na1)
      nbind=ind(nb1)
      do 60 ia1=1,na1
         ia=ia1-1
         if ((na-ia).eq.0) then
            caxp=one
         else 
            caxp=cax**(na-ia)
         endif
         afac=binom(naind+ia)*caxp
         do 55 ib1=1,nb1
            ib=ib1-1
            if ((nb-ib).eq.0) then
               cbxp=one
            else 
               cbxp=cbx**(nb-ib)
            endif
            xab(ia+ib+1)=xab(ia+ib+1)+afac*binom(nbind+ib)*cbxp
   55    continue
   60 continue
c
c
      return
      end
