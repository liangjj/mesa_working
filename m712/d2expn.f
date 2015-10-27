*deck @(#)d2expn.f	5.1  11/6/94
      subroutine d2expn(i2,di,d2i,imax,jmax,kmax,lmax,mmax,nv,ndcen,
     $     nd2,d1exp,d2exp,cen)
c
      implicit integer (a-z)
c
      real*8 i2(nv,3,0:imax,0:jmax,0:mmax,0:lmax)
      real*8 di(nv,3,0:imax,0:jmax,0:mmax,0:lmax,ndcen)
      real*8 d2i(nv,3,0:imax,0:jmax,0:mmax,0:lmax,nd2)
      real*8 d1exp(nv,3,ndcen)
      real*8 d2exp(nv,nd2)
      integer cen(2,nd2)
c
c     ----- add derivatives of exponential terms -----
c
      do 7 i=0,imax
         do 6 j=0,jmax
            do 5 k=0,kmax
               do 4 l=0,lmax
                  do 3 coord=1,3
                     do 2 d2=1,nd2
                        c1=cen(1,d2)
                        c2=cen(2,d2)
                        do 1 m=1,nv
                           d2i(m,coord,i,j,k,l,d2)=
     #                          d2i(m,coord,i,j,k,l,d2)+
     #                       d1exp(m,coord,c1)*di(m,coord,i,j,k,l,c2)+
     #                       d1exp(m,coord,c2)*di(m,coord,i,j,k,l,c1)+
     $                          (d1exp(m,coord,c1)*d1exp(m,coord,c2)+
     $                          d2exp(m,d2))*i2(m,coord,i,j,k,l)
    1                   continue
    2                continue
    3             continue
    4          continue
    5       continue
    6    continue
    7 continue
c
c
      return
      end
