*deck @(#)derexp.f	1.1  11/30/90
      subroutine derexp(i2,di,imax,jmax,kmax,lmax,mmax,nv,ndcen,d1exp)
c
      implicit integer (a-z)
c
      real*8 i2(nv,0:imax,0:jmax,0:mmax,0:lmax,3)
      real*8 di(nv,0:imax,0:jmax,0:mmax,0:lmax,3,ndcen)
      real*8 d1exp(nv,3,ndcen)
c
c     ----- add derivatives of exponential terms -----
c
      do 7 i=0,imax
         do 6 j=0,jmax
            do 5 k=0,kmax
               do 4 l=0,lmax
                  do 3 coord=1,3
                     do 2 centre=1,ndcen
                        do 1 m=1,nv
                           di(m,i,j,k,l,coord,centre)=
     #                              di(m,i,j,k,l,coord,centre)+
     #                       d1exp(m,coord,centre)*i2(m,i,j,k,l,coord)
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
