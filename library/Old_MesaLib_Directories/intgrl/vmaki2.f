*deck @(#)vmaki2.f	5.1  11/6/94
      subroutine vmaki2(i2,h,c,nat,symcen,nmax,mmax,imax,jmax,kmax,
     #                  lmax,nv)
c
c     module to form the four-centre two-dimensional integral
c                              paul saxe   23 april 1984  brl
c
c     vector version
c  paul saxe                      13 july 1984              lanl
c
      implicit integer (a-z)
c
      real*8 i2(nv,3,0:imax,0:jmax,0:mmax,0:lmax)
      real*8 h(nv,3,0:nmax,0:jmax,0:mmax)
      real*8 c(3,nat),scalar
      integer imax,jmax,kmax,lmax,mmax,symcen(4)
c
c     ----- transfer elements from h to i -----
c
      do 4 m=0,mmax
         do 3 j=0,jmax
            do 2 i=0,imax
               do 1 coord=1,3
                  call vmove(i2(1,coord,i,j,m,0),h(1,coord,i,j,m),nv)
    1          continue
    2       continue
    3    continue
    4 continue
c
c     ----- i(i,j,k,l)=i(i,j,k+1,l-1)+ckmcl*i(i,j,k,l-1) -----
c
      do 9 l=1,lmax
         do 8 k=0,mmax-l
            do 7 j=0,jmax
               do 6 i=0,imax
                  do 5 coord=1,3
                     scalar=c(coord,symcen(3))-c(coord,symcen(4))
                     call vwxs(i2(1,coord,i,j,k,l),
     #                         i2(1,coord,i,j,k+1,l-1),
     #                         i2(1,coord,i,j,k,l-1),
     #                         scalar,1,nv)
    5             continue
    6          continue
    7       continue
    8    continue
    9 continue
c
c
      return
      end
