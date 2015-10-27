*deck @(#)vmakh.f	5.1  11/6/94
      subroutine vmakh(h,g,c,nat,symcen,nmax,mmax,imax,jmax,nv)
c
c     module to form the three-centre two-dimensional integrals
c                               paul saxe   23 april 1984   brl
c
c       vector version.
c  paul saxe                         13 july 1984                   lanl
c
      implicit integer (a-z)
c
      real*8 h(nv,0:nmax,0:jmax,0:mmax,3)
      real*8 g(nv,0:nmax,0:mmax,3)
      real*8 c(3,nat)
      real*8 scalar
      integer i,j,n,m,coord
      integer imax,jmax,nmax,mmax,nv,symcen(4)
c
c     ----- transfer elements from g to h -----
c
      do 3 m=0,mmax
         do 2 n=0,nmax
            do 1 coord=1,3
               call vmove(h(1,n,0,m,coord),g(1,n,m,coord),nv)
    1       continue
    2    continue
    3 continue
c
c     ----- h(i,j,m)=h(i+1,j-1,m)+cimcj*h(i,j-1,m) -----
c
      do 7 m=0,mmax
         do 6 j=1,jmax
            do 5 i=0,nmax-j
               do 4 coord=1,3
                  scalar=c(coord,symcen(1))-c(coord,symcen(2))
                  call vwxs(h(1,i,j,m,coord),h(1,i+1,j-1,m,coord),
     #                      h(1,i,j-1,m,coord),scalar,1,nv)
    4          continue
    5       continue
    6    continue
    7 continue
c
c
      return
      end
