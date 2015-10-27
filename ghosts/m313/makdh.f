*deck @(#)makdh.f	1.1  11/30/90
      subroutine makdh(h,dh,dg,c,nat,symcen,nmax,mmax,imax,jmax,nv,
     #                 ncentr,npass)
c
c***begin prologue     makdh
c***date written       851106   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c***keywords           derivative integrals, rys polynomials,
c                      derivative three-center integrals
c
c***author             saxe, paul,    (lanl)
c***source             @(#)makdh.f	1.1   11/30/90
c***purpose            to form the derivatives of the three-center,
c                      two-dimensional integrals h.
c
c***description        #
c
c
c***references         (none)
c
c***routines called    (none)
c
c   common blocks      (none)
c
c***end prologue       makdh
c
c
      implicit integer (a-z)
c
      real*8 h(nv,0:nmax,0:jmax,0:mmax,3)
      real*8 dh(nv,0:nmax,0:jmax,0:mmax,3,ncentr)
      real*8 dg(nv,0:nmax,0:mmax,3,ncentr)
      real*8 c(3,nat)
      real*8 scalar
      integer i,j,n,m,coord
      integer imax,jmax,nmax,mmax,nv,symcen(4)
c
c     ----- transfer elements from g to h -----
c
      do 3 centr=1,ncentr
      do 2 m=0,mmax
         do 1 n=0,nmax
            call vmove(dh(1,n,0,m,1,centr),dg(1,n,m,1,centr),nv)
            call vmove(dh(1,n,0,m,2,centr),dg(1,n,m,2,centr),nv)
            call vmove(dh(1,n,0,m,3,centr),dg(1,n,m,3,centr),nv)
    1    continue
    2 continue
    3 continue
c
c     ----- h(i,j,m)=h(i+1,j-1,m)+cimcj*h(i,j-1,m) -----
c
      do 17 m=0,mmax
         do 16 j=1,jmax
            do 15 i=0,nmax-j
               do 14 coord=1,3
                  scalar=c(coord,symcen(1))-c(coord,symcen(2))
                  do 10 ii=1,nv
                     dh(ii,i,j,m,coord,1)=dh(ii,i+1,j-1,m,coord,1)+
     #                            scalar*dh(ii,i,j-1,m,coord,1)+
     #                            h(ii,i,j-1,m,coord)
   10             continue
                  if (npass.le.2) then
                     do 11 ii=1,nv
                        dh(ii,i,j,m,coord,2)=dh(ii,i+1,j-1,m,coord,2)+
     #                               scalar*dh(ii,i,j-1,m,coord,2)-
     #                               h(ii,i,j-1,m,coord)
   11                continue
                  else if (npass.eq.3) then
                     do 12 ii=1,nv
                        dh(ii,i,j,m,coord,2)=dh(ii,i+1,j-1,m,coord,2)+
     #                               scalar*dh(ii,i,j-1,m,coord,2)
   12                continue
                  end if
                  if (npass.eq.1) then
                     do 13 ii=1,nv
                        dh(ii,i,j,m,coord,3)=dh(ii,i+1,j-1,m,coord,3)+
     #                               scalar*dh(ii,i,j-1,m,coord,3)
   13                continue
                  end if
   14          continue
   15       continue
   16    continue
   17 continue
c
c
      return
      end
