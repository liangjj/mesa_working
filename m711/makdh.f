*deck @(#)makdh.f	5.1  11/6/94
      subroutine makdh(h,dh,dg,c,nat,symcen,nmax,mmax,imax,jmax,nv,
     $                 ncentr,npass)
c***begin prologue     makdh
c***date written       851106   (yymmdd)
c***revision date      11/6/94
c***keywords           derivative integrals, rys polynomials,
c                      derivative three-center integrals
c
c***author             saxe, paul,    (lanl)
c***source             @(#)makdh.f	5.1   11/6/94
c***purpose            to form the derivatives of the three-center,
c                      two-dimensional integrals h.
c
c***description         
c
c
c***references         (none)
c
c***routines called    (none)
c***end prologue       makdh.f
      implicit none
c     --- input variables -----
      integer nat,imax,jmax,nmax,mmax,nv,ncentr,npass
c     --- input arrays (unmodified) ---
      integer symcen(4)
      real*8 c(3,nat)
      real*8 h(nv,3,0:nmax,0:jmax,0:mmax)
      real*8 dg(nv,3,0:nmax,0:mmax,ncentr)
c     --- input arrays (scratch) ---
c     --- output arrays ---
      real*8 dh(nv,3,0:nmax,0:jmax,0:mmax,ncentr)
c     --- output variables ---
c     --- scratch arrays ---
c     --- local variables ---
      integer centr,m,n,i,j,coord,ii
      real*8 scalar
c
c     --- transfer elements from g to h ---
      do 3 centr=1,ncentr
         do 2 m=0,mmax
            do 1 n=0,nmax
               call vmove(dh(1,1,n,0,m,centr),dg(1,1,n,m,centr),nv*3)
    1       continue
    2    continue
    3 continue
c
c     --- h(i,j,m)=h(i+1,j-1,m)+cimcj*h(i,j-1,m) ---
      do 17 m=0,mmax
         do 16 j=1,jmax
            do 15 i=0,nmax-j
               do 14 coord=1,3
                  scalar=c(coord,symcen(1))-c(coord,symcen(2))
                  do 10 ii=1,nv
                     dh(ii,coord,i,j,m,1)=dh(ii,coord,i+1,j-1,m,1)+
     $                            scalar*dh(ii,coord,i,j-1,m,1)+
     $                            h(ii,coord,i,j-1,m)
   10             continue
                  if (npass.le.2) then
                     do 11 ii=1,nv
                        dh(ii,coord,i,j,m,2)=dh(ii,coord,i+1,j-1,m,2)+
     $                               scalar*dh(ii,coord,i,j-1,m,2)-
     $                               h(ii,coord,i,j-1,m)
   11                continue
                  else if (npass.eq.3) then
                     do 12 ii=1,nv
                        dh(ii,coord,i,j,m,2)=dh(ii,coord,i+1,j-1,m,2)+
     $                               scalar*dh(ii,coord,i,j-1,m,2)
   12                continue
                  end if
                  if (npass.eq.1) then
                     do 13 ii=1,nv
                        dh(ii,coord,i,j,m,3)=dh(ii,coord,i+1,j-1,m,3)+
     $                               scalar*dh(ii,coord,i,j-1,m,3)
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
