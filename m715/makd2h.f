*deck @(#)makd2h.f	5.1  11/28/95
      subroutine makd2h(dh,d2h,d2g,c,nat,symcen,nmax,mmax,imax,jmax,nv,
     $                 ncentr,nd2,npass)
c***begin prologue     makd2h.f
c***date written       851106   (yymmdd)
c***revision date      11/6/94
c***keywords           derivative integrals, rys polynomials,
c                      derivative three-center integrals
c
c***author             saxe, paul,    (lanl)
c***source             @(#)makd2h.f	5.1   11/28/95
c***purpose            to form the derivatives of the three-center,
c                      two-dimensional integrals h.
c
c***description     
c
c
c***references         (none)
c
c***routines called    (none)
c***end prologue       makd2h.f
      implicit none
c     --- input variables -----
      integer nat,nmax,mmax,imax,jmax,nv,ncentr,nd2,npass
      integer symcen(4)
c     --- input arrays (unmodified) ---
      real*8 dh(nv,3,0:nmax,0:jmax,0:mmax,ncentr)
      real*8 d2g(nv,3,0:nmax,0:mmax,nd2)
      real*8 c(3,nat)
c     --- input arrays (scratch) ---
c     --- output arrays ---
      real*8 d2h(nv,3,0:nmax,0:jmax,0:mmax,nd2)
c     --- output variables ---
c     --- scratch arrays ---
c     --- local variables ---
      integer d2,m,n,i,j,coord,ii
      real*8 scalar
c
c     --- transfer elements from g to h ---
      do 3 d2=1,nd2
         do 2 m=0,mmax
            do 1 n=0,nmax
               call vmove(d2h(1,1,n,0,m,d2),d2g(1,1,n,m,d2),nv*3)
    1       continue
    2    continue
    3    continue
c
c     --- h(i,j,m)=h(i+1,j-1,m)+cimcj*h(i,j-1,m) ---
      do 17 m=0,mmax
         do 16 j=1,jmax
            do 15 i=0,nmax-j
               do 14 coord=1,3
                  scalar=c(coord,symcen(1))-c(coord,symcen(2))
                  do 10 ii=1,nv
                     d2h(ii,coord,i,j,m,1)=d2h(ii,coord,i+1,j-1,m,1)+
     $                    scalar*d2h(ii,coord,i,j-1,m,1)+
     $                    dh(ii,coord,i,j-1,m,1)+
     $                    dh(ii,coord,i,j-1,m,1)
   10             continue
                  if (npass.le.2) then
                     do 11 ii=1,nv
                        d2h(ii,coord,i,j,m,2)=
     $                       d2h(ii,coord,i+1,j-1,m,2)+
     $                       scalar*d2h(ii,coord,i,j-1,m,2)-
     $                       dh(ii,coord,i,j-1,m,1)+
     $                       dh(ii,coord,i,j-1,m,2)
                        d2h(ii,coord,i,j,m,3)=
     $                       d2h(ii,coord,i+1,j-1,m,3)+
     $                       scalar*d2h(ii,coord,i,j-1,m,3)-
     $                       dh(ii,coord,i,j-1,m,2)-
     $                       dh(ii,coord,i,j-1,m,2)
   11                continue
                  else if (npass.eq.3) then
                     do 12 ii=1,nv
                        d2h(ii,coord,i,j,m,2)=
     $                       d2h(ii,coord,i+1,j-1,m,2)+
     $                       scalar*d2h(ii,coord,i,j-1,m,2)+
     $                       dh(ii,coord,i,j-1,m,2)
                        d2h(ii,coord,i,j,m,3)=
     $                       d2h(ii,coord,i+1,j-1,m,3)+
     $                       scalar*d2h(ii,coord,i,j-1,m,3)
   12                continue
                  end if
                  if (npass.eq.1) then
                     do 13 ii=1,nv
                        d2h(ii,coord,i,j,m,4)=
     $                       d2h(ii,coord,i+1,j-1,m,4)+
     $                       scalar*d2h(ii,coord,i,j-1,m,4)+
     $                       dh(ii,coord,i,j-1,m,3)
                        d2h(ii,coord,i,j,m,5)=
     $                       d2h(ii,coord,i+1,j-1,m,5)+
     $                       scalar*d2h(ii,coord,i,j-1,m,5)-
     $                       dh(ii,coord,i,j-1,m,3)
                        d2h(ii,coord,i,j,m,6)=
     $                       d2h(ii,coord,i+1,j-1,m,6)+
     $                       scalar*d2h(ii,coord,i,j-1,m,6)
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
