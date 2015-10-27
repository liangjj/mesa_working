*deck %W%  %G%
      subroutine makd2i(di,d2i,d2h,c,nat,symcen,nmax,mmax,imax,jmax,
     $     kmax,lmax,nv,ncentr,nd2,npass)
c***begin prologue     makd2i.f
c***date written       851106   (yymmdd)
c***revision date      11/6/94
c***keywords           derivative integrals, rys polynomials,
c                      derivative four-center two-dimensional integrals
c
c***author             saxe, paul,    (lanl)
c***source             %W%   %G%
c***purpose            to calculate the derivatives of the four-center,
c                      two-dimensional integrals i (i2).
c
c***description         
c
c
c***references         (none)
c
c***routines called    (none)
c***end prologue       makd2i.f
      implicit none
c     --- input variables -----
      integer nat,nmax,mmax,imax,jmax,kmax,lmax
      integer nv,ncentr,nd2,npass
c     --- input arrays (unmodified) ---
      integer symcen(4)
      real*8 di(nv,3,0:imax,0:jmax,0:mmax,0:lmax,ncentr)
      real*8 d2h(nv,3,0:nmax,0:jmax,0:mmax,nd2)
      real*8 c(3,nat)
c     --- input arrays (scratch) ---
c     --- output arrays ---
      real*8 d2i(nv,3,0:imax,0:jmax,0:mmax,0:lmax,nd2)
c     --- output variables ---
c     --- scratch arrays ---
c     --- local variables ---
      integer d2,m,i,j,k,l,coord,ii
      integer inp,iout
      logical debug
      real*8 scalar
c
      parameter (debug=.false.)
      common/io/inp,iout
c
 1000 format (/,' d2i:',/,(1x,3f15.9))
c
c     --- transfer elements from h to i ---
c
      do 5 d2=1,nd2
         do 4 m=0,mmax
            do 3 j=0,jmax
               do 2 i=0,imax
                  call vmove(d2i(1,1,i,j,m,0,d2),d2h(1,1,i,j,m,d2),
     $                                                       nv*3)
    2          continue
    3       continue
    4    continue
    5 continue
c
c     --- i(i,j,k,l)=i(i,j,k+1,l-1)+ckmcl*i(i,j,k,l-1) ---
      do 19 l=1,lmax
         do 18 k=0,mmax-l
            do 17 j=0,jmax
               do 16 i=0,imax
                  do 15 coord=1,3
                     scalar=c(coord,symcen(3))-c(coord,symcen(4))
                     do 10 ii=1,nv
                        d2i(ii,coord,i,j,k,l,1)=
     $                       d2i(ii,coord,i,j,k+1,l-1,1)+
     $                       scalar*d2i(ii,coord,i,j,k,l-1,1)
   10               continue
                    if (npass.le.2) then
                       do 11 ii=1,nv
                           d2i(ii,coord,i,j,k,l,2)=
     $                         d2i(ii,coord,i,j,k+1,l-1,2)+
     $                         scalar*d2i(ii,coord,i,j,k,l-1,2)
                           d2i(ii,coord,i,j,k,l,3)=
     $                         d2i(ii,coord,i,j,k+1,l-1,3)+
     $                         scalar*d2i(ii,coord,i,j,k,l-1,3)
   11                  continue
                    else if (npass.eq.3) then
                       do 12 ii=1,nv
                           d2i(ii,coord,i,j,k,l,2)=
     $                         d2i(ii,coord,i,j,k+1,l-1,2)+
     $                         scalar*d2i(ii,coord,i,j,k,l-1,2)+
     $                         di(ii,coord,i,j,k,l-1,1)
                           d2i(ii,coord,i,j,k,l,3)=
     $                         d2i(ii,coord,i,j,k+1,l-1,3)+
     $                         scalar*d2i(ii,coord,i,j,k,l-1,3)+
     $                         di(ii,coord,i,j,k,l-1,2)+
     $                         di(ii,coord,i,j,k,l-1,2)
   12                  continue
                    end if
                    if (npass.eq.1) then
                       do 13 ii=1,nv
                           d2i(ii,coord,i,j,k,l,4)=
     $                         d2i(ii,coord,i,j,k+1,l-1,4)+
     $                         scalar*d2i(ii,coord,i,j,k,l-1,4)+
     $                         di(ii,coord,i,j,k,l-1,1)
                           d2i(ii,coord,i,j,k,l,5)=
     $                         d2i(ii,coord,i,j,k+1,l-1,5)+
     $                         scalar*d2i(ii,coord,i,j,k,l-1,5)+
     $                         di(ii,coord,i,j,k,l-1,2)
                           d2i(ii,coord,i,j,k,l,6)=
     $                         d2i(ii,coord,i,j,k+1,l-1,6)+
     $                         scalar*d2i(ii,coord,i,j,k,l-1,6)+
     $                         di(ii,coord,i,j,k,l-1,3)+
     $                         di(ii,coord,i,j,k,l-1,3)
   13                  continue
                    end if
   15             continue
   16          continue
   17       continue
   18    continue
   19 continue
c
      if(debug) then
         write (iout,1000) d2i
      endif
c
c
      return
      end
