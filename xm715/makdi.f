*deck %W%  %G%
      subroutine makdi(i2,di,dh,c,nat,symcen,nmax,mmax,imax,jmax,kmax,
     $                  lmax,nv,ncentr,npass)
c***begin prologue     %W%   %G%
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
c***end prologue       makdi.f
      implicit none
c     --- input variables -----
      integer nat,nmax,mmax,imax,jmax,kmax,lmax
      integer nv,ncentr,npass
c     --- input arrays (unmodified) ---
      integer symcen(4)
      real*8 dh(nv,3,0:nmax,0:jmax,0:mmax,ncentr)
      real*8 i2(nv,3,0:imax,0:jmax,0:mmax,0:lmax)
      real*8 c(3,nat)
c     --- input arrays (scratch) ---
c     --- output arrays ---
      real*8 di(nv,3,0:imax,0:jmax,0:mmax,0:lmax,ncentr)
c     --- output variables ---
c     --- scratch arrays ---
c     --- local variables ---
      integer centr,m,i,j,k,l,coord,ii
      integer inp,iout
      logical debug
      real*8 scalar
c
      parameter (debug=.false.)
c
      common/io/inp,iout
c
 1000 format (//,' i2:',/,(1x,3f15.9))
c
c     --- transfer elements from h to i ---
      do 5 centr=1,ncentr
         do 4 m=0,mmax
            do 3 j=0,jmax
               do 2 i=0,imax
                  call vmove(di(1,1,i,j,m,0,centr),
     $                       dh(1,1,i,j,m,centr),nv*3)
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
                        di(ii,coord,i,j,k,l,1)=
     $                               di(ii,coord,i,j,k+1,l-1,1)+
     $                           scalar*di(ii,coord,i,j,k,l-1,1)
   10               continue
                    if (npass.le.2) then
                       do 11 ii=1,nv
                           di(ii,coord,i,j,k,l,2)=
     $                                  di(ii,coord,i,j,k+1,l-1,2)+
     $                              scalar*di(ii,coord,i,j,k,l-1,2)
   11                  continue
                    else if (npass.eq.3) then
                       do 12 ii=1,nv
                           di(ii,coord,i,j,k,l,2)=
     $                                  di(ii,coord,i,j,k+1,l-1,2)+
     $                              scalar*di(ii,coord,i,j,k,l-1,2)+
     $                                  i2(ii,coord,i,j,k,l-1)
   12                  continue
                    end if
                    if (npass.eq.1) then
                       do 13 ii=1,nv
                           di(ii,coord,i,j,k,l,3)=
     $                                  di(ii,coord,i,j,k+1,l-1,3)+
     $                              scalar*di(ii,coord,i,j,k,l-1,3)+
     $                                  i2(ii,coord,i,j,k,l-1)
   13                  continue
                    end if
   15             continue
   16          continue
   17       continue
   18    continue
   19 continue
c
      if(debug) then
         write (iout,1000) i2
      endif
c
c
      return
      end