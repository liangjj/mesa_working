*deck @(#)makdi.f	5.1  11/6/94
      subroutine makdi(i2,di,dh,c,nat,symcen,nmax,mmax,imax,jmax,kmax,
     #                  lmax,nv,ncentr,npass)
c
c***begin prologue     makdi
c***date written       851106   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c***keywords           derivative integrals, rys polynomials,
c                      derivative four-center two-dimensional integrals
c
c***author             saxe, paul,    (lanl)
c***source             @(#)makdi.f	5.1   11/6/94
c***purpose            to calculate the derivatives of the four-center,
c                      two-dimensional integrals i (i2).
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
c***end prologue       makdi
c
      implicit integer (a-z)
c
      real*8 i2(nv,0:imax,0:jmax,0:mmax,0:lmax,3)
      real*8 di(nv,0:imax,0:jmax,0:mmax,0:lmax,3,ncentr)
      real*8 dh(nv,0:nmax,0:jmax,0:mmax,3,ncentr)
      real*8 c(3,nat),scalar
      integer imax,jmax,kmax,lmax,mmax,symcen(4)
c
c     ----- transfer elements from h to i -----
c
      do 5 centr=1,ncentr
      do 4 m=0,mmax
         do 3 j=0,jmax
            do 2 i=0,imax
              call vmove(di(1,i,j,m,0,1,centr),dh(1,i,j,m,1,centr),nv)
              call vmove(di(1,i,j,m,0,2,centr),dh(1,i,j,m,2,centr),nv)
              call vmove(di(1,i,j,m,0,3,centr),dh(1,i,j,m,3,centr),nv)
    2       continue
    3    continue
    4 continue
    5 continue
c
c     ----- i(i,j,k,l)=i(i,j,k+1,l-1)+ckmcl*i(i,j,k,l-1) -----
c
      do 19 l=1,lmax
         do 18 k=0,mmax-l
            do 17 j=0,jmax
               do 16 i=0,imax
                  do 15 coord=1,3
                     scalar=c(coord,symcen(3))-c(coord,symcen(4))
                     do 10 ii=1,nv
                        di(ii,i,j,k,l,coord,1)=
     #                               di(ii,i,j,k+1,l-1,coord,1)+
     #                           scalar*di(ii,i,j,k,l-1,coord,1)
   10               continue
                    if (npass.le.2) then
                       do 11 ii=1,nv
                           di(ii,i,j,k,l,coord,2)=
     #                                  di(ii,i,j,k+1,l-1,coord,2)+
     #                              scalar*di(ii,i,j,k,l-1,coord,2)
   11                  continue
                    else if (npass.eq.3) then
                       do 12 ii=1,nv
                           di(ii,i,j,k,l,coord,2)=
     #                                  di(ii,i,j,k+1,l-1,coord,2)+
     #                              scalar*di(ii,i,j,k,l-1,coord,2)+
     #                                  i2(ii,i,j,k,l-1,coord)
   12                  continue
                    end if
                    if (npass.eq.1) then
                       do 13 ii=1,nv
                           di(ii,i,j,k,l,coord,3)=
     #                                  di(ii,i,j,k+1,l-1,coord,3)+
     #                              scalar*di(ii,i,j,k,l-1,coord,3)+
     #                                  i2(ii,i,j,k,l-1,coord)
   13                  continue
                    end if
   15             continue
   16          continue
   17       continue
   18    continue
   19 continue
c
c      write (6,10) i2
c   10 format (//,' i2:',/,(1x,3f15.9))
c
      return
      end
