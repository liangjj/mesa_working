*deck @(#)rotfff.f	5.1  11/6/94
      subroutine rotfff(natoms,tr,fffin,fffout)
c***begin prologue     rotfff.f
c***date written       yymmdd  
c***revision date      11/6/94      
c
c***keywords           
c***author             binkley, et al. (g82)
c***source             @(#)rotfff.f	5.1   11/6/94
c***purpose            
c***description
c     rotate third derivatives from one coordinate system to another.
c     tr is the transformation matrix.  fffin and fffout can be the
c     same.
c
c***references
c
c***routines called
c
c***end prologue       rotfff.f
      implicit none
c     --- input variables -----
      integer natoms
c     --- input arrays (unmodified) ---
c     --- input arrays (scratch) ---
      real*8 fffin(*)
c     --- output arrays ---
      real*8 fffout(*)
c     --- output variables ---
c     --- scratch arrays ---
      real*8 tr(3,3)
c     --- local variables ---
      integer i,j,k,l,ii,jj,kk,i1,i2,i3
      integer iat,jat,kat,locat
      real*8 t(3,3,3),t1(3,3,3)
      real*8 zero
      parameter (zero=0.0d+00)
c
c     --- loop over triplets of atoms.
      do 60 iat=1,natoms
         do 60 jat=1,iat
            do 60 kat=1,jat
c
c              --- pluck out all x,y,z for the current atoms.
               do 10 i=1,3
                  ii=i+3*(iat-1)
                  do 10 j=1,3
                     jj=j+3*(jat-1)
                     do 10 k=1,3
                        kk=k+3*(kat-1)
                        i1=max(ii,jj,kk)
                        i3=min(ii,jj,kk)
                        i2=ii+jj+kk -i1-i3
                        locat=((i1-1)*i1*(i1+1)/6) + (i2*(i2-1)/2) + i3
   10                   t(i,j,k)=fffin(locat)
c
c              --- transform.
               do 20 i=1,3
                  do 20 j=1,3
                     do 20 k=1,3
                        t1(i,j,k)=zero
                        do 20 l=1,3
   20                      t1(i,j,k)=t1(i,j,k) + tr(l,i)*t(l,j,k)
               do 30 i=1,3
                  do 30 j=1,3
                     do 30 k=1,3
                        t(i,j,k)=zero
                        do 30 l=1,3
   30                      t(i,j,k)=t(i,j,k) + tr(l,j)*t1(i,l,k)
               do 40 i=1,3
                  do 40 j=1,3
                     do 40 k=1,3
                        t1(i,j,k)=zero
                        do 40 l=1,3
   40                      t1(i,j,k)=t1(i,j,k) + tr(l,k)*t(i,j,l)
c
c              --- store in output array.
               do 50 i=1,3
                  ii=i+3*(iat-1)
                  do 50 j=1,3
                     jj=j+3*(jat-1)
                     do 50 k=1,3
                        kk=k+3*(kat-1)
                        i1=max(ii,jj,kk)
                        i3=min(ii,jj,kk)
                        i2=ii+jj+kk -i1-i3
                        locat=((i1-1)*i1*(i1+1)/6) + (i2*(i2-1)/2) + i3
   50                   fffout(locat)=t1(i,j,k)
   60 continue
c
c
      return
      end
