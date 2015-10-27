*deck @(#)fermi.f	5.1  11/6/94
      subroutine fermi(iatom,jatom,imax,jmax,nprimi,nprimj,
     $                 i1,j1,ex,alpha,c,expon,xyz,
     $                 nv,nprim,nat,nbtype,ctest)
c***begin prologue     fermi.f
c***date written       860904  
c***revision date      11/6/94      
c
c***keywords           
c***author             martin, richard(lanl)
c***source             @(#)fermi.f	5.1   11/6/94
c***purpose            forms delta functions integrals over primitives 
c***description
c   the point at which the operator is to be evaluated enters in ctest.   
c
c***references
c
c***routines called
c
c***end prologue       fermi.f
      implicit none
c     --- input variables -----
      integer iatom,jatom,imax,jmax,nprimi,nprimj,i1,j1
      integer nv,nprim,nat,nbtype
      real*8 ctest(3)
c     --- input arrays (unmodified) ---
      real*8 ex(nprim),c(3,nat)
c     --- input arrays (scratch) ---
      real*8 expon(nv),alpha(nv,2)
c     --- output arrays ---
      real*8 xyz(nv,0:imax,0:jmax,3)
c     --- output variables ---
c     --- scratch arrays ---
c     --- local variables ---
      integer inp,iout
      integer i,j,i2,j2,n,iprim,jprim,coord,ni,nj
      real*8 ix,jx,one
      real*8 disi,disj,vpoly
c
      parameter (one=1.0d+00)
c
      common/io/inp,iout
c
c     --- load the array containing the exponential argument.
      disi=-((ctest(1)-c(1,iatom))**2+(ctest(2)-c(2,iatom))**2+
     $         (ctest(3)-c(3,iatom))**2)
      disj=-((ctest(1)-c(1,jatom))**2+(ctest(2)-c(2,jatom))**2+
     $         (ctest(3)-c(3,jatom))**2)
      i2=i1+nprimi-1
      j2=j1+nprimj-1
c
      n=0
      do 20 jprim=j1,j2
         jx=ex(jprim)
         do 10 iprim=i1,i2
            ix=ex(iprim)
            n=n+1
            alpha(n,1)=ix*disi
            alpha(n,2)=jx*disj
   10    continue
   20 continue
c
      if(n.ne.nv) then
         call lnkerr('m1902:fermi;  error in vector lengths.')
      endif
c
c     --- evaluate the exponential protion.
      call vadd(expon,alpha(1,1),alpha(1,2),n)
      call vexp(expon,expon,n)
c
c     --- form the polynomial portion of the one-dimensional integrals.
      call rzero(xyz,n*(imax+1)*(jmax+1)*3)
      do 70 ni=0,imax
         do 60 nj=0,jmax
            do 50 coord=1,3
               if(ni.gt.0) then
                  vpoly=(ctest(coord)-c(coord,iatom))**ni
               else
                  vpoly=one
               endif
               if(nj.gt.0) then
                  vpoly=vpoly*((ctest(coord)-c(coord,jatom))**nj)
               else
                  vpoly=vpoly*one
               endif
               call vfill(xyz(1,ni,nj,coord),vpoly,n)
   50       continue
   60    continue
   70 continue
c
c     --- multiply the z 2-d integrals by exponential.
      do 90 j=0,jmax
         do 80 i=0,imax
               call vmul(xyz(1,i,j,3),xyz(1,i,j,3),expon,n)
   80    continue
   90 continue
c
c
      return
      end
