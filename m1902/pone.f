*deck @(#)pone.f	5.1  11/6/94
      subroutine pone(iatom,jatom,imax,jmax,nprimi,nprimj,
     $                 i1,j1,ex,alpha,c,ainv,xyza,expon,xyz,
     $                 nv,nprim,nat,nbtype,h,wt,mxpts,a,b,bp1,bm1)
c
c***begin prologue     pone.f
c***date written       840724  
c***revision date      11/6/94      
c
c***keywords           
c***author             martin,richard, and saxe,paul(lanl)
c***source             @(#)pone.f	5.1   11/6/94
c***purpose            forms integrals over primitives for momentum
c                      integrals
c***description
c     
c    
c
c***references
c
c***routines called
c
c***end prologue       pone.f
      implicit none
c     --- input variables -----
      integer iatom,jatom,imax,jmax,nprimi,nprimj,i1,j1
      integer nv,nprim,nat,nbtype,mxpts
c     --- input arrays (unmodified) ---
      real*8 ex(nprim),c(3,nat)
      real*8 h(mxpts),wt(mxpts)
c     --- input arrays (scratch) ---
      real*8 expon(nv),alpha(nv,2)
      real*8 a(nv),b(nv),bp1(nv),bm1(nv)
      real*8 ainv(nv),xyza(nv,3)
c     --- output arrays ---
      real*8 xyz(nv,0:imax,0:jmax,3,2)
c     --- output variables ---
c     --- scratch arrays ---
c     --- local variables ---
      integer inp,iout
      integer i,j,i2,j2,n,iprim,jprim,coord,ni,nj
      integer npts,maxpts,minpts,point
      real*8 ix,jx,scalar,one,half
c
      parameter (one=1.0d+00,half=0.5d+00)
c
      common/io/inp,iout
c
 1000 format(1x,'m1902:pone, error, n and nv:',2i8)
c
c     --- assemble the exponents
      i2=i1+nprimi-1
      j2=j1+nprimj-1
c
      n=0
      do 2 jprim=j1,j2
         jx=ex(jprim)
         do 1 iprim=i1,i2
            ix=ex(iprim)
            n=n+1
            alpha(n,1)=ix
            alpha(n,2)=jx
    1    continue
    2 continue
c
      if (n.ne.nv) then
         write (iout,1000) n,nv
         call lnkerr('m1902:pone error in vector lengths')
      end if
c
      call vadd(ainv,alpha(1,1),alpha(1,2),n)
      call vinv(ainv,ainv,n)
c
c     --- form xa, ya and za=(ix*xi + jx*xj)/(ix+jx) ---
      do 3 coord=1,3
         call smul(xyza(1,coord),alpha(1,1),c(coord,iatom),n)
         call saxpy(n,c(coord,jatom),alpha(1,2),1,xyza(1,coord),1)
         call vmul(xyza(1,coord),xyza(1,coord),ainv,n)
    3 continue
c
c     --- form exponential prefactor ---
      call vmul(expon,alpha(1,1),alpha(1,2),n)
      scalar=-((c(1,iatom)-c(1,jatom))**2+(c(2,iatom)-c(2,jatom))**2+
     $         (c(3,iatom)-c(3,jatom))**2)
      call smul(expon,expon,scalar,n)
      call vmul(expon,expon,ainv,n)
      call vexp(expon,expon,n)
      call vmul(expon,expon,ainv,n)
      call vsqrt(ainv,ainv,n)
      call vmul(expon,expon,ainv,n)
c
c     --- form two-dimensional integrals ---
      call rzero(xyz,n*(imax+1)*(jmax+1)*3*2)
c
c     --- double the values of alpha(j) ---
      call vadd(alpha(1,2),alpha(1,2),alpha(1,2),n)
c
c     --- evaluate the one-dimensional integrals
      do 7 ni=0,imax
         do 6 nj=0,jmax
            npts=(ni+nj)/2+2
            minpts=(npts-1)*npts/2+1
            maxpts=minpts+npts-1
            do 5 coord=1,3
               do 4 point=minpts,maxpts
                  call vwxs(b,xyza(1,coord),ainv,h(point),1,n)
                  call ssub(a,b,c(coord,iatom),n)
                  call ssub(b,b,c(coord,jatom),n)
                  if (ni.gt.0) then
                     call vpower(a,a,ni,n)
                  else
                     call vfill(a,one,n)
                  end if
c
c                 --- form del on j primitive ---
                  call vpower(bp1,b,nj+1,n)
                  if (nj-1.gt.0) then
                     call vpower(bm1,b,nj-1,n)
                  else if (nj-1.eq.0) then
                     call vfill(bm1,one,n)
                  else
                     call rzero(bm1,n)
                  end if
                  if (nj.gt.0) then
                     call vpower(b,b,nj,n)
                  else
                     call vfill(b,one,n)
                  end if
                  scalar=nj
                  call smul(bm1,bm1,scalar,n)
                  call vmul(bp1,bp1,alpha(1,2),n)
                  call vsub(bm1,bm1,bp1,n)
                  call vmul(bm1,bm1,a,n)
                  call saxpy(n,wt(point),bm1,1,xyz(1,ni,nj,coord,2),1)
                  call vmul(a,a,b,n)
                  call saxpy(n,wt(point),a,1,xyz(1,ni,nj,coord,1),1)
    4          continue
    5       continue
    6    continue
    7 continue
c
c     --- multiply z 2-d integrals by exponential ---
      do 9 j=0,jmax
         do 8 i=0,imax
            call vmul(xyz(1,i,j,3,1),xyz(1,i,j,3,1),expon,n)
            call vmul(xyz(1,i,j,3,2),xyz(1,i,j,3,2),expon,n)
    8    continue
    9 continue
c
c
      return
      end
