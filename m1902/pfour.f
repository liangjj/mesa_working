*deck @(#)pfour.f	5.1  11/6/94
      subroutine pfour(iatom,jatom,imax,jmax,nprimi,nprimj,
     $                 i1,j1,ex,alpha,c,ainv,xyza,expon,xyz,
     $                 nv,nprim,nat,nbtype,h,wt,mxpts,a,b,
     $                 ap2,am2,bp2,bm2)
c***begin prologue     pfour.f
c***date written       840724  
c***revision date      11/6/94      
c
c***keywords           
c***author             martin,richard, and saxe,paul(lanl)
c***source             @(#)pfour.f	5.1   11/6/94
c***purpose            forms primitive integrals over the del**4 operator
c***description
c     
c    
c
c***references
c
c***routines called
c
c***end prologue       pfour.f
      implicit none
c     --- input variables -----
      integer iatom,jatom,imax,jmax,nprimi,nprimj,i1,j1
      integer nv,nprim,nat,nbtype,mxpts
c     --- input arrays (unmodified) ---
      real*8 ex(nprim),c(3,nat)
      real*8 h(mxpts),wt(mxpts)
c     --- input arrays (scratch) ---
      real*8 expon(nv),alpha(nv,2)
      real*8 a(nv),b(nv),ap2(nv),am2(nv),bp2(nv),bm2(nv)
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
      common /io/     inp,iout
c
 1000 format('m1902: pfour, error, n and nv:',2i8)
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
         call lnkerr('m1902: pfour error in vector lengths')
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
c     --- double the values of alpha(i) and alpha(j).
      call vadd(alpha(1,1),alpha(1,1),alpha(1,1),n)
      call vadd(alpha(1,2),alpha(1,2),alpha(1,2),n)
c
c     --- evaluate the one-dimensional integrals
      do 7 ni=0,imax
         do 6 nj=0,jmax
            npts=(ni+nj)/2+2
            minpts=(npts-1)*npts/2+1
            maxpts=minpts+npts-1
            if(maxpts.gt.mxpts)
     $         call lnkerr('m1902:pfour; not enough rys points.')
            do 5 coord=1,3
               do 4 point=minpts,maxpts
                  call vwxs(b,xyza(1,coord),ainv,h(point),1,n)
                  call ssub(a,b,c(coord,iatom),n)
                  call ssub(b,b,c(coord,jatom),n)
c
c                 --- form del**2 on i primitive.
                  call vpower(ap2,a,ni+2,n)
                  if (ni-2.gt.0) then
                     call vpower(am2,a,ni-2,n)
                  else if (ni-2.eq.0) then
                     call vfill(am2,one,n)
                  else
                     call rzero(am2,n)
                  end if
                  if (ni.gt.0) then
                     call vpower(a,a,ni,n)
                  else
                     call vfill(a,one,n)
                  end if
c
c                 --- form del**2 on j primitive.
                  call vpower(bp2,b,nj+2,n)
                  if (nj-2.gt.0) then
                     call vpower(bm2,b,nj-2,n)
                  else if (nj-2.eq.0) then
                     call vfill(bm2,one,n)
                  else
                     call rzero(bm2,n)
                  end if
                  if (nj.gt.0) then
                     call vpower(b,b,nj,n)
                  else
                     call vfill(b,one,n)
                  end if
c
                  scalar=ni*(ni-1)
                  call smul(am2,am2,scalar,n)
                  call vmul(ap2,ap2,alpha(1,1),n)
                  scalar=-(2*ni+1)
                  call saxpy(n,scalar,a,1,ap2,1)
                  call vmul(ap2,ap2,alpha(1,1),n)
                  call vadd(am2,am2,ap2,n)
                  call smul(am2,am2,-half,n)
c
                  scalar=nj*(nj-1)
                  call smul(bm2,bm2,scalar,n)
                  call vmul(bp2,bp2,alpha(1,2),n)
                  scalar=-(2*nj+1)
                  call saxpy(n,scalar,b,1,bp2,1)
                  call vmul(bp2,bp2,alpha(1,2),n)
                  call vadd(bm2,bm2,bp2,n)
                  call smul(bm2,bm2,-half,n)
c
                  call vmul(bm2,bm2,am2,n)
                  call saxpy(n,wt(point),bm2,1,xyz(1,ni,nj,coord,2),1)
c
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
