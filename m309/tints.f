*deck @(#)tints.f	5.1  11/6/94
      subroutine tints(iatom,jatom,imax,jmax,nprimi,nprimj,
     #                 i1,j1,ex,alpha,c,ainv,xyza,expon,xyz,
     #                 nv,nprim,nat,nbtype,h,wt,mxpts,a,b,bp2,bm2)
c
c***module to form the two-dimensional integrals over primitives
c   for the kinetic-energy integrals
c
c paul saxe                24 july 1984                  lanl
c
      implicit integer (a-z)
c
      real*8 ex(nprim),alpha(nv,2),c(3,nat),ainv(nv),xyza(nv,3)
      real*8 expon(nv),xyz(nv,0:imax,0:jmax,3,2),h(mxpts),wt(mxpts)
      real*8 a(nv),b(nv),bp2(nv),bm2(nv)
      real*8 ix,jx,scalar,one,half
c
      common /io/     inp,iout
c
      parameter (one=1.0d+00,half=0.5d+00)
c
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
         write (iout,99) n,nv
   99    format (//,' ##### m302: tints, error, n and nv:',2i8,//)
         call lnkerr('m302: tints error in vector lengths')
      end if
c
      call vadd(ainv,alpha(1,1),alpha(1,2),n)
      call vinv(ainv,ainv,n)
c
c     ----- form xa, ya and za=(ix*xi + jx*xj)/(ix+jx) -----
c
      do 3 coord=1,3
         call smul(xyza(1,coord),alpha(1,1),c(coord,iatom),n)
         call saxpy(n,c(coord,jatom),alpha(1,2),1,xyza(1,coord),1)
c        call saxpy(xyza(1,coord),alpha(1,2),c(coord,jatom),n)
         call vmul(xyza(1,coord),xyza(1,coord),ainv,n)
    3 continue
c
c     ----- form exponential prefactor -----
c
      call vmul(expon,alpha(1,1),alpha(1,2),n)
      scalar=-((c(1,iatom)-c(1,jatom))**2+(c(2,iatom)-c(2,jatom))**2+
     #         (c(3,iatom)-c(3,jatom))**2)
      call smul(expon,expon,scalar,n)
      call vmul(expon,expon,ainv,n)
      call vexp(expon,expon,n)
      call vmul(expon,expon,ainv,n)
      call vsqrt(ainv,ainv,n)
      call vmul(expon,expon,ainv,n)
c
c     ----- form two-dimensional integrals -----
c
      call rzero(xyz,n*(imax+1)*(jmax+1)*3*2)
c
c     ----- double the values of alpha(j) ------
c
      call vadd(alpha(1,2),alpha(1,2),alpha(1,2),n)
c
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
c     ----- form del**2 on j primitive -----
c
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
                  scalar=nj*(nj-1)
                  call smul(bm2,bm2,scalar,n)
                  call vmul(bp2,bp2,alpha(1,2),n)
                  scalar=-(2*nj+1)
                  call saxpy(n,scalar,b,1,bp2,1)
c                 call saxpy(bp2,b,scalar,n)
                  call vmul(bp2,bp2,alpha(1,2),n)
                  call vadd(bm2,bm2,bp2,n)
                  call smul(bm2,bm2,-half,n)
                  call vmul(bm2,bm2,a,n)
                  call saxpy(n,wt(point),bm2,1,xyz(1,ni,nj,coord,2),1)
c                 call saxpy(xyz(1,ni,nj,coord,2),bm2,wt(point),n)
c
                  call vmul(a,a,b,n)
                  call saxpy(n,wt(point),a,1,xyz(1,ni,nj,coord,1),1)
c                 call saxpy(xyz(1,ni,nj,coord,1),a,wt(point),n)
    4          continue
    5       continue
    6    continue
    7 continue
c
c     ----- multiply z 2-d integrals by exponential -----
c
      do 9 j=0,jmax
         do 8 i=0,imax
            call vmul(xyz(1,i,j,3,1),xyz(1,i,j,3,1),expon,n)
            call vmul(xyz(1,i,j,3,2),xyz(1,i,j,3,2),expon,n)
    8    continue
    9 continue
c
c
c
      return
      end
