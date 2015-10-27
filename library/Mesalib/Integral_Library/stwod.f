*deck @(#)stwod.f	5.1  11/6/94
      subroutine stwod(imax,jmax,nprimi,nprimj,
     $                 exi,exj,alpha,ci,cj,ainv,xyza,expon,xyz,
     $                 nv,nbtype,h,wt,mxpts,a,b)
c***begin prologue     stwod
c***date written       850601   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c***author             saxe, paul (lanl)
c***source
c***keywords           one-electron, integrals, overlap
c***purpose            forms the two-dimensional overlap integrals over
c                      primitives.
c***description
c                      call stwod(iatom,jatom,imax,jmax,nprimi,nprimj,
c                                i1,j1,ex,alpha,c,ainv,xyza,expon,xyz,
c                                nv,nprim,nat,nbtype,h,wt,mxpts,a,b)
c***routines called    lnkerr(mdutil),vadd(math),vinv(math),smul(math),
c                      saxpy(clams),vmul(math),vexp(math),vsqrt(math),
c                      zero(math),vwxs(math),ssub(math),vpower(math),
c                      vfill(math)
c***end prologue       stwod
      implicit integer (a-z)
c
      real*8 exi(nprimi),exj(nprimj),alpha(nv,2),ci(3),cj(3)
      real*8 ainv(nv),xyza(nv,3)
      real*8 expon(nv),xyz(nv,0:imax,0:jmax,3),h(mxpts),wt(mxpts)
      real*8 a(nv),b(nv)
      real*8 ix,jx,scalar,one
c
      parameter (one=1.0d+00)
c
c     ----- start timing -----
c
c
c
      n=0
      do 2 jprim=1,nprimj
         jx=exj(jprim)
         do 1 iprim=1,nprimi
            ix=exi(iprim)
            n=n+1
            alpha(n,1)=ix
            alpha(n,2)=jx
    1    continue
    2 continue
c
      if (n.ne.nv) then
         call lnkerr('l302: stwod error in vector lengths')
      end if
c
      call vadd(ainv,alpha(1,1),alpha(1,2),n)
      call vinv(ainv,ainv,n)
c
c     ----- form xa, ya and za=(ix*xi + jx*xj)/(ix+jx) -----
c
      do 3 coord=1,3
         call smul(xyza(1,coord),alpha(1,1),ci(coord),n)
c        call saxpy(xyza(1,coord),alpha(1,2),cj(coord),n)
         call saxpy(n,cj(coord),alpha(1,2),1,xyza(1,coord),1)
         call vmul(xyza(1,coord),xyza(1,coord),ainv,n)
    3 continue
c
c     ----- form exponential prefactor -----
c
      call vmul(expon,alpha(1,1),alpha(1,2),n)
      scalar=-((ci(1)-cj(1))**2+(ci(2)-cj(2))**2+(ci(3)-cj(3))**2)
      call smul(expon,expon,scalar,n)
      call vmul(expon,expon,ainv,n)
      call vexp(expon,expon,n)
      call vmul(expon,expon,ainv,n)
      call vsqrt(ainv,ainv,n)
      call vmul(expon,expon,ainv,n)
c
c     ----- form two-dimensional integrals -----
c
      call rzero(xyz,n*(imax+1)*(jmax+1)*3)
c
      do 7 ni=0,imax
         do 6 nj=0,jmax
            npts=(ni+nj)/2+1
            minpts=(npts-1)*npts/2+1
            maxpts=minpts+npts-1
            do 5 coord=1,3
               do 4 point=minpts,maxpts
                  call vwxs(b,xyza(1,coord),ainv,h(point),1,n)
                  call ssub(a,b,ci(coord),n)
                  call ssub(b,b,cj(coord),n)
                  if (ni.gt.0) then
                     call vpower(a,a,ni,n)
                  else
                     call vfill(a,one,n)
                  end if
                  if (nj.gt.0) then
                     call vpower(b,b,nj,n)
                  else
                     call vfill(b,one,n)
                  end if
                  call vmul(a,a,b,n)
                  call saxpy(n,wt(point),a,1,xyz(1,ni,nj,coord),1)
c                 call saxpy(xyz(1,ni,nj,coord),a,wt(point),n)
    4          continue
    5       continue
    6    continue
    7 continue
c
c     ----- multiply z 2-d integrals by exponential -----
c
      do 9 j=0,jmax
         do 8 i=0,imax
            call vmul(xyz(1,i,j,3),xyz(1,i,j,3),expon,n)
    8    continue
    9 continue
c
c     ----- stop timing -----
c
c
c
      return
      end
