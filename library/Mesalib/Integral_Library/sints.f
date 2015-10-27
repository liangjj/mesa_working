*deck  @(#)sints.f	5.1 11/6/94
      subroutine sints(iatom,jatom,imax,jmax,nprimi,nprimj,
     #                 i1,j1,ex,alpha,c,ainv,xyza,expon,xyz,
     #                 nv,nprim,nat,nbtype,h,wt,mxpts,a,b)
c***begin prologue     sints.f
c***date written       840719  
c***revision date      11/6/94      
c
c***keywords           
c***author             saxe, paul(lanl) 
c***source             @(#)sints.f	5.1   11/6/94
c***purpose            
c      module to form the two-dimensional integrals over primitives
c      for the overlap integrals
c***description
c     
c    
c
c***references
c
c***routines called
c
c***end prologue       sints.f
      implicit integer (a-z)
c
      real*8 ex(nprim),alpha(nv,2),c(3,nat),ainv(nv),xyza(nv,3)
      real*8 expon(nv),xyz(nv,0:imax,0:jmax,3),h(mxpts),wt(mxpts)
      real*8 a(nv),b(nv)
      real*8 ix,jx,scalar,one
c
      parameter (one=1.0d+00)
c
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
         call lnkerr('sints: error in vector lengths')
      end if
c
      call vadd(ainv,alpha(1,1),alpha(1,2),n)
      call vinv(ainv,ainv,n)
c
c     ----- form xa, ya and za=(ix*xi + jx*xj)/(ix+jx) -----
c
      do 3 coord=1,3
         call smul(xyza(1,coord),alpha(1,1),c(coord,iatom),n)
c        call saxpy(xyza(1,coord),alpha(1,2),c(coord,jatom),n)
         call saxpy(n,c(coord,jatom),alpha(1,2),1,xyza(1,coord),1)
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
      call rzero(xyz,n*(imax+1)*(jmax+1)*3)
c
      do 7 ni=0,imax
         do 6 nj=0,jmax
            npts=(ni+nj)/2+1
            minpts=(npts-1)*npts/2+1
            maxpts=minpts+npts-1
            if(maxpts.gt.mxpts)
     $         call lnkerr('sints; not enough rys points.')
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
c
      return
      end
