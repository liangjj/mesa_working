*deck @(#)sints.f	5.1  11/6/94
      subroutine sints(iatom,jatom,imax,jmax,nprimi,nprimj,
     #                 i1,j1,ex,alpha,c,ainv,xyza,expon,xyz,
     #                 nv,nprim,nat,nbtype,h,wt,mxpts,a,b,
     #                 nderiv)
c
c***module to form the two-dimensional integrals over primitives
c   for the overlap integrals
c
c paul saxe                19 july 1984                  lanl
c
      implicit integer (a-z)
c
      real*8 ex(nprim),alpha(nv,2),c(3,nat),ainv(nv),xyza(nv,3)
      real*8 expon(nv),xyz(nv,0:imax,0:jmax,3,0:nderiv),h(mxpts)
      real*8 a(nv,-nderiv:nderiv),b(nv),wt(mxpts)
      real*8 ix,jx,scalar,one
c
      parameter (one=1.0d+00)
c
c     ----- start timing -----
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
         call lnkerr('m303: sints error in vector lengths')
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
      call rzero(xyz,n*(imax+1)*(jmax+1)*3*(nderiv+1))
c
      do 7 ni=0,imax
         do 6 nj=0,jmax
            npts=(ni+nj+nderiv)/2+1
            minpts=(npts-1)*npts/2+1
            maxpts=minpts+npts-1
            do 5 coord=1,3
               do 4 point=minpts,maxpts
                  call vwxs(b,xyza(1,coord),ainv,h(point),1,n)
                  call ssub(a(1,0),b,c(coord,iatom),n)
                  call ssub(b,b,c(coord,jatom),n)
                  if (nderiv.eq.1) then
                     if (ni.gt.1) then
                        call vpower(a(1,-1),a(1,0),ni-1,n)
                        call vpower(a(1, 1),a(1,0),ni+1,n)
                        call vpower(a(1, 0),a(1,0),ni  ,n)
                     else if (ni.eq.1) then
                        call vpower(a(1, 1),a(1,0),2   ,n)
                        call vfill(a(1,-1),one,n)
                     else
                        call vmove(a(1,1),a(1,0),n)
                        call vfill(a(1,0),one,n)
                        call rzero(a(1,-1),n)
                     end if
                  else if (nderiv.eq.0) then
                     if (ni.gt.0) then
                        call vpower(a(1,0),a(1,0),ni,n)
                     else
                        call vfill(a(1,0),one,n)
                     end if
                  end if
                  if (nj.gt.0) then
                     call vpower(b,b,nj,n)
                  else
                     call vfill(b,one,n)
                  end if
c                 call vmul(a,a,b,n)
c                 call saxpy(n,wt(point),a,1,xyz(1,ni,nj,coord,1),1)
                  do 20 i=1,n
                     xyz(i,ni,nj,coord,0)=xyz(i,ni,nj,coord,0)+
     #                      a(i,0)*b(i)*wt(point)
                     xyz(i,ni,nj,coord,1)=xyz(i,ni,nj,coord,1)+
     #                      (2*alpha(i,1)*a(i,+1)-ni*a(i,-1))*
     #                            b(i)*wt(point)
   20             continue
    4          continue
    5       continue
    6    continue
    7 continue
c
c     ----- multiply z 2-d integrals by exponential -----
c
      do 9 j=0,jmax
         do 8 i=0,imax
            call vmul(xyz(1,i,j,3,0),xyz(1,i,j,3,0),expon,n)
            if (nderiv.ge.1) then
               call vmul(xyz(1,i,j,3,1),xyz(1,i,j,3,1),expon,n)
            end if
    8    continue
    9 continue
c
c     ----- stop timing -----
c
c
c
      return
      end
