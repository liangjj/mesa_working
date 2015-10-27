*deck @(#)vints.f	5.1  11/6/94
      subroutine vints(iatom,jatom,imax,jmax,nprimi,nprimj,
     #                 i1,j1,ex,alpha,c,ainv,xyza,expon,xyz,
     #                 nv,nprim,nat,nbtype,h,wt,mxpts,a,b,
     #                 aiaj,xyz0,rysrt,ryswt,nroots,prmint,
     #                 lenblk,zan,xyz1,t1,mini,maxi,minj,maxj,
     #                 nx,ny,nz)
c
c***module to form the two-dimensional integrals over primitives
c   for the potential-energy integrals
c
c paul saxe                26 july 1984                  lanl
c
      implicit integer (a-z)
c
      real*8 ex(nprim),alpha(nv,2),c(3,nat),ainv(nv),xyza(nv,3)
      real*8 expon(nv),xyz(nv,0:imax,0:jmax,3),h(mxpts),wt(mxpts)
      real*8 a(nv),b(nv),aiaj(nv),xyz0(nv,3),rysrt(nv,nroots)
      real*8 ryswt(nv,nroots),prmint(nv,lenblk),zan(nat)
      real*8 xyz1(nv,3),t1(nv)
      real*8 iex,jex,scalar,one,temprt(9),tempwt(9),pi212,xx
      integer nx(*),ny(*),nz(*)
c
      common/io/inp,iout
c
      parameter (pi212=1.1283791670955d+00, one=1.0d+00)
c
      i2=i1+nprimi-1
      j2=j1+nprimj-1
c
      n=0
      do 2 jprim=j1,j2
         jex=ex(jprim)
         do 1 iprim=i1,i2
            iex=ex(iprim)
            n=n+1
            alpha(n,1)=iex
            alpha(n,2)=jex
    1    continue
    2 continue
c
      if (n.ne.nv) then
         write (iout,99) n,nv
   99    format (//,' ##### m302: vints, error, n and nv:',2i8,//)
         call lnkerr('m302: vints error in vector lengths')
      end if
c
      call vadd(aiaj,alpha(1,1),alpha(1,2),n)
c
c     ----- form xa, ya and za=(ix*xi + jx*xj) -----
c
      do 3 coord=1,3
         call smul(xyza(1,coord),alpha(1,1),c(coord,iatom),n)
         call saxpy(n,c(coord,jatom),alpha(1,2),1,xyza(1,coord),1)
c        call saxpy(xyza(1,coord),alpha(1,2),c(coord,jatom),n)
         call vdiv(xyz1(1,coord),xyza(1,coord),aiaj,n)
    3 continue
c
c     ----- form exponential prefactor -----
c
      call vmul(expon,alpha(1,1),alpha(1,2),n)
      scalar=-((c(1,iatom)-c(1,jatom))**2+(c(2,iatom)-c(2,jatom))**2+
     #         (c(3,iatom)-c(3,jatom))**2)
      call smul(expon,expon,scalar,n)
      call vdiv(expon,expon,aiaj,n)
      call vexp(expon,expon,n)
      call vdiv(expon,expon,aiaj,n)
c
c     ----- zero space for the accumulation of primitive integrals -----
c
      call rzero(prmint,n*lenblk)
c
c     ----- form two-dimensional integrals -----
c
      do 16 catom=1,nat
c
c     ----- calculate roots and weights of rys polynomials -----
c
         do 5 i=1,n
            xx=aiaj(i)*((xyz1(i,1)-c(1,catom))**2+
     #                  (xyz1(i,2)-c(2,catom))**2+
     #                  (xyz1(i,3)-c(3,catom))**2)
            call roots(nroots,xx,temprt,tempwt,1,junk,junk,junk,
     #                 junk,junk,junk)
            do 4 root=1,nroots
               rysrt(i,root)=temprt(root)
               ryswt(i,root)=tempwt(root)
    4       continue
    5    continue
c
         do 15 root=1,nroots
            call vmul(rysrt(1,root),rysrt(1,root),aiaj,n)
            scalar=-pi212*zan(catom)
            call smul(ryswt(1,root),ryswt(1,root),scalar,n)
            call vmul(ryswt(1,root),ryswt(1,root),expon,n)
            call vadd(ainv,aiaj,rysrt(1,root),n)
            call vinv(ainv,ainv,n)
            do 6 coord=1,3
               call vwxs(xyz0(1,coord),xyza(1,coord),rysrt(1,root),
     #                                      c(coord,catom),1,n)
               call vmul(xyz0(1,coord),xyz0(1,coord),ainv,n)
    6       continue
c
c     ----- form the two-dimensional integrals -----
c
            call vsqrt(ainv,ainv,n)
            call rzero(xyz,n*(imax+1)*(jmax+1)*3)
c
            do 10 ni=0,imax
               do 9 nj=0,jmax
                  npts=(ni+nj)/2+1
                  minpts=(npts-1)*npts/2+1
                  maxpts=minpts+npts-1
                  do 8 coord=1,3
                     do 7 point=minpts,maxpts
                        call vwxs(b,xyz0(1,coord),ainv,h(point),1,n)
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
c                       call saxpy(xyz(1,ni,nj,coord),a,wt(point),n)
    7                continue
    8             continue
    9          continue
   10       continue
c
c     ----- multiply z 2-d integrals by exponential -----
c
            do 12 j=0,jmax
               do 11 i=0,imax
                  call vmul(xyz(1,i,j,3),xyz(1,i,j,3),ryswt(1,root),n)
   11          continue
   12       continue
c
c     ----- form the primitive integrals -----
c
            intgrl=0
            do 14 i=mini,maxi
               ix=nx(i)
               iy=ny(i)
               iz=nz(i)
               do 13 j=minj,maxj
                  jx=nx(j)
                  jy=ny(j)
                  jz=nz(j)
                  intgrl=intgrl+1
c
                  call vmul(t1,xyz(1,ix,jx,1),xyz(1,iy,jy,2),n)
                  call vmul(t1,t1,xyz(1,iz,jz,3),n)
                  call vadd(prmint(1,intgrl),prmint(1,intgrl),t1,n)
   13          continue
   14       continue
   15    continue
   16 continue
c
c
c
      return
      end
