*deck @(#)vints.f	5.1  11/6/94
      subroutine vints(iatom,jatom,imax,jmax,nprimi,nprimj,
     #                 i1,j1,ex,alpha,c,ainv,xyza,expon,xyz,
     #                 nv,nprim,nat,nbtype,h,wt,mxpts,a,b,
     #                 aiaj,xyz0,rysrt,ryswt,nroots,prmint,
     #                 lenblk,zan,xyz1,t1,mini,maxi,minj,maxj,
     #                 nx,ny,nz,nderiv,conint,nconti,ncontj,
     #                 itype,jtype,acf,bcf,tmp1,len1,imin,jmin,nocart,
     #                 ds,nnp,start,nobf)
c
c***module to form the two-dimensional integrals over primitives
c   for the potential-energy integrals
c
c paul saxe                26 july 1984                  lanl
c
      implicit integer (a-z)
c
      real*8 ex(nprim),alpha(nv,2),c(3,nat),ainv(nv),xyza(nv,3)
      real*8 expon(nv),xyz(nv,0:imax,0:jmax,3,0:2)
      real*8 a(nv,-nderiv:nderiv),b(nv,-nderiv:nderiv),aiaj(nv)
      real*8 ryswt(nv,nroots),prmint(nv,lenblk,7),zan(nat)
      real*8 xyz1(nv,3),t1(nv),h(mxpts),wt(mxpts),rysrt(nv,nroots)
      real*8 xyz0(nv,3)
      real*8 conint(nconti,ncontj,lenblk),tmp1(len1),ds(nnp,3,nat)
      real*8 acf(nprimi,nprimj,imin:imax),bcf(nprimj,ncontj,jmin:jmax)
      real*8 iex,jex,scalar,one,temprt(9),tempwt(9),pi212,xx
      integer nx(*),ny(*),nz(*)
      integer nocart(nbtype),start(nat,nbtype),nobf(nbtype)
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
   99    format (//,' ##### m303: vints, error, n and nv:',2i8,//)
         call lnkerr('m303: vints error in vector lengths')
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
c     ----- zero space for the accumulation of primitive integrals ---
c
      call rzero(prmint,n*lenblk*7)
c
c     ----- form two-dimensional integrals -----
c
      do 16 catom=1,nat
c
c     ----- calculate roots and weights of rys polynomials -----
c
         call rzero(prmint(1,1,2),n*lenblk*6)
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
            call rzero(xyz,n*(imax+1)*(jmax+1)*3*3)
c
            do 10 ni=0,imax
               do 9 nj=0,jmax
                  npts=(ni+nj+nderiv)/2+1
                  minpts=(npts-1)*npts/2+1
                  maxpts=minpts+npts-1
                  do 8 coord=1,3
                     do 7 point=minpts,maxpts
                        call vwxs(b(1,0),xyz0(1,coord),ainv,h(point),
     #                                     1,n)
                        call ssub(a(1,0),b(1,0),c(coord,iatom),n)
                        call ssub(b(1,0),b(1,0),c(coord,jatom),n)
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
                        if (nderiv.eq.1) then
                           if (nj.gt.1) then
                              call vpower(b(1,-1),b(1,0),nj-1,n)
                              call vpower(b(1, 1),b(1,0),nj+1,n)
                              call vpower(b(1, 0),b(1,0),nj  ,n)
                           else if (nj.eq.1) then
                              call vpower(b(1, 1),b(1,0),2   ,n)
                              call vfill(b(1,-1),one,n)
                           else
                              call vmove(b(1,1),b(1,0),n)
                              call vfill(b(1,0),one,n)
                              call rzero(b(1,-1),n)
                           end if
                        else if (nderiv.eq.0) then
                           if (nj.gt.0) then
                              call vpower(b(1,0),b(1,0),nj,n)
                           else
                              call vfill(b(1,0),one,n)
                           end if
                        end if
                        do 50 i=1,n
                           xyz(i,ni,nj,coord,0)=xyz(i,ni,nj,coord,0)+
     #                               a(i,0)*b(i,0)*wt(point)
                           if (iatom.ne.catom) then
                             xyz(i,ni,nj,coord,1)=xyz(i,ni,nj,coord,1)+
     #                             (2*alpha(i,1)*a(i,+1)-ni*a(i,-1))*
     #                                     b(i,0)*wt(point)
                           end if
                           if (jatom.ne.catom) then
                             xyz(i,ni,nj,coord,2)=xyz(i,ni,nj,coord,2)+
     #                             (2*alpha(i,2)*b(i,+1)-nj*b(i,-1))*
     #                                     a(i,0)*wt(point)
                           end if
   50                   continue
    7                continue
    8             continue
    9          continue
   10       continue
c
c     ----- multiply z 2-d integrals by exponential -----
c
            do 12 j=0,jmax
               do 11 i=0,imax
                  call vmul(xyz(1,i,j,3,0),xyz(1,i,j,3,0),
     #                            ryswt(1,root),n)
                  call vmul(xyz(1,i,j,3,1),xyz(1,i,j,3,1),
     #                            ryswt(1,root),n)
                  call vmul(xyz(1,i,j,3,2),xyz(1,i,j,3,2),
     #                            ryswt(1,root),n)
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
                  do 51 k=1,n
                     prmint(k,intgrl,1)=prmint(k,intgrl,1)+
     #                                    xyz(k,ix,jx,1,0)*
     #                                    xyz(k,iy,jy,2,0)*
     #                                    xyz(k,iz,jz,3,0)
                     prmint(k,intgrl,2)=prmint(k,intgrl,2)+
     #                                    xyz(k,ix,jx,1,1)*
     #                                    xyz(k,iy,jy,2,0)*
     #                                    xyz(k,iz,jz,3,0)
                     prmint(k,intgrl,3)=prmint(k,intgrl,3)+
     #                                    xyz(k,ix,jx,1,0)*
     #                                    xyz(k,iy,jy,2,1)*
     #                                    xyz(k,iz,jz,3,0)
                     prmint(k,intgrl,4)=prmint(k,intgrl,4)+
     #                                    xyz(k,ix,jx,1,0)*
     #                                    xyz(k,iy,jy,2,0)*
     #                                    xyz(k,iz,jz,3,1)
                     prmint(k,intgrl,5)=prmint(k,intgrl,5)+
     #                                    xyz(k,ix,jx,1,2)*
     #                                    xyz(k,iy,jy,2,0)*
     #                                    xyz(k,iz,jz,3,0)
                     prmint(k,intgrl,6)=prmint(k,intgrl,6)+
     #                                    xyz(k,ix,jx,1,0)*
     #                                    xyz(k,iy,jy,2,2)*
     #                                    xyz(k,iz,jz,3,0)
                     prmint(k,intgrl,7)=prmint(k,intgrl,7)+
     #                                    xyz(k,ix,jx,1,0)*
     #                                    xyz(k,iy,jy,2,0)*
     #                                    xyz(k,iz,jz,3,2)
   51             continue
   13          continue
   14       continue
   15    continue
c
c     ----- transform the primitive derivative integrals, and
c            then sum into the correct place in the final array
c
         do 100 junk=1,6
            if (junk.le.3) then
               coord=junk
               atom=iatom
            else
               coord=junk-3
               atom=jatom
            end if
c
            if (atom.eq.catom) go to 100
c
            call trans1(prmint(1,1,junk+1),conint,nprimi,nprimj,
     #                 nconti,ncontj,acf,bcf,tmp1,len1,lenblk,imin,
     #                 imax,jmin,jmax,nocart)
c
c           ----- put this angular momentum block in the right places
c
            numi=nobf(itype)
            numj=nobf(jtype)
            istart=start(iatom,itype)
            jstart=start(jatom,jtype)
            intgrl=0
c
            do 64 if=1,numi
               do 63 jf=1,numj
                  intgrl=intgrl+1
                  do 62 jc=1,ncontj
                     j=jstart+(jc-1)*numj+jf
                     do 61 ic=1,nconti
                        i=istart+(ic-1)*numi+if
                        if (i.ge.j) then
                           ij=i*(i-1)/2+j
                        else
                           if (iatom.eq.jatom.and.itype.eq.jtype)
     #                                   go to 61
                           ij=j*(j-1)/2+i
                        end if
c
                        ds(ij,coord,atom)=ds(ij,coord,atom)+
     #                                             conint(ic,jc,intgrl)
                        ds(ij,coord,catom)=ds(ij,coord,catom)-
     #                                             conint(ic,jc,intgrl)
   61                continue
   62             continue
   63          continue
   64       continue
c
  100    continue
c
   16 continue
c
c
c
      return
      end