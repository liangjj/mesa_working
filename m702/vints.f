*deck @(#)vints.f	5.1  11/6/94
      subroutine vints(iatom,jatom,imax,jmax,nprimi,nprimj,
     #                 i1,j1,ex,alpha,c,ainv,xyza,expon,xyz,
     #                 nv,nprim,nat,nbtype,h,wt,mxpts,a,b,
     #                 aiaj,xyz0,rysrt,ryswt,nroots,prmint,
     #                 lenblk,zan,xyz1,t1,mini,maxi,minj,maxj,
     #                 nx,ny,nz,nderiv,conint,nconti,ncontj,
     #                 itype,jtype,acf,bcf,tmp1,len1,imin,jmin,nocart,
     #                 d,nnp,start,nobf,grad,nxyz,ld2e,nd2e)
c
c***module to form the two-dimensional integrals over primitives
c   for the potential-energy integrals
c
c paul saxe                26 july 1984                  lanl
c
      implicit integer (a-z)
c
      real*8 d2e(45)
      real*8 ex(nprim),alpha(nv,2),c(3,nat),ainv(nv),xyza(nv,3)
      real*8 expon(nv),xyz(nv,0:imax,0:jmax,3,nxyz)
      real*8 ld2e(nd2e)
      real*8 a(nv,-nderiv:nderiv),b(nv,-nderiv:nderiv),aiaj(nv)
      real*8 ryswt(nv,nroots),prmint(nv,lenblk,*),zan(nat)
      real*8 xyz1(nv,3),t1(nv),h(mxpts),wt(mxpts),rysrt(nv,nroots)
      real*8 xyz0(nv,3)
      real*8 conint(nconti,ncontj,lenblk),tmp1(len1),d(nnp)
      real*8 grad(3,nat)
      real*8 acf(nprimi,nprimj,imin:imax),bcf(nprimj,ncontj,jmin:jmax)
      real*8 iex,jex,scalar,one,temprt(9),tempwt(9),pi212,xx
      integer nx(*),ny(*),nz(*)
      integer nocart(nbtype),start(nat,nbtype),nobf(nbtype)
      integer centre(3)
      integer idx(28),idy(28),idz(28)
c
      common/io/inp,iout
c
      parameter (pi212=1.1283791670955d+00, one=1.0d+00)
c
c   derivative order defined by the following arrays:
c                0
c          x y z x'y'z'
c          xx xy yy xz yz zz xx' yx' zx' x'x' xy' yy' zy' x'y'
c              y'y' xz' yz' zz' x'z' y'z' z'z'
c
      data idx / 1,
     $     2,1,1,3,1,1,
     $     4,2,1,2,1,1,5,3,3,6,2,1,1,3,1,2,1,1,3,1,1/
      data idy / 1,
     $     1,2,1,1,3,1,
     $     1,2,4,1,2,1,1,2,1,1,3,5,3,3,6,1,2,1,1,3,1/
      data idz / 1,
     $     1,1,2,1,1,3,
     $     1,1,1,2,2,4,1,1,2,1,1,1,2,1,1,3,3,5,3,3,6/
      save idx,idy,idz
c
      centre(1)=iatom
      centre(2)=jatom
c
c     ----- how mny integrals derivatives to form -----
c
      if (nderiv.eq.0) then
         ndxyz=1
      else if (nderiv.eq.1) then
         ndxyz=1+6
      else if (nderiv.eq.2) then
         ndxyz=1+6+21
      end if
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
   99    format (//,' ##### m702: vints, error, n and nv:',2i8,//)
         call lnkerr('m702: vints error in vector lengths')
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
      call rzero(prmint,n*lenblk*ndxyz)
c
c     ----- form two-dimensional integrals -----
c
      do 16 catom=1,nat
         centre(3)=catom
c
c     ----- calculate roots and weights of rys polynomials -----
c
         call rzero(prmint(1,1,2),n*lenblk*(ndxyz-1))
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
            call rzero(xyz,n*(imax+1)*(jmax+1)*3*nxyz)
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
                        if (nderiv.eq.2) then
                           if (ni.gt.2) then
                              call vpower(a(1,-2),a(1,0),ni-2,n)
                              call vpower(a(1,-1),a(1,0),ni-1,n)
                              call vpower(a(1, 1),a(1,0),ni+1,n)
                              call vpower(a(1, 2),a(1,0),ni+2,n)
                              call vpower(a(1, 0),a(1,0),ni  ,n)
                           else if (ni.eq.2) then
                              call vfill(a(1,-2),one,n)
                              call vpower(a(1,-1),a(1,0),ni-1,n)
                              call vpower(a(1, 1),a(1,0),3   ,n)
                              call vpower(a(1, 2),a(1,0),4   ,n)
                              call vpower(a(1, 0),a(1,0),2   ,n)
                           else if (ni.eq.1) then
                              call rzero(a(1,-2),n)
                              call vfill(a(1,-1),one,n)
                              call vpower(a(1,1),a(1,0),2,n)
                              call vpower(a(1,2),a(1,0),3,n)
                           else if (ni.eq.0) then
                              call rzero(a(1,-2),n)
                              call rzero(a(1,-1),n)
                              call vmove(a(1,1),a(1,0),n)
                              call vpower(a(1,2),a(1,0),2,n)
                              call vfill(a(1,0),one,n)
                           end if
                        else if (nderiv.eq.1) then
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
                        if (nderiv.eq.2) then
                           if (nj.gt.2) then
                              call vpower(b(1,-2),b(1,0),nj-2,n)
                              call vpower(b(1,-1),b(1,0),nj-1,n)
                              call vpower(b(1, 1),b(1,0),nj+1,n)
                              call vpower(b(1, 2),b(1,0),nj+2,n)
                              call vpower(b(1, 0),b(1,0),nj  ,n)
                           else if (nj.eq.2) then
                              call vfill(b(1,-2),one,n)
                              call vpower(b(1,-1),b(1,0),nj-1,n)
                              call vpower(b(1, 1),b(1,0),3   ,n)
                              call vpower(b(1, 2),b(1,0),4   ,n)
                              call vpower(b(1, 0),b(1,0),2   ,n)
                           else if (nj.eq.1) then
                              call rzero(b(1,-2),n)
                              call vfill(b(1,-1),one,n)
                              call vpower(b(1,1),b(1,0),2,n)
                              call vpower(b(1,2),b(1,0),3,n)
                           else if (nj.eq.0) then
                              call rzero(b(1,-2),n)
                              call rzero(b(1,-1),n)
                              call vmove(b(1,1),b(1,0),n)
                              call vpower(b(1,2),b(1,0),2,n)
                              call vfill(b(1,0),one,n)
                           end if
                        else if (nderiv.eq.1) then
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
                           xyz(i,ni,nj,coord,1)=xyz(i,ni,nj,coord,1)+
     #                          a(i,0)*b(i,0)*wt(point)
                           if (iatom.ne.catom) then
                              xyz(i,ni,nj,coord,2)=
     $                             xyz(i,ni,nj,coord,2)+
     #                             (2*alpha(i,1)*a(i,+1)-ni*a(i,-1))*
     #                             b(i,0)*wt(point)
                           end if
                           if (jatom.ne.catom) then
                              xyz(i,ni,nj,coord,3)=
     $                             xyz(i,ni,nj,coord,3)+
     #                             (2*alpha(i,2)*b(i,+1)-nj*b(i,-1))*
     #                             a(i,0)*wt(point)
                           end if
                           if (nderiv.ge.2) then
                              if (iatom.ne.catom) then
                                 xyz(i,ni,nj,coord,4)=
     $                                xyz(i,ni,nj,coord,4)+
     $                                (4*alpha(i,1)**2*a(i,+2)-
     $                                2*alpha(i,1)*(2*ni+1)*a(i,0)+
     $                                ni*(ni-1)*a(i,-2))*
     $                                b(i,0)*wt(point)
                              end if
                              if (iatom.ne.catom.and.jatom.ne.catom)
     $                             then
                                 xyz(i,ni,nj,coord,5)=
     $                                xyz(i,ni,nj,coord,5)+
     $                                (2*alpha(i,1)*a(i,+1)-ni*a(i,-1))*
     $                                (2*alpha(i,2)*b(i,+1)-nj*b(i,-1))*
     $                                wt(point)
                              end if
                              if (jatom.ne.catom) then
                                 xyz(i,ni,nj,coord,6)=
     $                                xyz(i,ni,nj,coord,6)+
     $                                (4*alpha(i,2)**2*b(i,+2)-
     $                                2*alpha(i,2)*(2*nj+1)*b(i,0)+
     $                                nj*(nj-1)*b(i,-2))*
     $                                a(i,0)*wt(point)
                              end if
                           end if
   50                   continue
    7                continue
    8             continue
    9          continue
   10       continue
c
c     ----- multiply z 2-d integrals by exponential -----
c
            do 113 ixyz=1,nxyz
               do 12 j=0,jmax
                  do 11 i=0,imax
                     call vmul(xyz(1,i,j,3,ixyz),xyz(1,i,j,3,ixyz),
     #                    ryswt(1,root),n)
 11               continue
 12            continue
 113         continue
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
                  do 150 dxyz=1,ndxyz
                     do 51 k=1,n
                        prmint(k,intgrl,dxyz)=prmint(k,intgrl,dxyz)+
     $                       xyz(k,ix,jx,1,idx(dxyz))*
     $                       xyz(k,iy,jy,2,idy(dxyz))*
     $                       xyz(k,iz,jz,3,idz(dxyz))
 51                  continue
 150              continue
c
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
                        grad(coord,atom)=grad(coord,atom)+
     #                           conint(ic,jc,intgrl)*d(ij)
                        grad(coord,catom)=grad(coord,catom)-
     #                           conint(ic,jc,intgrl)*d(ij)
   61                continue
   62             continue
   63          continue
   64       continue
c
  100    continue
c
c        ----- second derivatives -----
c
         if (nderiv.ge.2) then
            call rzero(d2e,45)
            do 200 junk=1,21
c
cps               if (atom.eq.catom) go to 200
c
               call trans1(prmint(1,1,junk+7),conint,nprimi,nprimj,
     #              nconti,ncontj,acf,bcf,tmp1,len1,lenblk,imin,
     #              imax,jmin,jmax,nocart)
c
c           ----- put this angular momentum block in the right places
c
               numi=nobf(itype)
               numj=nobf(jtype)
               istart=start(iatom,itype)
               jstart=start(jatom,jtype)
               intgrl=0
c
               do 164 if=1,numi
                  do 163 jf=1,numj
                     intgrl=intgrl+1
                     do 162 jc=1,ncontj
                        j=jstart+(jc-1)*numj+jf
                        do 161 ic=1,nconti
                           i=istart+(ic-1)*numi+if
                           if (i.ge.j) then
                              ij=i*(i-1)/2+j
                           else
                              if (iatom.eq.jatom.and.itype.eq.jtype)
     #                             go to 161
                              ij=j*(j-1)/2+i
                           end if
c
                           d2e(junk)=d2e(junk)+conint(ic,jc,intgrl)*
     $                          d(ij)
 161                    continue
 162                 continue
 163              continue
 164           continue
c
 200        continue
c
c           ---- and fill out using translational invariance -----
c
            call movd2e(d2e,ld2e,nd2e,centre,2,2)
         end if
c
   16 continue
c
c
c
      return
      end
