*deck @(#)vfmint.f	5.1  11/6/94
      subroutine vfmint(i2,d2,
     $                  lenblk,rv,angmom,imax,jmax,
     $                  mmax,lmax,nroots,px,py,pz,lenxyz,mintyp,
     $                  nocart,nbtype,nv,energy)
c***begin prologue     vfmint.f
c***date written       840713   (yymmdd)
c***revision date      11/6/94
c   june 3, 1993       rlm at lanl
c          adding loops to handle higher angular momentum. this choked on
c          second derivatives of f functions.
c   11 november 1985     pws  at lanl
c          modifying original version from m312 to sum up energy by directly
c          multiplying primitive integrals by two-particle density matrix.
c***keywords
c***author             saxe, paul    (lanl)
c***source             @(#)vfmint.f	5.1   11/6/94
c***purpose
c***description
c***references         (none)
c
c***routines called    (none)
c
c***end prologue       vfmint.f
      implicit none
c     --- input variables -----
      integer lenblk,lenxyz,nbtype
      integer imax,jmax,mmax,lmax,rv,nv
      integer nroots
c     --- input arrays (unmodified) ---
      integer angmom(4)
      integer px(lenxyz),py(lenxyz),pz(lenxyz)
      integer mintyp(nbtype),nocart(nbtype)
      real*8 i2(*),d2(nv,lenblk)
c     --- input arrays (scratch) ---
c     --- output arrays ---
c     --- output variables ---
      real*8 energy
c     --- scratch arrays ---
c     --- local variables ---
      integer inp,iout
      integer nindex
      integer mini,maxi,minj,maxj,mink,maxk,minl,maxl
      integer i,j,k,l,m,n,ix,iy,iz,jx,jy,jz,kx,ky,kz,lx,ly,lz
      integer junk,prim
c
      parameter (nindex=56)
      integer ipx(nindex),ipy(nindex),ipz(nindex)
      integer jpx(nindex),jpy(nindex),jpz(nindex)
      integer kpx(nindex),kpy(nindex),kpz(nindex)
      integer lpx(nindex),lpy(nindex),lpz(nindex)
c
      common/io/inp,iout
c
c     --- set up minimum and maximum functions in each momentum group ---
      mini=mintyp(angmom(1))
      maxi=mini+nocart(angmom(1))-1
      minj=mintyp(angmom(2))
      maxj=minj+nocart(angmom(2))-1
      mink=mintyp(angmom(3))
      maxk=mink+nocart(angmom(3))-1
      minl=mintyp(angmom(4))
      maxl=minl+nocart(angmom(4))-1
c
c     --- set up offsets ---
      junk=3*rv
      do 1 i=mini,maxi
         ipx(i)=     junk*px(i)
         ipy(i)=  rv+junk*py(i)
         ipz(i)=2*rv+junk*pz(i)
    1 continue
      junk=junk*(imax+1)
      do 2 j=minj,maxj
         jpx(j)=junk*px(j)
         jpy(j)=junk*py(j)
         jpz(j)=junk*pz(j)
    2 continue
      junk=junk*(jmax+1)
      do 3 k=mink,maxk
         kpx(k)=junk*px(k)
         kpy(k)=junk*py(k)
         kpz(k)=junk*pz(k)
    3 continue
      junk=junk*(mmax+1)
      do 4 l=minl,maxl
         lpx(l)=junk*px(l)
         lpy(l)=junk*py(l)
         lpz(l)=junk*pz(l)
    4 continue
c
      n=0
      do 400 i=mini,maxi
         ix=ipx(i)
         iy=ipy(i)
         iz=ipz(i)
         do 300 j=minj,maxj
            jx=ix+jpx(j)
            jy=iy+jpy(j)
            jz=iz+jpz(j)
c
            do 200 k=mink,maxk
               kx=jx+kpx(k)
               ky=jy+kpy(k)
               kz=jz+kpz(k)
               do 100 l=minl,maxl
                  lx=kx+lpx(l)
                  ly=ky+lpy(l)
                  lz=kz+lpz(l)
c
                  n=n+1
c
                  go to (5,6,7,8,9,10,11,12,13),nroots
                  write(iout,*) 'vfmint:712,nroots',nroots
                  call lnkerr('bad number of rys roots')
c
    5             continue
                     do 50 prim=1,nv
                        energy=energy+i2(prim+lx)*i2(prim+ly)*
     $                                i2(prim+lz)*d2(prim,n)
   50                continue
                     go to 15
    6             continue
                     m=1
                     do 60 prim=1,nv
                        energy=energy+(i2(m+lx)*i2(m+ly)*i2(m+lz)+
     $                               i2(m+1+lx)*i2(m+1+ly)*i2(m+1+lz))*
     $                                     d2(prim,n)
                        m=m+2
   60                continue
                     go to 15
    7             continue
                     m=1
                     do 70 prim=1,nv
                        energy=energy+(i2(m+lx)*i2(m+ly)*i2(m+lz)+
     $                               i2(m+1+lx)*i2(m+1+ly)*i2(m+1+lz)+
     $                               i2(m+2+lx)*i2(m+2+ly)*i2(m+2+lz))*
     $                                     d2(prim,n)
                        m=m+3
   70                continue
                     go to 15
    8             continue
                     m=1
                     do 80 prim=1,nv
                        energy=energy+(i2(m+lx)*i2(m+ly)*i2(m+lz)+
     $                               i2(m+1+lx)*i2(m+1+ly)*i2(m+1+lz)+
     $                               i2(m+2+lx)*i2(m+2+ly)*i2(m+2+lz)+
     $                               i2(m+3+lx)*i2(m+3+ly)*i2(m+3+lz))*
     $                                     d2(prim,n)
                        m=m+4
   80                continue
                     go to 15
    9             continue
                     m=1
                     do 90 prim=1,nv
                        energy=energy+(i2(m+lx)*i2(m+ly)*i2(m+lz)+
     $                               i2(m+1+lx)*i2(m+1+ly)*i2(m+1+lz)+
     $                               i2(m+2+lx)*i2(m+2+ly)*i2(m+2+lz)+
     $                               i2(m+3+lx)*i2(m+3+ly)*i2(m+3+lz)+
     $                               i2(m+4+lx)*i2(m+4+ly)*i2(m+4+lz))*
     $                                     d2(prim,n)
                        m=m+5
   90                continue
                     go to 15
   10             continue
                     m=1
                     do 110 prim=1,nv
                        energy=energy+(i2(m+lx)*i2(m+ly)*i2(m+lz)+
     $                               i2(m+1+lx)*i2(m+1+ly)*i2(m+1+lz)+
     $                               i2(m+2+lx)*i2(m+2+ly)*i2(m+2+lz)+
     $                               i2(m+3+lx)*i2(m+3+ly)*i2(m+3+lz)+
     $                               i2(m+4+lx)*i2(m+4+ly)*i2(m+4+lz)+
     $                               i2(m+5+lx)*i2(m+5+ly)*i2(m+5+lz))*
     $                                     d2(prim,n)
                        m=m+6
  110                continue
                     go to 15
   11             continue
                     m=1
                     do 120 prim=1,nv
                        energy=energy+(i2(m+lx)*i2(m+ly)*i2(m+lz)+
     $                               i2(m+1+lx)*i2(m+1+ly)*i2(m+1+lz)+
     $                               i2(m+2+lx)*i2(m+2+ly)*i2(m+2+lz)+
     $                               i2(m+3+lx)*i2(m+3+ly)*i2(m+3+lz)+
     $                               i2(m+4+lx)*i2(m+4+ly)*i2(m+4+lz)+
     $                               i2(m+5+lx)*i2(m+5+ly)*i2(m+5+lz)+
     $                               i2(m+6+lx)*i2(m+6+ly)*i2(m+6+lz))*
     $                                     d2(prim,n)
                        m=m+7
  120                continue
                     go to 15
   12                continue
                     m=1
                     do 130 prim=1,nv
                        energy=energy+(i2(m+lx)*i2(m+ly)*i2(m+lz)+
     $                                i2(m+1+lx)*i2(m+1+ly)*i2(m+1+lz)+
     $                                i2(m+2+lx)*i2(m+2+ly)*i2(m+2+lz)+
     $                                i2(m+3+lx)*i2(m+3+ly)*i2(m+3+lz)+
     $                                i2(m+4+lx)*i2(m+4+ly)*i2(m+4+lz)+
     $                                i2(m+5+lx)*i2(m+5+ly)*i2(m+5+lz)+
     $                                i2(m+6+lx)*i2(m+6+ly)*i2(m+6+lz)+
     $                                i2(m+7+lx)*i2(m+7+ly)*i2(m+7+lz))*
     $                                d2(prim,n) 
                        m=m+8
  130                continue
                     go to 15
   13                continue
                     m=1
                     do 140 prim=1,nv
                        energy=energy+(i2(m+lx)*i2(m+ly)*i2(m+lz)+
     $                                i2(m+1+lx)*i2(m+1+ly)*i2(m+1+lz)+
     $                                i2(m+2+lx)*i2(m+2+ly)*i2(m+2+lz)+
     $                                i2(m+3+lx)*i2(m+3+ly)*i2(m+3+lz)+
     $                                i2(m+4+lx)*i2(m+4+ly)*i2(m+4+lz)+
     $                                i2(m+5+lx)*i2(m+5+ly)*i2(m+5+lz)+
     $                                i2(m+6+lx)*i2(m+6+ly)*i2(m+6+lz)+
     $                                i2(m+7+lx)*i2(m+7+ly)*i2(m+7+lz)+
     $                                i2(m+8+lx)*i2(m+8+ly)*i2(m+8+lz))*
     $                                d2(prim,n)
                        m=m+9
  140                continue
   15             continue
c
  100          continue
  200       continue
  300    continue
  400 continue
c
c
      return
      end
