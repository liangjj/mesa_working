*deck %W%  %G%
      subroutine fmdint(i2,d2,lenblk,rv,angmom,imax,jmax,
     $                  mmax,lmax,nroots,px,py,pz,lenxyz,mintyp,
     $                  nocart,nbtype,nv,energy,di,ndcen,der,leni)
c***begin prologue     fmdint.f
c***date written       840713   (yymmdd)
c***revision date      11/6/94  
c   11 november 1985     pws  at lanl
c          modifying original version (vfmint/m312) to sum up energy by directly
c          multiplying primitive integrals by two-particle density matrix.
c***keywords
c***author             saxe, paul    (lanl)
c***source             %W%   %G%
c***purpose
c***description
c***references         (none)
c
c***routines called    (none)
c
c***end prologue       fmdint.f
      implicit none
c     --- input variables -----
      integer lenblk,rv,imax,jmax,mmax,lmax,nroots
      integer lenxyz,nbtype,nv,ndcen,leni
c     --- input arrays (unmodified) ---
      integer px(lenxyz),py(lenxyz),pz(lenxyz)
      integer mintyp(nbtype),nocart(nbtype)
      integer angmom(4)
      real*8 i2(leni),d2(nv,lenblk)
      real*8 di(leni,ndcen)
c     --- input arrays (scratch) ---
c     --- output arrays ---
      real*8 der(3,ndcen)
c     --- output variables ---
      real*8 energy
c     --- scratch arrays ---
c     --- local variables ---
      integer nindex
      parameter (nindex=56)
      integer ipx(nindex),ipy(nindex),ipz(nindex)
      integer jpx(nindex),jpy(nindex),jpz(nindex)
      integer kpx(nindex),kpy(nindex),kpz(nindex)
      integer lpx(nindex),lpy(nindex),lpz(nindex)
      integer ix,iy,iz,jx,jy,jz,kx,ky,kz,lx,ly,lz
      integer mini,maxi,minj,maxj,mink,maxk,minl,maxl
      integer junk,i,j,k,l,n,m,prim,c
c     real*8 i2(rv,3,0:imax,0:jmax,0:mmax,0:lmax)
c
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
                  call plnkerr('bad number of rys roots',301)
c
c
    5             continue
                     do 50 prim=1,nv
                        energy=energy+i2(prim+lx)*i2(prim+ly)*
     $                                i2(prim+lz)*d2(prim,n)
   50                continue
                     do 55 c=1,ndcen
                        do 54 m=1,nv
                           der(1,c)=der(1,c)+d2(m,n)*
     $                         di(m+lx,c)*i2(m+ly)*i2(m+lz)
                           der(2,c)=der(2,c)+d2(m,n)*
     $                         i2(m+lx)*di(m+ly,c)*i2(m+lz)
                           der(3,c)=der(3,c)+d2(m,n)*
     $                         i2(m+lx)*i2(m+ly)*di(m+lz,c)
   54                   continue
   55                continue
                     go to 15
c
c
    6             continue
                     m=1
                     do 60 prim=1,nv
                        energy=energy+(i2(m+lx)*i2(m+ly)*i2(m+lz)+
     $                             i2(m+1+lx)*i2(m+1+ly)*i2(m+1+lz))*
     $                                     d2(prim,n)
                        m=m+2
   60                continue
                     do 65 c=1,ndcen
                        m=1
                        do 64 prim=1,nv
                           der(1,c)=der(1,c)+d2(prim,n)*
     $                        (di(m+lx  ,c)*i2(m+ly  )*i2(m+lz  )+
     $                         di(m+lx+1,c)*i2(m+ly+1)*i2(m+lz+1))
                           der(2,c)=der(2,c)+d2(prim,n)*
     $                        (i2(m+lx  )*di(m+ly  ,c)*i2(m+lz  )+
     $                         i2(m+lx+1)*di(m+ly+1,c)*i2(m+lz+1))
                           der(3,c)=der(3,c)+d2(prim,n)*
     $                        (i2(m+lx  )*i2(m+ly  )*di(m+lz  ,c)+
     $                         i2(m+lx+1)*i2(m+ly+1)*di(m+lz+1,c))
                           m=m+2
   64                   continue
   65                continue
                     go to 15
c
c
    7             continue
                     m=1
                     do 70 prim=1,nv
                        energy=energy+(i2(m+lx)*i2(m+ly)*i2(m+lz)+
     $                               i2(m+1+lx)*i2(m+1+ly)*i2(m+1+lz)+
     $                               i2(m+2+lx)*i2(m+2+ly)*i2(m+2+lz))*
     $                                    d2(prim,n)
                        m=m+3
   70                continue
                     do 75 c=1,ndcen
                        m=1
                        do 74 prim=1,nv
                           der(1,c)=der(1,c)+d2(prim,n)*
     $                        (di(m+lx  ,c)*i2(m+ly  )*i2(m+lz  )+
     $                         di(m+lx+1,c)*i2(m+ly+1)*i2(m+lz+1)+
     $                         di(m+lx+2,c)*i2(m+ly+2)*i2(m+lz+2))
                           der(2,c)=der(2,c)+d2(prim,n)*
     $                        (i2(m+lx  )*di(m+ly  ,c)*i2(m+lz  )+
     $                         i2(m+lx+1)*di(m+ly+1,c)*i2(m+lz+1)+
     $                         i2(m+lx+2)*di(m+ly+2,c)*i2(m+lz+2))
                           der(3,c)=der(3,c)+d2(prim,n)*
     $                        (i2(m+lx  )*i2(m+ly  )*di(m+lz  ,c)+
     $                         i2(m+lx+1)*i2(m+ly+1)*di(m+lz+1,c)+
     $                         i2(m+lx+2)*i2(m+ly+2)*di(m+lz+2,c))
                           m=m+3
   74                   continue
   75                continue
                     go to 15
c
c
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
                     do 85 c=1,ndcen
                        m=1
                        do 84 prim=1,nv
                           der(1,c)=der(1,c)+d2(prim,n)*
     $                        (di(m+lx  ,c)*i2(m+ly  )*i2(m+lz  )+
     $                         di(m+lx+1,c)*i2(m+ly+1)*i2(m+lz+1)+
     $                         di(m+lx+2,c)*i2(m+ly+2)*i2(m+lz+2)+
     $                         di(m+lx+3,c)*i2(m+ly+3)*i2(m+lz+3))
                           der(2,c)=der(2,c)+d2(prim,n)*
     $                        (i2(m+lx  )*di(m+ly  ,c)*i2(m+lz  )+
     $                         i2(m+lx+1)*di(m+ly+1,c)*i2(m+lz+1)+
     $                         i2(m+lx+2)*di(m+ly+2,c)*i2(m+lz+2)+
     $                         i2(m+lx+3)*di(m+ly+3,c)*i2(m+lz+3))
                           der(3,c)=der(3,c)+d2(prim,n)*
     $                        (i2(m+lx  )*i2(m+ly  )*di(m+lz  ,c)+
     $                         i2(m+lx+1)*i2(m+ly+1)*di(m+lz+1,c)+
     $                         i2(m+lx+2)*i2(m+ly+2)*di(m+lz+2,c)+
     $                         i2(m+lx+3)*i2(m+ly+3)*di(m+lz+3,c))
                           m=m+4
   84                   continue
   85                continue
                     go to 15
c
c
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
                     do 95 c=1,ndcen
                        m=1
                        do 94 prim=1,nv
                           der(1,c)=der(1,c)+d2(prim,n)*
     $                        (di(m+lx  ,c)*i2(m+ly  )*i2(m+lz  )+
     $                         di(m+lx+1,c)*i2(m+ly+1)*i2(m+lz+1)+
     $                         di(m+lx+2,c)*i2(m+ly+2)*i2(m+lz+2)+
     $                         di(m+lx+3,c)*i2(m+ly+3)*i2(m+lz+3)+
     $                         di(m+lx+4,c)*i2(m+ly+4)*i2(m+lz+4))
                           der(2,c)=der(2,c)+d2(prim,n)*
     $                        (i2(m+lx  )*di(m+ly  ,c)*i2(m+lz  )+
     $                         i2(m+lx+1)*di(m+ly+1,c)*i2(m+lz+1)+
     $                         i2(m+lx+2)*di(m+ly+2,c)*i2(m+lz+2)+
     $                         i2(m+lx+3)*di(m+ly+3,c)*i2(m+lz+3)+
     $                         i2(m+lx+4)*di(m+ly+4,c)*i2(m+lz+4))
                           der(3,c)=der(3,c)+d2(prim,n)*
     $                        (i2(m+lx  )*i2(m+ly  )*di(m+lz  ,c)+
     $                         i2(m+lx+1)*i2(m+ly+1)*di(m+lz+1,c)+
     $                         i2(m+lx+2)*i2(m+ly+2)*di(m+lz+2,c)+
     $                         i2(m+lx+3)*i2(m+ly+3)*di(m+lz+3,c)+
     $                         i2(m+lx+4)*i2(m+ly+4)*di(m+lz+4,c))
                           m=m+5
   94                   continue
   95                continue
                     go to 15
c
c
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
                     do 115 c=1,ndcen
                        m=1
                        do 114 prim=1,nv
                           der(1,c)=der(1,c)+d2(prim,n)*
     $                        (di(m+lx  ,c)*i2(m+ly  )*i2(m+lz  )+
     $                         di(m+lx+1,c)*i2(m+ly+1)*i2(m+lz+1)+
     $                         di(m+lx+2,c)*i2(m+ly+2)*i2(m+lz+2)+
     $                         di(m+lx+3,c)*i2(m+ly+3)*i2(m+lz+3)+
     $                         di(m+lx+4,c)*i2(m+ly+4)*i2(m+lz+4)+
     $                         di(m+lx+5,c)*i2(m+ly+5)*i2(m+lz+5))
                           der(2,c)=der(2,c)+d2(prim,n)*
     $                        (i2(m+lx  )*di(m+ly  ,c)*i2(m+lz  )+
     $                         i2(m+lx+1)*di(m+ly+1,c)*i2(m+lz+1)+
     $                         i2(m+lx+2)*di(m+ly+2,c)*i2(m+lz+2)+
     $                         i2(m+lx+3)*di(m+ly+3,c)*i2(m+lz+3)+
     $                         i2(m+lx+4)*di(m+ly+4,c)*i2(m+lz+4)+
     $                         i2(m+lx+5)*di(m+ly+5,c)*i2(m+lz+5))
                           der(3,c)=der(3,c)+d2(prim,n)*
     $                        (i2(m+lx  )*i2(m+ly  )*di(m+lz  ,c)+
     $                         i2(m+lx+1)*i2(m+ly+1)*di(m+lz+1,c)+
     $                         i2(m+lx+2)*i2(m+ly+2)*di(m+lz+2,c)+
     $                         i2(m+lx+3)*i2(m+ly+3)*di(m+lz+3,c)+
     $                         i2(m+lx+4)*i2(m+ly+4)*di(m+lz+4,c)+
     $                         i2(m+lx+5)*i2(m+ly+5)*di(m+lz+5,c))
                           m=m+6
  114                   continue
  115                continue
                     go to 15
c
c
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
                     do 125 c=1,ndcen
                        m=1
                        do 124 prim=1,nv
                           der(1,c)=der(1,c)+d2(prim,n)*
     $                        (di(m+lx  ,c)*i2(m+ly  )*i2(m+lz  )+
     $                         di(m+lx+1,c)*i2(m+ly+1)*i2(m+lz+1)+
     $                         di(m+lx+2,c)*i2(m+ly+2)*i2(m+lz+2)+
     $                         di(m+lx+3,c)*i2(m+ly+3)*i2(m+lz+3)+
     $                         di(m+lx+4,c)*i2(m+ly+4)*i2(m+lz+4)+
     $                         di(m+lx+5,c)*i2(m+ly+5)*i2(m+lz+5)+
     $                         di(m+lx+6,c)*i2(m+ly+6)*i2(m+lz+6))
                           der(2,c)=der(2,c)+d2(prim,n)*
     $                        (i2(m+lx  )*di(m+ly  ,c)*i2(m+lz  )+
     $                         i2(m+lx+1)*di(m+ly+1,c)*i2(m+lz+1)+
     $                         i2(m+lx+2)*di(m+ly+2,c)*i2(m+lz+2)+
     $                         i2(m+lx+3)*di(m+ly+3,c)*i2(m+lz+3)+
     $                         i2(m+lx+4)*di(m+ly+4,c)*i2(m+lz+4)+
     $                         i2(m+lx+5)*di(m+ly+5,c)*i2(m+lz+5)+
     $                         i2(m+lx+6)*di(m+ly+6,c)*i2(m+lz+6))
                           der(3,c)=der(3,c)+d2(prim,n)*
     $                        (i2(m+lx  )*i2(m+ly  )*di(m+lz  ,c)+
     $                         i2(m+lx+1)*i2(m+ly+1)*di(m+lz+1,c)+
     $                         i2(m+lx+2)*i2(m+ly+2)*di(m+lz+2,c)+
     $                         i2(m+lx+3)*i2(m+ly+3)*di(m+lz+3,c)+
     $                         i2(m+lx+4)*i2(m+ly+4)*di(m+lz+4,c)+
     $                         i2(m+lx+5)*i2(m+ly+5)*di(m+lz+5,c)+
     $                         i2(m+lx+6)*i2(m+ly+6)*di(m+lz+6,c))
                           m=m+7
  124                   continue
  125                continue
                     go to 15
c
c
 12                  continue
                     m=1
                     do 130 prim=1,nv
                        energy=energy+(i2(m+lx)*i2(m+ly)*i2(m+lz)+
     $                               i2(m+1+lx)*i2(m+1+ly)*i2(m+1+lz)+
     $                               i2(m+2+lx)*i2(m+2+ly)*i2(m+2+lz)+
     $                               i2(m+3+lx)*i2(m+3+ly)*i2(m+3+lz)+
     $                               i2(m+4+lx)*i2(m+4+ly)*i2(m+4+lz)+
     $                               i2(m+5+lx)*i2(m+5+ly)*i2(m+5+lz)+
     $                               i2(m+6+lx)*i2(m+6+ly)*i2(m+6+lz)+
     $                               i2(m+7+lx)*i2(m+7+ly)*i2(m+7+lz))*
     $                                     d2(prim,n)
                        m=m+8
  130                continue
                     do 135 c=1,ndcen
                        m=1
                        do 134 prim=1,nv
                           der(1,c)=der(1,c)+d2(prim,n)*
     $                        (di(m+lx  ,c)*i2(m+ly  )*i2(m+lz  )+
     $                         di(m+lx+1,c)*i2(m+ly+1)*i2(m+lz+1)+
     $                         di(m+lx+2,c)*i2(m+ly+2)*i2(m+lz+2)+
     $                         di(m+lx+3,c)*i2(m+ly+3)*i2(m+lz+3)+
     $                         di(m+lx+4,c)*i2(m+ly+4)*i2(m+lz+4)+
     $                         di(m+lx+5,c)*i2(m+ly+5)*i2(m+lz+5)+
     $                         di(m+lx+6,c)*i2(m+ly+6)*i2(m+lz+6)+
     $                         di(m+lx+7,c)*i2(m+ly+7)*i2(m+lz+7))
                           der(2,c)=der(2,c)+d2(prim,n)*
     $                        (i2(m+lx  )*di(m+ly  ,c)*i2(m+lz  )+
     $                         i2(m+lx+1)*di(m+ly+1,c)*i2(m+lz+1)+
     $                         i2(m+lx+2)*di(m+ly+2,c)*i2(m+lz+2)+
     $                         i2(m+lx+3)*di(m+ly+3,c)*i2(m+lz+3)+
     $                         i2(m+lx+4)*di(m+ly+4,c)*i2(m+lz+4)+
     $                         i2(m+lx+5)*di(m+ly+5,c)*i2(m+lz+5)+
     $                         i2(m+lx+6)*di(m+ly+6,c)*i2(m+lz+6)+
     $                         i2(m+lx+7)*di(m+ly+7,c)*i2(m+lz+7))
                           der(3,c)=der(3,c)+d2(prim,n)*
     $                        (i2(m+lx  )*i2(m+ly  )*di(m+lz  ,c)+
     $                         i2(m+lx+1)*i2(m+ly+1)*di(m+lz+1,c)+
     $                         i2(m+lx+2)*i2(m+ly+2)*di(m+lz+2,c)+
     $                         i2(m+lx+3)*i2(m+ly+3)*di(m+lz+3,c)+
     $                         i2(m+lx+4)*i2(m+ly+4)*di(m+lz+4,c)+
     $                         i2(m+lx+5)*i2(m+ly+5)*di(m+lz+5,c)+
     $                         i2(m+lx+6)*i2(m+ly+6)*di(m+lz+6,c)+
     $                         i2(m+lx+7)*i2(m+ly+7)*di(m+lz+7,c))
                           m=m+8
  134                   continue
  135                continue
                     go to 15
c
c
 13                  continue
                     m=1
                     do 140 prim=1,nv
                        energy=energy+(i2(m+lx)*i2(m+ly)*i2(m+lz)+
     $                               i2(m+1+lx)*i2(m+1+ly)*i2(m+1+lz)+
     $                               i2(m+2+lx)*i2(m+2+ly)*i2(m+2+lz)+
     $                               i2(m+3+lx)*i2(m+3+ly)*i2(m+3+lz)+
     $                               i2(m+4+lx)*i2(m+4+ly)*i2(m+4+lz)+
     $                               i2(m+5+lx)*i2(m+5+ly)*i2(m+5+lz)+
     $                               i2(m+6+lx)*i2(m+6+ly)*i2(m+6+lz)+
     $                               i2(m+7+lx)*i2(m+7+ly)*i2(m+7+lz)+
     $                               i2(m+8+lx)*i2(m+8+ly)*i2(m+8+lz))*
     $                                     d2(prim,n)
                        m=m+9
  140                continue
                     do 145 c=1,ndcen
                        m=1
                        do 144 prim=1,nv
                           der(1,c)=der(1,c)+d2(prim,n)*
     $                        (di(m+lx  ,c)*i2(m+ly  )*i2(m+lz  )+
     $                         di(m+lx+1,c)*i2(m+ly+1)*i2(m+lz+1)+
     $                         di(m+lx+2,c)*i2(m+ly+2)*i2(m+lz+2)+
     $                         di(m+lx+3,c)*i2(m+ly+3)*i2(m+lz+3)+
     $                         di(m+lx+4,c)*i2(m+ly+4)*i2(m+lz+4)+
     $                         di(m+lx+5,c)*i2(m+ly+5)*i2(m+lz+5)+
     $                         di(m+lx+6,c)*i2(m+ly+6)*i2(m+lz+6)+
     $                         di(m+lx+7,c)*i2(m+ly+7)*i2(m+lz+7)+
     $                         di(m+lx+8,c)*i2(m+ly+8)*i2(m+lz+8))
                           der(2,c)=der(2,c)+d2(prim,n)*
     $                        (i2(m+lx  )*di(m+ly  ,c)*i2(m+lz  )+
     $                         i2(m+lx+1)*di(m+ly+1,c)*i2(m+lz+1)+
     $                         i2(m+lx+2)*di(m+ly+2,c)*i2(m+lz+2)+
     $                         i2(m+lx+3)*di(m+ly+3,c)*i2(m+lz+3)+
     $                         i2(m+lx+4)*di(m+ly+4,c)*i2(m+lz+4)+
     $                         i2(m+lx+5)*di(m+ly+5,c)*i2(m+lz+5)+
     $                         i2(m+lx+6)*di(m+ly+6,c)*i2(m+lz+6)+
     $                         i2(m+lx+7)*di(m+ly+7,c)*i2(m+lz+7)+
     $                         i2(m+lx+8)*di(m+ly+8,c)*i2(m+lz+8))
                           der(3,c)=der(3,c)+d2(prim,n)*
     $                        (i2(m+lx  )*i2(m+ly  )*di(m+lz  ,c)+
     $                         i2(m+lx+1)*i2(m+ly+1)*di(m+lz+1,c)+
     $                         i2(m+lx+2)*i2(m+ly+2)*di(m+lz+2,c)+
     $                         i2(m+lx+3)*i2(m+ly+3)*di(m+lz+3,c)+
     $                         i2(m+lx+4)*i2(m+ly+4)*di(m+lz+4,c)+
     $                         i2(m+lx+5)*i2(m+ly+5)*di(m+lz+5,c)+
     $                         i2(m+lx+6)*i2(m+ly+6)*di(m+lz+6,c)+
     $                         i2(m+lx+7)*i2(m+ly+7)*di(m+lz+7,c)+
     $                         i2(m+lx+8)*i2(m+ly+8)*di(m+lz+8,c))
                           m=m+9
  144                   continue
  145                continue
                     go to 15
c
c
   15             continue

  100          continue
  200       continue
  300    continue
  400 continue
c
c
      return
      end
