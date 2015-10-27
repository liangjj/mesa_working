*deck @(#)fd2ant.f	5.1  11/28/95
      subroutine fd2ant(i2,dens,lenblk,rv,angmom,imax,jmax,
     $     mmax,lmax,nroots,px,py,pz,lenxyz,mintyp,
     $     nocart,nbtype,nv,energy,di,ndcen,der,leni,d2i,nd2,d2,
     $     npass)
c***begin prologue     fd2ant.f
c***date written       840713   (yymmdd)
c***revision date      11/6/94
c   11 november 1985     pws  at lanl
c          modifying original version from m312 to sum up energy by directly
c          multiplying primitive integrals by two-particle density matrix.
c***keywords
c***author             saxe, paul    (lanl)
c***source             @(#)fd2ant.f	5.1   11/28/95
c***purpose
c***description
c***references         (none)
c
c***routines called    (none)
c
c***end prologue       fd2ant.f
      implicit none
c     --- input variables -----
      integer lenblk,rv,imax,jmax,mmax,lmax,nroots
      integer nbtype,lenxyz,nv,ndcen,leni,nd2,npass
c     --- input arrays (unmodified) ---
      integer angmom(4)
      integer px(lenxyz),py(lenxyz),pz(lenxyz)
      integer mintyp(nbtype),nocart(nbtype)
      real*8 i2(leni)
      real*8 dens(nv,lenblk)
      real*8 di(leni,ndcen)
      real*8 d2i(leni,nd2)
c     --- input arrays (scratch) ---
c     --- output arrays ---
      real*8 der(3,ndcen)
      real*8 d2(78)
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
c
      integer mini,maxi,minj,maxj,mink,maxk,minl,maxl
      integer junk,i,j,k,l,m,n,prim
      integer ix,iy,iz,jx,jy,jz,kx,ky,kz,lx,ly,lz
c     real*8 i2(rv,3,0:imax,0:jmax,0:mmax,0:lmax)
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
                  call plnkerr('bad number of rys roots',501)
c
c
    5             continue
                  call plnkerr(' fd2ant: error ',502)
                     go to 15
    6             continue
                  call plnkerr(' fd2ant: error ',503)
                     go to 15
    7             continue
                  call plnkerr(' fd2ant: error ',504)
                     go to 15
    8             continue
                  call plnkerr(' fd2ant: error ',505)
                     go to 15
    9             continue
                  call plnkerr(' fd2ant: error ',506)
                     go to 15
   10             continue
                     m=1
                     do 1000 prim=1,nv
                        energy=energy+dens(prim,n)*
     $                      (i2(m+lx  )*i2(m+ly  )*i2(m+lz  )+
     $                       i2(m+lx+1)*i2(m+ly+1)*i2(m+lz+1)+
     $                       i2(m+lx+2)*i2(m+ly+2)*i2(m+lz+2)+
     $                       i2(m+lx+3)*i2(m+ly+3)*i2(m+lz+3)+
     $                       i2(m+lx+4)*i2(m+ly+4)*i2(m+lz+4)+
     $                       i2(m+lx+5)*i2(m+ly+5)*i2(m+lz+5))
                        der(1,1)=der(1,1)+dens(prim,n)*
     $                      (di(m+lx  ,1)*i2(m+ly  )*i2(m+lz  )+
     $                       di(m+lx+1,1)*i2(m+ly+1)*i2(m+lz+1)+
     $                       di(m+lx+2,1)*i2(m+ly+2)*i2(m+lz+2)+
     $                       di(m+lx+3,1)*i2(m+ly+3)*i2(m+lz+3)+
     $                       di(m+lx+4,1)*i2(m+ly+4)*i2(m+lz+4)+
     $                       di(m+lx+5,1)*i2(m+ly+5)*i2(m+lz+5))
                        der(2,1)=der(2,1)+dens(prim,n)*
     $                      (i2(m+lx  )*di(m+ly  ,1)*i2(m+lz  )+
     $                       i2(m+lx+1)*di(m+ly+1,1)*i2(m+lz+1)+
     $                       i2(m+lx+2)*di(m+ly+2,1)*i2(m+lz+2)+
     $                       i2(m+lx+3)*di(m+ly+3,1)*i2(m+lz+3)+
     $                       i2(m+lx+4)*di(m+ly+4,1)*i2(m+lz+4)+
     $                       i2(m+lx+5)*di(m+ly+5,1)*i2(m+lz+5))
                        der(3,1)=der(3,1)+dens(prim,n)*
     $                      (i2(m+lx  )*i2(m+ly  )*di(m+lz  ,1)+
     $                       i2(m+lx+1)*i2(m+ly+1)*di(m+lz+1,1)+
     $                       i2(m+lx+2)*i2(m+ly+2)*di(m+lz+2,1)+
     $                       i2(m+lx+3)*i2(m+ly+3)*di(m+lz+3,1)+
     $                       i2(m+lx+4)*i2(m+ly+4)*di(m+lz+4,1)+
     $                       i2(m+lx+5)*i2(m+ly+5)*di(m+lz+5,1))
                        d2(1)=d2(1)+dens(prim,n)*
     $                      (d2i(m+lx  ,1)*i2(m+ly  )*i2(m+lz  )+
     $                       d2i(m+lx+1,1)*i2(m+ly+1)*i2(m+lz+1)+
     $                       d2i(m+lx+2,1)*i2(m+ly+2)*i2(m+lz+2)+
     $                       d2i(m+lx+3,1)*i2(m+ly+3)*i2(m+lz+3)+
     $                       d2i(m+lx+4,1)*i2(m+ly+4)*i2(m+lz+4)+
     $                       d2i(m+lx+5,1)*i2(m+ly+5)*i2(m+lz+5))
                        d2(2)=d2(2)+dens(prim,n)*
     $                      (di(m+lx  ,1)*di(m+ly  ,1)*i2(m+lz  )+
     $                       di(m+lx+1,1)*di(m+ly+1,1)*i2(m+lz+1)+
     $                       di(m+lx+2,1)*di(m+ly+2,1)*i2(m+lz+2)+
     $                       di(m+lx+3,1)*di(m+ly+3,1)*i2(m+lz+3)+
     $                       di(m+lx+4,1)*di(m+ly+4,1)*i2(m+lz+4)+
     $                       di(m+lx+5,1)*di(m+ly+5,1)*i2(m+lz+5))
                        d2(3)=d2(3)+dens(prim,n)*
     $                      (i2(m+lx  )*d2i(m+ly  ,1)*i2(m+lz  )+
     $                       i2(m+lx+1)*d2i(m+ly+1,1)*i2(m+lz+1)+
     $                       i2(m+lx+2)*d2i(m+ly+2,1)*i2(m+lz+2)+
     $                       i2(m+lx+3)*d2i(m+ly+3,1)*i2(m+lz+3)+
     $                       i2(m+lx+4)*d2i(m+ly+4,1)*i2(m+lz+4)+
     $                       i2(m+lx+5)*d2i(m+ly+5,1)*i2(m+lz+5))
                        d2(4)=d2(4)+dens(prim,n)*
     $                      (di(m+lx  ,1)*i2(m+ly  )*di(m+lz  ,1)+
     $                       di(m+lx+1,1)*i2(m+ly+1)*di(m+lz+1,1)+
     $                       di(m+lx+2,1)*i2(m+ly+2)*di(m+lz+2,1)+
     $                       di(m+lx+3,1)*i2(m+ly+3)*di(m+lz+3,1)+
     $                       di(m+lx+4,1)*i2(m+ly+4)*di(m+lz+4,1)+
     $                       di(m+lx+5,1)*i2(m+ly+5)*di(m+lz+5,1))
                        d2(5)=d2(5)+dens(prim,n)*
     $                      (i2(m+lx  )*di(m+ly  ,1)*di(m+lz  ,1)+
     $                       i2(m+lx+1)*di(m+ly+1,1)*di(m+lz+1,1)+
     $                       i2(m+lx+2)*di(m+ly+2,1)*di(m+lz+2,1)+
     $                       i2(m+lx+3)*di(m+ly+3,1)*di(m+lz+3,1)+
     $                       i2(m+lx+4)*di(m+ly+4,1)*di(m+lz+4,1)+
     $                       i2(m+lx+5)*di(m+ly+5,1)*di(m+lz+5,1))
                        d2(6)=d2(6)+dens(prim,n)*
     $                      (i2(m+lx  )*i2(m+ly  )*d2i(m+lz  ,1)+
     $                       i2(m+lx+1)*i2(m+ly+1)*d2i(m+lz+1,1)+
     $                       i2(m+lx+2)*i2(m+ly+2)*d2i(m+lz+2,1)+
     $                       i2(m+lx+3)*i2(m+ly+3)*d2i(m+lz+3,1)+
     $                       i2(m+lx+4)*i2(m+ly+4)*d2i(m+lz+4,1)+
     $                       i2(m+lx+5)*i2(m+ly+5)*d2i(m+lz+5,1))
                        m=m+6
 1000                continue
                     if (npass.lt.4) then
                        m=1
                        do 1001 prim=1,nv
                           der(1,2)=der(1,2)+dens(prim,n)*
     $                        (di(m+lx  ,2)*i2(m+ly  )*i2(m+lz  )+
     $                         di(m+lx+1,2)*i2(m+ly+1)*i2(m+lz+1)+
     $                         di(m+lx+2,2)*i2(m+ly+2)*i2(m+lz+2)+
     $                         di(m+lx+3,2)*i2(m+ly+3)*i2(m+lz+3)+
     $                         di(m+lx+4,2)*i2(m+ly+4)*i2(m+lz+4)+
     $                         di(m+lx+5,2)*i2(m+ly+5)*i2(m+lz+5))
                           der(2,2)=der(2,2)+dens(prim,n)*
     $                        (i2(m+lx  )*di(m+ly  ,2)*i2(m+lz  )+
     $                         i2(m+lx+1)*di(m+ly+1,2)*i2(m+lz+1)+
     $                         i2(m+lx+2)*di(m+ly+2,2)*i2(m+lz+2)+
     $                         i2(m+lx+3)*di(m+ly+3,2)*i2(m+lz+3)+
     $                         i2(m+lx+4)*di(m+ly+4,2)*i2(m+lz+4)+
     $                         i2(m+lx+5)*di(m+ly+5,2)*i2(m+lz+5))
                           der(3,2)=der(3,2)+dens(prim,n)*
     $                        (i2(m+lx  )*i2(m+ly  )*di(m+lz  ,2)+
     $                         i2(m+lx+1)*i2(m+ly+1)*di(m+lz+1,2)+
     $                         i2(m+lx+2)*i2(m+ly+2)*di(m+lz+2,2)+
     $                         i2(m+lx+3)*i2(m+ly+3)*di(m+lz+3,2)+
     $                         i2(m+lx+4)*i2(m+ly+4)*di(m+lz+4,2)+
     $                         i2(m+lx+5)*i2(m+ly+5)*di(m+lz+5,2))
                           d2(7)=d2(7)+dens(prim,n)*
     $                         (d2i(m+lx  ,2)*i2(m+ly  )*i2(m+lz  )+
     $                          d2i(m+lx+1,2)*i2(m+ly+1)*i2(m+lz+1)+
     $                          d2i(m+lx+2,2)*i2(m+ly+2)*i2(m+lz+2)+
     $                          d2i(m+lx+3,2)*i2(m+ly+3)*i2(m+lz+3)+
     $                          d2i(m+lx+4,2)*i2(m+ly+4)*i2(m+lz+4)+
     $                          d2i(m+lx+5,2)*i2(m+ly+5)*i2(m+lz+5))
                           d2(8)=d2(8)+dens(prim,n)*
     $                         (di(m+lx  ,2)*di(m+ly  ,1)*i2(m+lz  )+
     $                          di(m+lx+1,2)*di(m+ly+1,1)*i2(m+lz+1)+
     $                          di(m+lx+2,2)*di(m+ly+2,1)*i2(m+lz+2)+
     $                          di(m+lx+3,2)*di(m+ly+3,1)*i2(m+lz+3)+
     $                          di(m+lx+4,2)*di(m+ly+4,1)*i2(m+lz+4)+
     $                          di(m+lx+5,2)*di(m+ly+5,1)*i2(m+lz+5))
                           d2(9)=d2(9)+dens(prim,n)*
     $                         (di(m+lx  ,2)*i2(m+ly  )*di(m+lz  ,1)+
     $                          di(m+lx+1,2)*i2(m+ly+1)*di(m+lz+1,1)+
     $                          di(m+lx+2,2)*i2(m+ly+2)*di(m+lz+2,1)+
     $                          di(m+lx+3,2)*i2(m+ly+3)*di(m+lz+3,1)+
     $                          di(m+lx+4,2)*i2(m+ly+4)*di(m+lz+4,1)+
     $                          di(m+lx+5,2)*i2(m+ly+5)*di(m+lz+5,1))
                           d2(10)=d2(10)+dens(prim,n)*
     $                         (d2i(m+lx  ,3)*i2(m+ly  )*i2(m+lz  )+
     $                          d2i(m+lx+1,3)*i2(m+ly+1)*i2(m+lz+1)+
     $                          d2i(m+lx+2,3)*i2(m+ly+2)*i2(m+lz+2)+
     $                          d2i(m+lx+3,3)*i2(m+ly+3)*i2(m+lz+3)+
     $                          d2i(m+lx+4,3)*i2(m+ly+4)*i2(m+lz+4)+
     $                          d2i(m+lx+5,3)*i2(m+ly+5)*i2(m+lz+5))
                           d2(11)=d2(11)+dens(prim,n)*
     $                         (di(m+lx  ,1)*di(m+ly  ,2)*i2(m+lz  )+
     $                          di(m+lx+1,1)*di(m+ly+1,2)*i2(m+lz+1)+
     $                          di(m+lx+2,1)*di(m+ly+2,2)*i2(m+lz+2)+
     $                          di(m+lx+3,1)*di(m+ly+3,2)*i2(m+lz+3)+
     $                          di(m+lx+4,1)*di(m+ly+4,2)*i2(m+lz+4)+
     $                          di(m+lx+5,1)*di(m+ly+5,2)*i2(m+lz+5))
                           d2(12)=d2(12)+dens(prim,n)*
     $                         (i2(m+lx  )*d2i(m+ly  ,2)*i2(m+lz  )+
     $                          i2(m+lx+1)*d2i(m+ly+1,2)*i2(m+lz+1)+
     $                          i2(m+lx+2)*d2i(m+ly+2,2)*i2(m+lz+2)+
     $                          i2(m+lx+3)*d2i(m+ly+3,2)*i2(m+lz+3)+
     $                          i2(m+lx+4)*d2i(m+ly+4,2)*i2(m+lz+4)+
     $                          i2(m+lx+5)*d2i(m+ly+5,2)*i2(m+lz+5))
                           d2(13)=d2(13)+dens(prim,n)*
     $                         (i2(m+lx  )*di(m+ly  ,2)*di(m+lz  ,1)+
     $                          i2(m+lx+1)*di(m+ly+1,2)*di(m+lz+1,1)+
     $                          i2(m+lx+2)*di(m+ly+2,2)*di(m+lz+2,1)+
     $                          i2(m+lx+3)*di(m+ly+3,2)*di(m+lz+3,1)+
     $                          i2(m+lx+4)*di(m+ly+4,2)*di(m+lz+4,1)+
     $                          i2(m+lx+5)*di(m+ly+5,2)*di(m+lz+5,1))
                           d2(14)=d2(14)+dens(prim,n)*
     $                         (di(m+lx  ,2)*di(m+ly  ,2)*i2(m+lz  )+
     $                          di(m+lx+1,2)*di(m+ly+1,2)*i2(m+lz+1)+
     $                          di(m+lx+2,2)*di(m+ly+2,2)*i2(m+lz+2)+
     $                          di(m+lx+3,2)*di(m+ly+3,2)*i2(m+lz+3)+
     $                          di(m+lx+4,2)*di(m+ly+4,2)*i2(m+lz+4)+
     $                          di(m+lx+5,2)*di(m+ly+5,2)*i2(m+lz+5))
                           d2(15)=d2(15)+dens(prim,n)*
     $                         (i2(m+lx  )*d2i(m+ly  ,3)*i2(m+lz  )+
     $                          i2(m+lx+1)*d2i(m+ly+1,3)*i2(m+lz+1)+
     $                          i2(m+lx+2)*d2i(m+ly+2,3)*i2(m+lz+2)+
     $                          i2(m+lx+3)*d2i(m+ly+3,3)*i2(m+lz+3)+
     $                          i2(m+lx+4)*d2i(m+ly+4,3)*i2(m+lz+4)+
     $                          i2(m+lx+5)*d2i(m+ly+5,3)*i2(m+lz+5))
                           d2(16)=d2(16)+dens(prim,n)*
     $                         (di(m+lx  ,1)*i2(m+ly  )*di(m+lz  ,2)+
     $                          di(m+lx+1,1)*i2(m+ly+1)*di(m+lz+1,2)+
     $                          di(m+lx+2,1)*i2(m+ly+2)*di(m+lz+2,2)+
     $                          di(m+lx+3,1)*i2(m+ly+3)*di(m+lz+3,2)+
     $                          di(m+lx+4,1)*i2(m+ly+4)*di(m+lz+4,2)+
     $                          di(m+lx+5,1)*i2(m+ly+5)*di(m+lz+5,2))
                           d2(17)=d2(17)+dens(prim,n)*
     $                         (i2(m+lx  )*di(m+ly  ,1)*di(m+lz  ,2)+
     $                          i2(m+lx+1)*di(m+ly+1,1)*di(m+lz+1,2)+
     $                          i2(m+lx+2)*di(m+ly+2,1)*di(m+lz+2,2)+
     $                          i2(m+lx+3)*di(m+ly+3,1)*di(m+lz+3,2)+
     $                          i2(m+lx+4)*di(m+ly+4,1)*di(m+lz+4,2)+
     $                          i2(m+lx+5)*di(m+ly+5,1)*di(m+lz+5,2))
                           d2(18)=d2(18)+dens(prim,n)*
     $                         (i2(m+lx  )*i2(m+ly  )*d2i(m+lz  ,2)+
     $                          i2(m+lx+1)*i2(m+ly+1)*d2i(m+lz+1,2)+
     $                          i2(m+lx+2)*i2(m+ly+2)*d2i(m+lz+2,2)+
     $                          i2(m+lx+3)*i2(m+ly+3)*d2i(m+lz+3,2)+
     $                          i2(m+lx+4)*i2(m+ly+4)*d2i(m+lz+4,2)+
     $                          i2(m+lx+5)*i2(m+ly+5)*d2i(m+lz+5,2))
                           d2(19)=d2(19)+dens(prim,n)*
     $                         (di(m+lx  ,2)*i2(m+ly  )*di(m+lz  ,2)+
     $                          di(m+lx+1,2)*i2(m+ly+1)*di(m+lz+1,2)+
     $                          di(m+lx+2,2)*i2(m+ly+2)*di(m+lz+2,2)+
     $                          di(m+lx+3,2)*i2(m+ly+3)*di(m+lz+3,2)+
     $                          di(m+lx+4,2)*i2(m+ly+4)*di(m+lz+4,2)+
     $                          di(m+lx+5,2)*i2(m+ly+5)*di(m+lz+5,2))
                           d2(20)=d2(20)+dens(prim,n)*
     $                         (i2(m+lx  )*di(m+ly  ,2)*di(m+lz  ,2)+
     $                          i2(m+lx+1)*di(m+ly+1,2)*di(m+lz+1,2)+
     $                          i2(m+lx+2)*di(m+ly+2,2)*di(m+lz+2,2)+
     $                          i2(m+lx+3)*di(m+ly+3,2)*di(m+lz+3,2)+
     $                          i2(m+lx+4)*di(m+ly+4,2)*di(m+lz+4,2)+
     $                          i2(m+lx+5)*di(m+ly+5,2)*di(m+lz+5,2))
                           d2(21)=d2(21)+dens(prim,n)*
     $                         (i2(m+lx  )*i2(m+ly  )*d2i(m+lz  ,3)+
     $                          i2(m+lx+1)*i2(m+ly+1)*d2i(m+lz+1,3)+
     $                          i2(m+lx+2)*i2(m+ly+2)*d2i(m+lz+2,3)+
     $                          i2(m+lx+3)*i2(m+ly+3)*d2i(m+lz+3,3)+
     $                          i2(m+lx+4)*i2(m+ly+4)*d2i(m+lz+4,3)+
     $                          i2(m+lx+5)*i2(m+ly+5)*d2i(m+lz+5,3))
                           m=m+6
 1001                   continue
                     end if
                     if (npass.eq.1) then
                        m=1
                        do 1002 prim=1,nv
                           der(1,3)=der(1,3)+dens(prim,n)*
     $                        (di(m+lx  ,3)*i2(m+ly  )*i2(m+lz  )+
     $                         di(m+lx+1,3)*i2(m+ly+1)*i2(m+lz+1)+
     $                         di(m+lx+2,3)*i2(m+ly+2)*i2(m+lz+2)+
     $                         di(m+lx+3,3)*i2(m+ly+3)*i2(m+lz+3)+
     $                         di(m+lx+4,3)*i2(m+ly+4)*i2(m+lz+4)+
     $                         di(m+lx+5,3)*i2(m+ly+5)*i2(m+lz+5))
                           der(2,3)=der(2,3)+dens(prim,n)*
     $                        (i2(m+lx  )*di(m+ly  ,3)*i2(m+lz  )+
     $                         i2(m+lx+1)*di(m+ly+1,3)*i2(m+lz+1)+
     $                         i2(m+lx+2)*di(m+ly+2,3)*i2(m+lz+2)+
     $                         i2(m+lx+3)*di(m+ly+3,3)*i2(m+lz+3)+
     $                         i2(m+lx+4)*di(m+ly+4,3)*i2(m+lz+4)+
     $                         i2(m+lx+5)*di(m+ly+5,3)*i2(m+lz+5))
                           der(3,3)=der(3,3)+dens(prim,n)*
     $                        (i2(m+lx  )*i2(m+ly  )*di(m+lz  ,3)+
     $                         i2(m+lx+1)*i2(m+ly+1)*di(m+lz+1,3)+
     $                         i2(m+lx+2)*i2(m+ly+2)*di(m+lz+2,3)+
     $                         i2(m+lx+3)*i2(m+ly+3)*di(m+lz+3,3)+
     $                         i2(m+lx+4)*i2(m+ly+4)*di(m+lz+4,3)+
     $                         i2(m+lx+5)*i2(m+ly+5)*di(m+lz+5,3))
                           d2(22)=d2(22)+dens(prim,n)*
     $                         (d2i(m+lx  ,4)*i2(m+ly  )*i2(m+lz  )+
     $                          d2i(m+lx+1,4)*i2(m+ly+1)*i2(m+lz+1)+
     $                          d2i(m+lx+2,4)*i2(m+ly+2)*i2(m+lz+2)+
     $                          d2i(m+lx+3,4)*i2(m+ly+3)*i2(m+lz+3)+
     $                          d2i(m+lx+4,4)*i2(m+ly+4)*i2(m+lz+4)+
     $                          d2i(m+lx+5,4)*i2(m+ly+5)*i2(m+lz+5))
                           d2(23)=d2(23)+dens(prim,n)*
     $                         (di(m+lx  ,3)*di(m+ly  ,1)*i2(m+lz  )+
     $                          di(m+lx+1,3)*di(m+ly+1,1)*i2(m+lz+1)+
     $                          di(m+lx+2,3)*di(m+ly+2,1)*i2(m+lz+2)+
     $                          di(m+lx+3,3)*di(m+ly+3,1)*i2(m+lz+3)+
     $                          di(m+lx+4,3)*di(m+ly+4,1)*i2(m+lz+4)+
     $                          di(m+lx+5,3)*di(m+ly+5,1)*i2(m+lz+5))
                           d2(24)=d2(24)+dens(prim,n)*
     $                         (di(m+lx  ,3)*i2(m+ly  )*di(m+lz  ,1)+
     $                          di(m+lx+1,3)*i2(m+ly+1)*di(m+lz+1,1)+
     $                          di(m+lx+2,3)*i2(m+ly+2)*di(m+lz+2,1)+
     $                          di(m+lx+3,3)*i2(m+ly+3)*di(m+lz+3,1)+
     $                          di(m+lx+4,3)*i2(m+ly+4)*di(m+lz+4,1)+
     $                          di(m+lx+5,3)*i2(m+ly+5)*di(m+lz+5,1))
                           d2(25)=d2(25)+dens(prim,n)*
     $                         (d2i(m+lx  ,5)*i2(m+ly  )*i2(m+lz  )+
     $                          d2i(m+lx+1,5)*i2(m+ly+1)*i2(m+lz+1)+
     $                          d2i(m+lx+2,5)*i2(m+ly+2)*i2(m+lz+2)+
     $                          d2i(m+lx+3,5)*i2(m+ly+3)*i2(m+lz+3)+
     $                          d2i(m+lx+4,5)*i2(m+ly+4)*i2(m+lz+4)+
     $                          d2i(m+lx+5,5)*i2(m+ly+5)*i2(m+lz+5))
                           d2(26)=d2(26)+dens(prim,n)*
     $                         (di(m+lx  ,3)*di(m+ly  ,2)*i2(m+lz  )+
     $                          di(m+lx+1,3)*di(m+ly+1,2)*i2(m+lz+1)+
     $                          di(m+lx+2,3)*di(m+ly+2,2)*i2(m+lz+2)+
     $                          di(m+lx+3,3)*di(m+ly+3,2)*i2(m+lz+3)+
     $                          di(m+lx+4,3)*di(m+ly+4,2)*i2(m+lz+4)+
     $                          di(m+lx+5,3)*di(m+ly+5,2)*i2(m+lz+5))
                           d2(27)=d2(27)+dens(prim,n)*
     $                         (di(m+lx  ,3)*i2(m+ly  )*di(m+lz  ,2)+
     $                          di(m+lx+1,3)*i2(m+ly+1)*di(m+lz+1,2)+
     $                          di(m+lx+2,3)*i2(m+ly+2)*di(m+lz+2,2)+
     $                          di(m+lx+3,3)*i2(m+ly+3)*di(m+lz+3,2)+
     $                          di(m+lx+4,3)*i2(m+ly+4)*di(m+lz+4,2)+
     $                          di(m+lx+5,3)*i2(m+ly+5)*di(m+lz+5,2))
                           d2(28)=d2(28)+dens(prim,n)*
     $                         (d2i(m+lx  ,6)*i2(m+ly  )*i2(m+lz  )+
     $                          d2i(m+lx+1,6)*i2(m+ly+1)*i2(m+lz+1)+
     $                          d2i(m+lx+2,6)*i2(m+ly+2)*i2(m+lz+2)+
     $                          d2i(m+lx+3,6)*i2(m+ly+3)*i2(m+lz+3)+
     $                          d2i(m+lx+4,6)*i2(m+ly+4)*i2(m+lz+4)+
     $                          d2i(m+lx+5,6)*i2(m+ly+5)*i2(m+lz+5))
                           d2(29)=d2(29)+dens(prim,n)*
     $                         (di(m+lx  ,1)*di(m+ly  ,3)*i2(m+lz  )+
     $                          di(m+lx+1,1)*di(m+ly+1,3)*i2(m+lz+1)+
     $                          di(m+lx+2,1)*di(m+ly+2,3)*i2(m+lz+2)+
     $                          di(m+lx+3,1)*di(m+ly+3,3)*i2(m+lz+3)+
     $                          di(m+lx+4,1)*di(m+ly+4,3)*i2(m+lz+4)+
     $                          di(m+lx+5,1)*di(m+ly+5,3)*i2(m+lz+5))
                           d2(30)=d2(30)+dens(prim,n)*
     $                         (i2(m+lx  )*d2i(m+ly  ,4)*i2(m+lz  )+
     $                          i2(m+lx+1)*d2i(m+ly+1,4)*i2(m+lz+1)+
     $                          i2(m+lx+2)*d2i(m+ly+2,4)*i2(m+lz+2)+
     $                          i2(m+lx+3)*d2i(m+ly+3,4)*i2(m+lz+3)+
     $                          i2(m+lx+4)*d2i(m+ly+4,4)*i2(m+lz+4)+
     $                          i2(m+lx+5)*d2i(m+ly+5,4)*i2(m+lz+5))
                           d2(31)=d2(31)+dens(prim,n)*
     $                         (i2(m+lx  )*di(m+ly  ,3)*di(m+lz  ,1)+
     $                          i2(m+lx+1)*di(m+ly+1,3)*di(m+lz+1,1)+
     $                          i2(m+lx+2)*di(m+ly+2,3)*di(m+lz+2,1)+
     $                          i2(m+lx+3)*di(m+ly+3,3)*di(m+lz+3,1)+
     $                          i2(m+lx+4)*di(m+ly+4,3)*di(m+lz+4,1)+
     $                          i2(m+lx+5)*di(m+ly+5,3)*di(m+lz+5,1))
                           d2(32)=d2(32)+dens(prim,n)*
     $                         (di(m+lx  ,2)*di(m+ly  ,3)*i2(m+lz  )+
     $                          di(m+lx+1,2)*di(m+ly+1,3)*i2(m+lz+1)+
     $                          di(m+lx+2,2)*di(m+ly+2,3)*i2(m+lz+2)+
     $                          di(m+lx+3,2)*di(m+ly+3,3)*i2(m+lz+3)+
     $                          di(m+lx+4,2)*di(m+ly+4,3)*i2(m+lz+4)+
     $                          di(m+lx+5,2)*di(m+ly+5,3)*i2(m+lz+5))
                           d2(33)=d2(33)+dens(prim,n)*
     $                         (i2(m+lx  )*d2i(m+ly  ,5)*i2(m+lz  )+
     $                          i2(m+lx+1)*d2i(m+ly+1,5)*i2(m+lz+1)+
     $                          i2(m+lx+2)*d2i(m+ly+2,5)*i2(m+lz+2)+
     $                          i2(m+lx+3)*d2i(m+ly+3,5)*i2(m+lz+3)+
     $                          i2(m+lx+4)*d2i(m+ly+4,5)*i2(m+lz+4)+
     $                          i2(m+lx+5)*d2i(m+ly+5,5)*i2(m+lz+5))
                           d2(34)=d2(34)+dens(prim,n)*
     $                         (i2(m+lx  )*di(m+ly  ,3)*di(m+lz  ,2)+
     $                          i2(m+lx+1)*di(m+ly+1,3)*di(m+lz+1,2)+
     $                          i2(m+lx+2)*di(m+ly+2,3)*di(m+lz+2,2)+
     $                          i2(m+lx+3)*di(m+ly+3,3)*di(m+lz+3,2)+
     $                          i2(m+lx+4)*di(m+ly+4,3)*di(m+lz+4,2)+
     $                          i2(m+lx+5)*di(m+ly+5,3)*di(m+lz+5,2))
                           d2(35)=d2(35)+dens(prim,n)*
     $                         (di(m+lx  ,3)*di(m+ly  ,3)*i2(m+lz  )+
     $                          di(m+lx+1,3)*di(m+ly+1,3)*i2(m+lz+1)+
     $                          di(m+lx+2,3)*di(m+ly+2,3)*i2(m+lz+2)+
     $                          di(m+lx+3,3)*di(m+ly+3,3)*i2(m+lz+3)+
     $                          di(m+lx+4,3)*di(m+ly+4,3)*i2(m+lz+4)+
     $                          di(m+lx+5,3)*di(m+ly+5,3)*i2(m+lz+5))
                           d2(36)=d2(36)+dens(prim,n)*
     $                         (i2(m+lx  )*d2i(m+ly  ,6)*i2(m+lz  )+
     $                          i2(m+lx+1)*d2i(m+ly+1,6)*i2(m+lz+1)+
     $                          i2(m+lx+2)*d2i(m+ly+2,6)*i2(m+lz+2)+
     $                          i2(m+lx+3)*d2i(m+ly+3,6)*i2(m+lz+3)+
     $                          i2(m+lx+4)*d2i(m+ly+4,6)*i2(m+lz+4)+
     $                          i2(m+lx+5)*d2i(m+ly+5,6)*i2(m+lz+5))
                           d2(37)=d2(37)+dens(prim,n)*
     $                         (di(m+lx  ,1)*i2(m+ly  )*di(m+lz  ,3)+
     $                          di(m+lx+1,1)*i2(m+ly+1)*di(m+lz+1,3)+
     $                          di(m+lx+2,1)*i2(m+ly+2)*di(m+lz+2,3)+
     $                          di(m+lx+3,1)*i2(m+ly+3)*di(m+lz+3,3)+
     $                          di(m+lx+4,1)*i2(m+ly+4)*di(m+lz+4,3)+
     $                          di(m+lx+5,1)*i2(m+ly+5)*di(m+lz+5,3))
                           d2(38)=d2(38)+dens(prim,n)*
     $                         (i2(m+lx  )*di(m+ly  ,1)*di(m+lz  ,3)+
     $                          i2(m+lx+1)*di(m+ly+1,1)*di(m+lz+1,3)+
     $                          i2(m+lx+2)*di(m+ly+2,1)*di(m+lz+2,3)+
     $                          i2(m+lx+3)*di(m+ly+3,1)*di(m+lz+3,3)+
     $                          i2(m+lx+4)*di(m+ly+4,1)*di(m+lz+4,3)+
     $                          i2(m+lx+5)*di(m+ly+5,1)*di(m+lz+5,3))
                           d2(39)=d2(39)+dens(prim,n)*
     $                         (i2(m+lx  )*i2(m+ly  )*d2i(m+lz  ,4)+
     $                          i2(m+lx+1)*i2(m+ly+1)*d2i(m+lz+1,4)+
     $                          i2(m+lx+2)*i2(m+ly+2)*d2i(m+lz+2,4)+
     $                          i2(m+lx+3)*i2(m+ly+3)*d2i(m+lz+3,4)+
     $                          i2(m+lx+4)*i2(m+ly+4)*d2i(m+lz+4,4)+
     $                          i2(m+lx+5)*i2(m+ly+5)*d2i(m+lz+5,4))
                           d2(40)=d2(40)+dens(prim,n)*
     $                         (di(m+lx  ,2)*i2(m+ly  )*di(m+lz  ,3)+
     $                          di(m+lx+1,2)*i2(m+ly+1)*di(m+lz+1,3)+
     $                          di(m+lx+2,2)*i2(m+ly+2)*di(m+lz+2,3)+
     $                          di(m+lx+3,2)*i2(m+ly+3)*di(m+lz+3,3)+
     $                          di(m+lx+4,2)*i2(m+ly+4)*di(m+lz+4,3)+
     $                          di(m+lx+5,2)*i2(m+ly+5)*di(m+lz+5,3))
                           d2(41)=d2(41)+dens(prim,n)*
     $                         (i2(m+lx  )*di(m+ly  ,2)*di(m+lz  ,3)+
     $                          i2(m+lx+1)*di(m+ly+1,2)*di(m+lz+1,3)+
     $                          i2(m+lx+2)*di(m+ly+2,2)*di(m+lz+2,3)+
     $                          i2(m+lx+3)*di(m+ly+3,2)*di(m+lz+3,3)+
     $                          i2(m+lx+4)*di(m+ly+4,2)*di(m+lz+4,3)+
     $                          i2(m+lx+5)*di(m+ly+5,2)*di(m+lz+5,3))
                           d2(42)=d2(42)+dens(prim,n)*
     $                         (i2(m+lx  )*i2(m+ly  )*d2i(m+lz  ,5)+
     $                          i2(m+lx+1)*i2(m+ly+1)*d2i(m+lz+1,5)+
     $                          i2(m+lx+2)*i2(m+ly+2)*d2i(m+lz+2,5)+
     $                          i2(m+lx+3)*i2(m+ly+3)*d2i(m+lz+3,5)+
     $                          i2(m+lx+4)*i2(m+ly+4)*d2i(m+lz+4,5)+
     $                          i2(m+lx+5)*i2(m+ly+5)*d2i(m+lz+5,5))
                           d2(43)=d2(43)+dens(prim,n)*
     $                         (di(m+lx  ,3)*i2(m+ly  )*di(m+lz  ,3)+
     $                          di(m+lx+1,3)*i2(m+ly+1)*di(m+lz+1,3)+
     $                          di(m+lx+2,3)*i2(m+ly+2)*di(m+lz+2,3)+
     $                          di(m+lx+3,3)*i2(m+ly+3)*di(m+lz+3,3)+
     $                          di(m+lx+4,3)*i2(m+ly+4)*di(m+lz+4,3)+
     $                          di(m+lx+5,3)*i2(m+ly+5)*di(m+lz+5,3))
                           d2(44)=d2(44)+dens(prim,n)*
     $                         (i2(m+lx  )*di(m+ly  ,3)*di(m+lz  ,3)+
     $                          i2(m+lx+1)*di(m+ly+1,3)*di(m+lz+1,3)+
     $                          i2(m+lx+2)*di(m+ly+2,3)*di(m+lz+2,3)+
     $                          i2(m+lx+3)*di(m+ly+3,3)*di(m+lz+3,3)+
     $                          i2(m+lx+4)*di(m+ly+4,3)*di(m+lz+4,3)+
     $                          i2(m+lx+5)*di(m+ly+5,3)*di(m+lz+5,3))
                           d2(45)=d2(45)+dens(prim,n)*
     $                         (i2(m+lx  )*i2(m+ly  )*d2i(m+lz  ,6)+
     $                          i2(m+lx+1)*i2(m+ly+1)*d2i(m+lz+1,6)+
     $                          i2(m+lx+2)*i2(m+ly+2)*d2i(m+lz+2,6)+
     $                          i2(m+lx+3)*i2(m+ly+3)*d2i(m+lz+3,6)+
     $                          i2(m+lx+4)*i2(m+ly+4)*d2i(m+lz+4,6)+
     $                          i2(m+lx+5)*i2(m+ly+5)*d2i(m+lz+5,6))
                           m=m+6
 1002                   continue
                     end if
                     go to 15
c
c
   11             continue
                     m=1
                     do 1100 prim=1,nv
                        energy=energy+dens(prim,n)*
     $                      (i2(m+lx  )*i2(m+ly  )*i2(m+lz  )+
     $                       i2(m+lx+1)*i2(m+ly+1)*i2(m+lz+1)+
     $                       i2(m+lx+2)*i2(m+ly+2)*i2(m+lz+2)+
     $                       i2(m+lx+3)*i2(m+ly+3)*i2(m+lz+3)+
     $                       i2(m+lx+4)*i2(m+ly+4)*i2(m+lz+4)+
     $                       i2(m+lx+5)*i2(m+ly+5)*i2(m+lz+5)+
     $                       i2(m+lx+6)*i2(m+ly+6)*i2(m+lz+6))
                        der(1,1)=der(1,1)+dens(prim,n)*
     $                      (di(m+lx  ,1)*i2(m+ly  )*i2(m+lz  )+
     $                       di(m+lx+1,1)*i2(m+ly+1)*i2(m+lz+1)+
     $                       di(m+lx+2,1)*i2(m+ly+2)*i2(m+lz+2)+
     $                       di(m+lx+3,1)*i2(m+ly+3)*i2(m+lz+3)+
     $                       di(m+lx+4,1)*i2(m+ly+4)*i2(m+lz+4)+
     $                       di(m+lx+5,1)*i2(m+ly+5)*i2(m+lz+5)+
     $                       di(m+lx+6,1)*i2(m+ly+6)*i2(m+lz+6))
                        der(2,1)=der(2,1)+dens(prim,n)*
     $                      (i2(m+lx  )*di(m+ly  ,1)*i2(m+lz  )+
     $                       i2(m+lx+1)*di(m+ly+1,1)*i2(m+lz+1)+
     $                       i2(m+lx+2)*di(m+ly+2,1)*i2(m+lz+2)+
     $                       i2(m+lx+3)*di(m+ly+3,1)*i2(m+lz+3)+
     $                       i2(m+lx+4)*di(m+ly+4,1)*i2(m+lz+4)+
     $                       i2(m+lx+5)*di(m+ly+5,1)*i2(m+lz+5)+
     $                       i2(m+lx+6)*di(m+ly+6,1)*i2(m+lz+6))
                        der(3,1)=der(3,1)+dens(prim,n)*
     $                      (i2(m+lx  )*i2(m+ly  )*di(m+lz  ,1)+
     $                       i2(m+lx+1)*i2(m+ly+1)*di(m+lz+1,1)+
     $                       i2(m+lx+2)*i2(m+ly+2)*di(m+lz+2,1)+
     $                       i2(m+lx+3)*i2(m+ly+3)*di(m+lz+3,1)+
     $                       i2(m+lx+4)*i2(m+ly+4)*di(m+lz+4,1)+
     $                       i2(m+lx+5)*i2(m+ly+5)*di(m+lz+5,1)+
     $                       i2(m+lx+6)*i2(m+ly+6)*di(m+lz+6,1))
                        d2(1)=d2(1)+dens(prim,n)*
     $                      (d2i(m+lx  ,1)*i2(m+ly  )*i2(m+lz  )+
     $                       d2i(m+lx+1,1)*i2(m+ly+1)*i2(m+lz+1)+
     $                       d2i(m+lx+2,1)*i2(m+ly+2)*i2(m+lz+2)+
     $                       d2i(m+lx+3,1)*i2(m+ly+3)*i2(m+lz+3)+
     $                       d2i(m+lx+4,1)*i2(m+ly+4)*i2(m+lz+4)+
     $                       d2i(m+lx+5,1)*i2(m+ly+5)*i2(m+lz+5)+
     $                       d2i(m+lx+6,1)*i2(m+ly+6)*i2(m+lz+6))
                        d2(2)=d2(2)+dens(prim,n)*
     $                      (di(m+lx  ,1)*di(m+ly  ,1)*i2(m+lz  )+
     $                       di(m+lx+1,1)*di(m+ly+1,1)*i2(m+lz+1)+
     $                       di(m+lx+2,1)*di(m+ly+2,1)*i2(m+lz+2)+
     $                       di(m+lx+3,1)*di(m+ly+3,1)*i2(m+lz+3)+
     $                       di(m+lx+4,1)*di(m+ly+4,1)*i2(m+lz+4)+
     $                       di(m+lx+5,1)*di(m+ly+5,1)*i2(m+lz+5)+
     $                       di(m+lx+6,1)*di(m+ly+6,1)*i2(m+lz+6))
                        d2(3)=d2(3)+dens(prim,n)*
     $                      (i2(m+lx  )*d2i(m+ly  ,1)*i2(m+lz  )+
     $                       i2(m+lx+1)*d2i(m+ly+1,1)*i2(m+lz+1)+
     $                       i2(m+lx+2)*d2i(m+ly+2,1)*i2(m+lz+2)+
     $                       i2(m+lx+3)*d2i(m+ly+3,1)*i2(m+lz+3)+
     $                       i2(m+lx+4)*d2i(m+ly+4,1)*i2(m+lz+4)+
     $                       i2(m+lx+5)*d2i(m+ly+5,1)*i2(m+lz+5)+
     $                       i2(m+lx+6)*d2i(m+ly+6,1)*i2(m+lz+6))
                        d2(4)=d2(4)+dens(prim,n)*
     $                      (di(m+lx  ,1)*i2(m+ly  )*di(m+lz  ,1)+
     $                       di(m+lx+1,1)*i2(m+ly+1)*di(m+lz+1,1)+
     $                       di(m+lx+2,1)*i2(m+ly+2)*di(m+lz+2,1)+
     $                       di(m+lx+3,1)*i2(m+ly+3)*di(m+lz+3,1)+
     $                       di(m+lx+4,1)*i2(m+ly+4)*di(m+lz+4,1)+
     $                       di(m+lx+5,1)*i2(m+ly+5)*di(m+lz+5,1)+
     $                       di(m+lx+6,1)*i2(m+ly+6)*di(m+lz+6,1))
                        d2(5)=d2(5)+dens(prim,n)*
     $                      (i2(m+lx  )*di(m+ly  ,1)*di(m+lz  ,1)+
     $                       i2(m+lx+1)*di(m+ly+1,1)*di(m+lz+1,1)+
     $                       i2(m+lx+2)*di(m+ly+2,1)*di(m+lz+2,1)+
     $                       i2(m+lx+3)*di(m+ly+3,1)*di(m+lz+3,1)+
     $                       i2(m+lx+4)*di(m+ly+4,1)*di(m+lz+4,1)+
     $                       i2(m+lx+5)*di(m+ly+5,1)*di(m+lz+5,1)+
     $                       i2(m+lx+6)*di(m+ly+6,1)*di(m+lz+6,1))
                        d2(6)=d2(6)+dens(prim,n)*
     $                      (i2(m+lx  )*i2(m+ly  )*d2i(m+lz  ,1)+
     $                       i2(m+lx+1)*i2(m+ly+1)*d2i(m+lz+1,1)+
     $                       i2(m+lx+2)*i2(m+ly+2)*d2i(m+lz+2,1)+
     $                       i2(m+lx+3)*i2(m+ly+3)*d2i(m+lz+3,1)+
     $                       i2(m+lx+4)*i2(m+ly+4)*d2i(m+lz+4,1)+
     $                       i2(m+lx+5)*i2(m+ly+5)*d2i(m+lz+5,1)+
     $                       i2(m+lx+6)*i2(m+ly+6)*d2i(m+lz+6,1))
                        m=m+7
 1100                continue
                     if (npass.lt.4) then
                        m=1
                        do 1101 prim=1,nv
                           der(1,2)=der(1,2)+dens(prim,n)*
     $                        (di(m+lx  ,2)*i2(m+ly  )*i2(m+lz  )+
     $                         di(m+lx+1,2)*i2(m+ly+1)*i2(m+lz+1)+
     $                         di(m+lx+2,2)*i2(m+ly+2)*i2(m+lz+2)+
     $                         di(m+lx+3,2)*i2(m+ly+3)*i2(m+lz+3)+
     $                         di(m+lx+4,2)*i2(m+ly+4)*i2(m+lz+4)+
     $                         di(m+lx+5,2)*i2(m+ly+5)*i2(m+lz+5)+
     $                         di(m+lx+6,2)*i2(m+ly+6)*i2(m+lz+6))
                           der(2,2)=der(2,2)+dens(prim,n)*
     $                        (i2(m+lx  )*di(m+ly  ,2)*i2(m+lz  )+
     $                         i2(m+lx+1)*di(m+ly+1,2)*i2(m+lz+1)+
     $                         i2(m+lx+2)*di(m+ly+2,2)*i2(m+lz+2)+
     $                         i2(m+lx+3)*di(m+ly+3,2)*i2(m+lz+3)+
     $                         i2(m+lx+4)*di(m+ly+4,2)*i2(m+lz+4)+
     $                         i2(m+lx+5)*di(m+ly+5,2)*i2(m+lz+5)+
     $                         i2(m+lx+6)*di(m+ly+6,2)*i2(m+lz+6))
                           der(3,2)=der(3,2)+dens(prim,n)*
     $                        (i2(m+lx  )*i2(m+ly  )*di(m+lz  ,2)+
     $                         i2(m+lx+1)*i2(m+ly+1)*di(m+lz+1,2)+
     $                         i2(m+lx+2)*i2(m+ly+2)*di(m+lz+2,2)+
     $                         i2(m+lx+3)*i2(m+ly+3)*di(m+lz+3,2)+
     $                         i2(m+lx+4)*i2(m+ly+4)*di(m+lz+4,2)+
     $                         i2(m+lx+5)*i2(m+ly+5)*di(m+lz+5,2)+
     $                         i2(m+lx+6)*i2(m+ly+6)*di(m+lz+6,2))
                           d2(7)=d2(7)+dens(prim,n)*
     $                         (d2i(m+lx  ,2)*i2(m+ly  )*i2(m+lz  )+
     $                          d2i(m+lx+1,2)*i2(m+ly+1)*i2(m+lz+1)+
     $                          d2i(m+lx+2,2)*i2(m+ly+2)*i2(m+lz+2)+
     $                          d2i(m+lx+3,2)*i2(m+ly+3)*i2(m+lz+3)+
     $                          d2i(m+lx+4,2)*i2(m+ly+4)*i2(m+lz+4)+
     $                          d2i(m+lx+5,2)*i2(m+ly+5)*i2(m+lz+5)+
     $                          d2i(m+lx+6,2)*i2(m+ly+6)*i2(m+lz+6))
                           d2(8)=d2(8)+dens(prim,n)*
     $                         (di(m+lx  ,2)*di(m+ly  ,1)*i2(m+lz  )+
     $                          di(m+lx+1,2)*di(m+ly+1,1)*i2(m+lz+1)+
     $                          di(m+lx+2,2)*di(m+ly+2,1)*i2(m+lz+2)+
     $                          di(m+lx+3,2)*di(m+ly+3,1)*i2(m+lz+3)+
     $                          di(m+lx+4,2)*di(m+ly+4,1)*i2(m+lz+4)+
     $                          di(m+lx+5,2)*di(m+ly+5,1)*i2(m+lz+5)+
     $                          di(m+lx+6,2)*di(m+ly+6,1)*i2(m+lz+6))
                           d2(9)=d2(9)+dens(prim,n)*
     $                         (di(m+lx  ,2)*i2(m+ly  )*di(m+lz  ,1)+
     $                          di(m+lx+1,2)*i2(m+ly+1)*di(m+lz+1,1)+
     $                          di(m+lx+2,2)*i2(m+ly+2)*di(m+lz+2,1)+
     $                          di(m+lx+3,2)*i2(m+ly+3)*di(m+lz+3,1)+
     $                          di(m+lx+4,2)*i2(m+ly+4)*di(m+lz+4,1)+
     $                          di(m+lx+5,2)*i2(m+ly+5)*di(m+lz+5,1)+
     $                          di(m+lx+6,2)*i2(m+ly+6)*di(m+lz+6,1))
                           d2(10)=d2(10)+dens(prim,n)*
     $                         (d2i(m+lx  ,3)*i2(m+ly  )*i2(m+lz  )+
     $                          d2i(m+lx+1,3)*i2(m+ly+1)*i2(m+lz+1)+
     $                          d2i(m+lx+2,3)*i2(m+ly+2)*i2(m+lz+2)+
     $                          d2i(m+lx+3,3)*i2(m+ly+3)*i2(m+lz+3)+
     $                          d2i(m+lx+4,3)*i2(m+ly+4)*i2(m+lz+4)+
     $                          d2i(m+lx+5,3)*i2(m+ly+5)*i2(m+lz+5)+
     $                          d2i(m+lx+6,3)*i2(m+ly+6)*i2(m+lz+6))
                           d2(11)=d2(11)+dens(prim,n)*
     $                         (di(m+lx  ,1)*di(m+ly  ,2)*i2(m+lz  )+
     $                          di(m+lx+1,1)*di(m+ly+1,2)*i2(m+lz+1)+
     $                          di(m+lx+2,1)*di(m+ly+2,2)*i2(m+lz+2)+
     $                          di(m+lx+3,1)*di(m+ly+3,2)*i2(m+lz+3)+
     $                          di(m+lx+4,1)*di(m+ly+4,2)*i2(m+lz+4)+
     $                          di(m+lx+5,1)*di(m+ly+5,2)*i2(m+lz+5)+
     $                          di(m+lx+6,1)*di(m+ly+6,2)*i2(m+lz+6))
                           d2(12)=d2(12)+dens(prim,n)*
     $                         (i2(m+lx  )*d2i(m+ly  ,2)*i2(m+lz  )+
     $                          i2(m+lx+1)*d2i(m+ly+1,2)*i2(m+lz+1)+
     $                          i2(m+lx+2)*d2i(m+ly+2,2)*i2(m+lz+2)+
     $                          i2(m+lx+3)*d2i(m+ly+3,2)*i2(m+lz+3)+
     $                          i2(m+lx+4)*d2i(m+ly+4,2)*i2(m+lz+4)+
     $                          i2(m+lx+5)*d2i(m+ly+5,2)*i2(m+lz+5)+
     $                          i2(m+lx+6)*d2i(m+ly+6,2)*i2(m+lz+6))
                           d2(13)=d2(13)+dens(prim,n)*
     $                         (i2(m+lx  )*di(m+ly  ,2)*di(m+lz  ,1)+
     $                          i2(m+lx+1)*di(m+ly+1,2)*di(m+lz+1,1)+
     $                          i2(m+lx+2)*di(m+ly+2,2)*di(m+lz+2,1)+
     $                          i2(m+lx+3)*di(m+ly+3,2)*di(m+lz+3,1)+
     $                          i2(m+lx+4)*di(m+ly+4,2)*di(m+lz+4,1)+
     $                          i2(m+lx+5)*di(m+ly+5,2)*di(m+lz+5,1)+
     $                          i2(m+lx+6)*di(m+ly+6,2)*di(m+lz+6,1))
                           d2(14)=d2(14)+dens(prim,n)*
     $                         (di(m+lx  ,2)*di(m+ly  ,2)*i2(m+lz  )+
     $                          di(m+lx+1,2)*di(m+ly+1,2)*i2(m+lz+1)+
     $                          di(m+lx+2,2)*di(m+ly+2,2)*i2(m+lz+2)+
     $                          di(m+lx+3,2)*di(m+ly+3,2)*i2(m+lz+3)+
     $                          di(m+lx+4,2)*di(m+ly+4,2)*i2(m+lz+4)+
     $                          di(m+lx+5,2)*di(m+ly+5,2)*i2(m+lz+5)+
     $                          di(m+lx+6,2)*di(m+ly+6,2)*i2(m+lz+6))
                           d2(15)=d2(15)+dens(prim,n)*
     $                         (i2(m+lx  )*d2i(m+ly  ,3)*i2(m+lz  )+
     $                          i2(m+lx+1)*d2i(m+ly+1,3)*i2(m+lz+1)+
     $                          i2(m+lx+2)*d2i(m+ly+2,3)*i2(m+lz+2)+
     $                          i2(m+lx+3)*d2i(m+ly+3,3)*i2(m+lz+3)+
     $                          i2(m+lx+4)*d2i(m+ly+4,3)*i2(m+lz+4)+
     $                          i2(m+lx+5)*d2i(m+ly+5,3)*i2(m+lz+5)+
     $                          i2(m+lx+6)*d2i(m+ly+6,3)*i2(m+lz+6))
                           d2(16)=d2(16)+dens(prim,n)*
     $                         (di(m+lx  ,1)*i2(m+ly  )*di(m+lz  ,2)+
     $                          di(m+lx+1,1)*i2(m+ly+1)*di(m+lz+1,2)+
     $                          di(m+lx+2,1)*i2(m+ly+2)*di(m+lz+2,2)+
     $                          di(m+lx+3,1)*i2(m+ly+3)*di(m+lz+3,2)+
     $                          di(m+lx+4,1)*i2(m+ly+4)*di(m+lz+4,2)+
     $                          di(m+lx+5,1)*i2(m+ly+5)*di(m+lz+5,2)+
     $                          di(m+lx+6,1)*i2(m+ly+6)*di(m+lz+6,2))
                           d2(17)=d2(17)+dens(prim,n)*
     $                         (i2(m+lx  )*di(m+ly  ,1)*di(m+lz  ,2)+
     $                          i2(m+lx+1)*di(m+ly+1,1)*di(m+lz+1,2)+
     $                          i2(m+lx+2)*di(m+ly+2,1)*di(m+lz+2,2)+
     $                          i2(m+lx+3)*di(m+ly+3,1)*di(m+lz+3,2)+
     $                          i2(m+lx+4)*di(m+ly+4,1)*di(m+lz+4,2)+
     $                          i2(m+lx+5)*di(m+ly+5,1)*di(m+lz+5,2)+
     $                          i2(m+lx+6)*di(m+ly+6,1)*di(m+lz+6,2))
                           d2(18)=d2(18)+dens(prim,n)*
     $                         (i2(m+lx  )*i2(m+ly  )*d2i(m+lz  ,2)+
     $                          i2(m+lx+1)*i2(m+ly+1)*d2i(m+lz+1,2)+
     $                          i2(m+lx+2)*i2(m+ly+2)*d2i(m+lz+2,2)+
     $                          i2(m+lx+3)*i2(m+ly+3)*d2i(m+lz+3,2)+
     $                          i2(m+lx+4)*i2(m+ly+4)*d2i(m+lz+4,2)+
     $                          i2(m+lx+5)*i2(m+ly+5)*d2i(m+lz+5,2)+
     $                          i2(m+lx+6)*i2(m+ly+6)*d2i(m+lz+6,2))
                           d2(19)=d2(19)+dens(prim,n)*
     $                         (di(m+lx  ,2)*i2(m+ly  )*di(m+lz  ,2)+
     $                          di(m+lx+1,2)*i2(m+ly+1)*di(m+lz+1,2)+
     $                          di(m+lx+2,2)*i2(m+ly+2)*di(m+lz+2,2)+
     $                          di(m+lx+3,2)*i2(m+ly+3)*di(m+lz+3,2)+
     $                          di(m+lx+4,2)*i2(m+ly+4)*di(m+lz+4,2)+
     $                          di(m+lx+5,2)*i2(m+ly+5)*di(m+lz+5,2)+
     $                          di(m+lx+6,2)*i2(m+ly+6)*di(m+lz+6,2))
                           d2(20)=d2(20)+dens(prim,n)*
     $                         (i2(m+lx  )*di(m+ly  ,2)*di(m+lz  ,2)+
     $                          i2(m+lx+1)*di(m+ly+1,2)*di(m+lz+1,2)+
     $                          i2(m+lx+2)*di(m+ly+2,2)*di(m+lz+2,2)+
     $                          i2(m+lx+3)*di(m+ly+3,2)*di(m+lz+3,2)+
     $                          i2(m+lx+4)*di(m+ly+4,2)*di(m+lz+4,2)+
     $                          i2(m+lx+5)*di(m+ly+5,2)*di(m+lz+5,2)+
     $                          i2(m+lx+6)*di(m+ly+6,2)*di(m+lz+6,2))
                           d2(21)=d2(21)+dens(prim,n)*
     $                         (i2(m+lx  )*i2(m+ly  )*d2i(m+lz  ,3)+
     $                          i2(m+lx+1)*i2(m+ly+1)*d2i(m+lz+1,3)+
     $                          i2(m+lx+2)*i2(m+ly+2)*d2i(m+lz+2,3)+
     $                          i2(m+lx+3)*i2(m+ly+3)*d2i(m+lz+3,3)+
     $                          i2(m+lx+4)*i2(m+ly+4)*d2i(m+lz+4,3)+
     $                          i2(m+lx+5)*i2(m+ly+5)*d2i(m+lz+5,3)+
     $                          i2(m+lx+6)*i2(m+ly+6)*d2i(m+lz+6,3))
                           m=m+7
 1101                   continue
                     end if
                     if (npass.eq.1) then
                        m=1
                        do 1102 prim=1,nv
                           der(1,3)=der(1,3)+dens(prim,n)*
     $                        (di(m+lx  ,3)*i2(m+ly  )*i2(m+lz  )+
     $                         di(m+lx+1,3)*i2(m+ly+1)*i2(m+lz+1)+
     $                         di(m+lx+2,3)*i2(m+ly+2)*i2(m+lz+2)+
     $                         di(m+lx+3,3)*i2(m+ly+3)*i2(m+lz+3)+
     $                         di(m+lx+4,3)*i2(m+ly+4)*i2(m+lz+4)+
     $                         di(m+lx+5,3)*i2(m+ly+5)*i2(m+lz+5)+
     $                         di(m+lx+6,3)*i2(m+ly+6)*i2(m+lz+6))
                           der(2,3)=der(2,3)+dens(prim,n)*
     $                        (i2(m+lx  )*di(m+ly  ,3)*i2(m+lz  )+
     $                         i2(m+lx+1)*di(m+ly+1,3)*i2(m+lz+1)+
     $                         i2(m+lx+2)*di(m+ly+2,3)*i2(m+lz+2)+
     $                         i2(m+lx+3)*di(m+ly+3,3)*i2(m+lz+3)+
     $                         i2(m+lx+4)*di(m+ly+4,3)*i2(m+lz+4)+
     $                         i2(m+lx+5)*di(m+ly+5,3)*i2(m+lz+5)+
     $                         i2(m+lx+6)*di(m+ly+6,3)*i2(m+lz+6))
                           der(3,3)=der(3,3)+dens(prim,n)*
     $                        (i2(m+lx  )*i2(m+ly  )*di(m+lz  ,3)+
     $                         i2(m+lx+1)*i2(m+ly+1)*di(m+lz+1,3)+
     $                         i2(m+lx+2)*i2(m+ly+2)*di(m+lz+2,3)+
     $                         i2(m+lx+3)*i2(m+ly+3)*di(m+lz+3,3)+
     $                         i2(m+lx+4)*i2(m+ly+4)*di(m+lz+4,3)+
     $                         i2(m+lx+5)*i2(m+ly+5)*di(m+lz+5,3)+
     $                         i2(m+lx+6)*i2(m+ly+6)*di(m+lz+6,3))
                           d2(22)=d2(22)+dens(prim,n)*
     $                         (d2i(m+lx  ,4)*i2(m+ly  )*i2(m+lz  )+
     $                          d2i(m+lx+1,4)*i2(m+ly+1)*i2(m+lz+1)+
     $                          d2i(m+lx+2,4)*i2(m+ly+2)*i2(m+lz+2)+
     $                          d2i(m+lx+3,4)*i2(m+ly+3)*i2(m+lz+3)+
     $                          d2i(m+lx+4,4)*i2(m+ly+4)*i2(m+lz+4)+
     $                          d2i(m+lx+5,4)*i2(m+ly+5)*i2(m+lz+5)+
     $                          d2i(m+lx+6,4)*i2(m+ly+6)*i2(m+lz+6))
                           d2(23)=d2(23)+dens(prim,n)*
     $                         (di(m+lx  ,3)*di(m+ly  ,1)*i2(m+lz  )+
     $                          di(m+lx+1,3)*di(m+ly+1,1)*i2(m+lz+1)+
     $                          di(m+lx+2,3)*di(m+ly+2,1)*i2(m+lz+2)+
     $                          di(m+lx+3,3)*di(m+ly+3,1)*i2(m+lz+3)+
     $                          di(m+lx+4,3)*di(m+ly+4,1)*i2(m+lz+4)+
     $                          di(m+lx+5,3)*di(m+ly+5,1)*i2(m+lz+5)+
     $                          di(m+lx+6,3)*di(m+ly+6,1)*i2(m+lz+6))
                           d2(24)=d2(24)+dens(prim,n)*
     $                         (di(m+lx  ,3)*i2(m+ly  )*di(m+lz  ,1)+
     $                          di(m+lx+1,3)*i2(m+ly+1)*di(m+lz+1,1)+
     $                          di(m+lx+2,3)*i2(m+ly+2)*di(m+lz+2,1)+
     $                          di(m+lx+3,3)*i2(m+ly+3)*di(m+lz+3,1)+
     $                          di(m+lx+4,3)*i2(m+ly+4)*di(m+lz+4,1)+
     $                          di(m+lx+5,3)*i2(m+ly+5)*di(m+lz+5,1)+
     $                          di(m+lx+6,3)*i2(m+ly+6)*di(m+lz+6,1))
                           d2(25)=d2(25)+dens(prim,n)*
     $                         (d2i(m+lx  ,5)*i2(m+ly  )*i2(m+lz  )+
     $                          d2i(m+lx+1,5)*i2(m+ly+1)*i2(m+lz+1)+
     $                          d2i(m+lx+2,5)*i2(m+ly+2)*i2(m+lz+2)+
     $                          d2i(m+lx+3,5)*i2(m+ly+3)*i2(m+lz+3)+
     $                          d2i(m+lx+4,5)*i2(m+ly+4)*i2(m+lz+4)+
     $                          d2i(m+lx+5,5)*i2(m+ly+5)*i2(m+lz+5)+
     $                          d2i(m+lx+6,5)*i2(m+ly+6)*i2(m+lz+6))
                           d2(26)=d2(26)+dens(prim,n)*
     $                         (di(m+lx  ,3)*di(m+ly  ,2)*i2(m+lz  )+
     $                          di(m+lx+1,3)*di(m+ly+1,2)*i2(m+lz+1)+
     $                          di(m+lx+2,3)*di(m+ly+2,2)*i2(m+lz+2)+
     $                          di(m+lx+3,3)*di(m+ly+3,2)*i2(m+lz+3)+
     $                          di(m+lx+4,3)*di(m+ly+4,2)*i2(m+lz+4)+
     $                          di(m+lx+5,3)*di(m+ly+5,2)*i2(m+lz+5)+
     $                          di(m+lx+6,3)*di(m+ly+6,2)*i2(m+lz+6))
                           d2(27)=d2(27)+dens(prim,n)*
     $                         (di(m+lx  ,3)*i2(m+ly  )*di(m+lz  ,2)+
     $                          di(m+lx+1,3)*i2(m+ly+1)*di(m+lz+1,2)+
     $                          di(m+lx+2,3)*i2(m+ly+2)*di(m+lz+2,2)+
     $                          di(m+lx+3,3)*i2(m+ly+3)*di(m+lz+3,2)+
     $                          di(m+lx+4,3)*i2(m+ly+4)*di(m+lz+4,2)+
     $                          di(m+lx+5,3)*i2(m+ly+5)*di(m+lz+5,2)+
     $                          di(m+lx+6,3)*i2(m+ly+6)*di(m+lz+6,2))
                           d2(28)=d2(28)+dens(prim,n)*
     $                         (d2i(m+lx  ,6)*i2(m+ly  )*i2(m+lz  )+
     $                          d2i(m+lx+1,6)*i2(m+ly+1)*i2(m+lz+1)+
     $                          d2i(m+lx+2,6)*i2(m+ly+2)*i2(m+lz+2)+
     $                          d2i(m+lx+3,6)*i2(m+ly+3)*i2(m+lz+3)+
     $                          d2i(m+lx+4,6)*i2(m+ly+4)*i2(m+lz+4)+
     $                          d2i(m+lx+5,6)*i2(m+ly+5)*i2(m+lz+5)+
     $                          d2i(m+lx+6,6)*i2(m+ly+6)*i2(m+lz+6))
                           d2(29)=d2(29)+dens(prim,n)*
     $                         (di(m+lx  ,1)*di(m+ly  ,3)*i2(m+lz  )+
     $                          di(m+lx+1,1)*di(m+ly+1,3)*i2(m+lz+1)+
     $                          di(m+lx+2,1)*di(m+ly+2,3)*i2(m+lz+2)+
     $                          di(m+lx+3,1)*di(m+ly+3,3)*i2(m+lz+3)+
     $                          di(m+lx+4,1)*di(m+ly+4,3)*i2(m+lz+4)+
     $                          di(m+lx+5,1)*di(m+ly+5,3)*i2(m+lz+5)+
     $                          di(m+lx+6,1)*di(m+ly+6,3)*i2(m+lz+6))
                           d2(30)=d2(30)+dens(prim,n)*
     $                         (i2(m+lx  )*d2i(m+ly  ,4)*i2(m+lz  )+
     $                          i2(m+lx+1)*d2i(m+ly+1,4)*i2(m+lz+1)+
     $                          i2(m+lx+2)*d2i(m+ly+2,4)*i2(m+lz+2)+
     $                          i2(m+lx+3)*d2i(m+ly+3,4)*i2(m+lz+3)+
     $                          i2(m+lx+4)*d2i(m+ly+4,4)*i2(m+lz+4)+
     $                          i2(m+lx+5)*d2i(m+ly+5,4)*i2(m+lz+5)+
     $                          i2(m+lx+6)*d2i(m+ly+6,4)*i2(m+lz+6))
                           d2(31)=d2(31)+dens(prim,n)*
     $                         (i2(m+lx  )*di(m+ly  ,3)*di(m+lz  ,1)+
     $                          i2(m+lx+1)*di(m+ly+1,3)*di(m+lz+1,1)+
     $                          i2(m+lx+2)*di(m+ly+2,3)*di(m+lz+2,1)+
     $                          i2(m+lx+3)*di(m+ly+3,3)*di(m+lz+3,1)+
     $                          i2(m+lx+4)*di(m+ly+4,3)*di(m+lz+4,1)+
     $                          i2(m+lx+5)*di(m+ly+5,3)*di(m+lz+5,1)+
     $                          i2(m+lx+6)*di(m+ly+6,3)*di(m+lz+6,1))
                           d2(32)=d2(32)+dens(prim,n)*
     $                         (di(m+lx  ,2)*di(m+ly  ,3)*i2(m+lz  )+
     $                          di(m+lx+1,2)*di(m+ly+1,3)*i2(m+lz+1)+
     $                          di(m+lx+2,2)*di(m+ly+2,3)*i2(m+lz+2)+
     $                          di(m+lx+3,2)*di(m+ly+3,3)*i2(m+lz+3)+
     $                          di(m+lx+4,2)*di(m+ly+4,3)*i2(m+lz+4)+
     $                          di(m+lx+5,2)*di(m+ly+5,3)*i2(m+lz+5)+
     $                          di(m+lx+6,2)*di(m+ly+6,3)*i2(m+lz+6))
                           d2(33)=d2(33)+dens(prim,n)*
     $                         (i2(m+lx  )*d2i(m+ly  ,5)*i2(m+lz  )+
     $                          i2(m+lx+1)*d2i(m+ly+1,5)*i2(m+lz+1)+
     $                          i2(m+lx+2)*d2i(m+ly+2,5)*i2(m+lz+2)+
     $                          i2(m+lx+3)*d2i(m+ly+3,5)*i2(m+lz+3)+
     $                          i2(m+lx+4)*d2i(m+ly+4,5)*i2(m+lz+4)+
     $                          i2(m+lx+5)*d2i(m+ly+5,5)*i2(m+lz+5)+
     $                          i2(m+lx+6)*d2i(m+ly+6,5)*i2(m+lz+6))
                           d2(34)=d2(34)+dens(prim,n)*
     $                         (i2(m+lx  )*di(m+ly  ,3)*di(m+lz  ,2)+
     $                          i2(m+lx+1)*di(m+ly+1,3)*di(m+lz+1,2)+
     $                          i2(m+lx+2)*di(m+ly+2,3)*di(m+lz+2,2)+
     $                          i2(m+lx+3)*di(m+ly+3,3)*di(m+lz+3,2)+
     $                          i2(m+lx+4)*di(m+ly+4,3)*di(m+lz+4,2)+
     $                          i2(m+lx+5)*di(m+ly+5,3)*di(m+lz+5,2)+
     $                          i2(m+lx+6)*di(m+ly+6,3)*di(m+lz+6,2))
                           d2(35)=d2(35)+dens(prim,n)*
     $                         (di(m+lx  ,3)*di(m+ly  ,3)*i2(m+lz  )+
     $                          di(m+lx+1,3)*di(m+ly+1,3)*i2(m+lz+1)+
     $                          di(m+lx+2,3)*di(m+ly+2,3)*i2(m+lz+2)+
     $                          di(m+lx+3,3)*di(m+ly+3,3)*i2(m+lz+3)+
     $                          di(m+lx+4,3)*di(m+ly+4,3)*i2(m+lz+4)+
     $                          di(m+lx+5,3)*di(m+ly+5,3)*i2(m+lz+5)+
     $                          di(m+lx+6,3)*di(m+ly+6,3)*i2(m+lz+6))
                           d2(36)=d2(36)+dens(prim,n)*
     $                         (i2(m+lx  )*d2i(m+ly  ,6)*i2(m+lz  )+
     $                          i2(m+lx+1)*d2i(m+ly+1,6)*i2(m+lz+1)+
     $                          i2(m+lx+2)*d2i(m+ly+2,6)*i2(m+lz+2)+
     $                          i2(m+lx+3)*d2i(m+ly+3,6)*i2(m+lz+3)+
     $                          i2(m+lx+4)*d2i(m+ly+4,6)*i2(m+lz+4)+
     $                          i2(m+lx+5)*d2i(m+ly+5,6)*i2(m+lz+5)+
     $                          i2(m+lx+6)*d2i(m+ly+6,6)*i2(m+lz+6))
                           d2(37)=d2(37)+dens(prim,n)*
     $                         (di(m+lx  ,1)*i2(m+ly  )*di(m+lz  ,3)+
     $                          di(m+lx+1,1)*i2(m+ly+1)*di(m+lz+1,3)+
     $                          di(m+lx+2,1)*i2(m+ly+2)*di(m+lz+2,3)+
     $                          di(m+lx+3,1)*i2(m+ly+3)*di(m+lz+3,3)+
     $                          di(m+lx+4,1)*i2(m+ly+4)*di(m+lz+4,3)+
     $                          di(m+lx+5,1)*i2(m+ly+5)*di(m+lz+5,3)+
     $                          di(m+lx+6,1)*i2(m+ly+6)*di(m+lz+6,3))
                           d2(38)=d2(38)+dens(prim,n)*
     $                         (i2(m+lx  )*di(m+ly  ,1)*di(m+lz  ,3)+
     $                          i2(m+lx+1)*di(m+ly+1,1)*di(m+lz+1,3)+
     $                          i2(m+lx+2)*di(m+ly+2,1)*di(m+lz+2,3)+
     $                          i2(m+lx+3)*di(m+ly+3,1)*di(m+lz+3,3)+
     $                          i2(m+lx+4)*di(m+ly+4,1)*di(m+lz+4,3)+
     $                          i2(m+lx+5)*di(m+ly+5,1)*di(m+lz+5,3)+
     $                          i2(m+lx+6)*di(m+ly+6,1)*di(m+lz+6,3))
                           d2(39)=d2(39)+dens(prim,n)*
     $                         (i2(m+lx  )*i2(m+ly  )*d2i(m+lz  ,4)+
     $                          i2(m+lx+1)*i2(m+ly+1)*d2i(m+lz+1,4)+
     $                          i2(m+lx+2)*i2(m+ly+2)*d2i(m+lz+2,4)+
     $                          i2(m+lx+3)*i2(m+ly+3)*d2i(m+lz+3,4)+
     $                          i2(m+lx+4)*i2(m+ly+4)*d2i(m+lz+4,4)+
     $                          i2(m+lx+5)*i2(m+ly+5)*d2i(m+lz+5,4)+
     $                          i2(m+lx+6)*i2(m+ly+6)*d2i(m+lz+6,4))
                           d2(40)=d2(40)+dens(prim,n)*
     $                         (di(m+lx  ,2)*i2(m+ly  )*di(m+lz  ,3)+
     $                          di(m+lx+1,2)*i2(m+ly+1)*di(m+lz+1,3)+
     $                          di(m+lx+2,2)*i2(m+ly+2)*di(m+lz+2,3)+
     $                          di(m+lx+3,2)*i2(m+ly+3)*di(m+lz+3,3)+
     $                          di(m+lx+4,2)*i2(m+ly+4)*di(m+lz+4,3)+
     $                          di(m+lx+5,2)*i2(m+ly+5)*di(m+lz+5,3)+
     $                          di(m+lx+6,2)*i2(m+ly+6)*di(m+lz+6,3))
                           d2(41)=d2(41)+dens(prim,n)*
     $                         (i2(m+lx  )*di(m+ly  ,2)*di(m+lz  ,3)+
     $                          i2(m+lx+1)*di(m+ly+1,2)*di(m+lz+1,3)+
     $                          i2(m+lx+2)*di(m+ly+2,2)*di(m+lz+2,3)+
     $                          i2(m+lx+3)*di(m+ly+3,2)*di(m+lz+3,3)+
     $                          i2(m+lx+4)*di(m+ly+4,2)*di(m+lz+4,3)+
     $                          i2(m+lx+5)*di(m+ly+5,2)*di(m+lz+5,3)+
     $                          i2(m+lx+6)*di(m+ly+6,2)*di(m+lz+6,3))
                           d2(42)=d2(42)+dens(prim,n)*
     $                         (i2(m+lx  )*i2(m+ly  )*d2i(m+lz  ,5)+
     $                          i2(m+lx+1)*i2(m+ly+1)*d2i(m+lz+1,5)+
     $                          i2(m+lx+2)*i2(m+ly+2)*d2i(m+lz+2,5)+
     $                          i2(m+lx+3)*i2(m+ly+3)*d2i(m+lz+3,5)+
     $                          i2(m+lx+4)*i2(m+ly+4)*d2i(m+lz+4,5)+
     $                          i2(m+lx+5)*i2(m+ly+5)*d2i(m+lz+5,5)+
     $                          i2(m+lx+6)*i2(m+ly+6)*d2i(m+lz+6,5))
                           d2(43)=d2(43)+dens(prim,n)*
     $                         (di(m+lx  ,3)*i2(m+ly  )*di(m+lz  ,3)+
     $                          di(m+lx+1,3)*i2(m+ly+1)*di(m+lz+1,3)+
     $                          di(m+lx+2,3)*i2(m+ly+2)*di(m+lz+2,3)+
     $                          di(m+lx+3,3)*i2(m+ly+3)*di(m+lz+3,3)+
     $                          di(m+lx+4,3)*i2(m+ly+4)*di(m+lz+4,3)+
     $                          di(m+lx+5,3)*i2(m+ly+5)*di(m+lz+5,3)+
     $                          di(m+lx+6,3)*i2(m+ly+6)*di(m+lz+6,3))
                           d2(44)=d2(44)+dens(prim,n)*
     $                         (i2(m+lx  )*di(m+ly  ,3)*di(m+lz  ,3)+
     $                          i2(m+lx+1)*di(m+ly+1,3)*di(m+lz+1,3)+
     $                          i2(m+lx+2)*di(m+ly+2,3)*di(m+lz+2,3)+
     $                          i2(m+lx+3)*di(m+ly+3,3)*di(m+lz+3,3)+
     $                          i2(m+lx+4)*di(m+ly+4,3)*di(m+lz+4,3)+
     $                          i2(m+lx+5)*di(m+ly+5,3)*di(m+lz+5,3)+
     $                          i2(m+lx+6)*di(m+ly+6,3)*di(m+lz+6,3))
                           d2(45)=d2(45)+dens(prim,n)*
     $                         (i2(m+lx  )*i2(m+ly  )*d2i(m+lz  ,6)+
     $                          i2(m+lx+1)*i2(m+ly+1)*d2i(m+lz+1,6)+
     $                          i2(m+lx+2)*i2(m+ly+2)*d2i(m+lz+2,6)+
     $                          i2(m+lx+3)*i2(m+ly+3)*d2i(m+lz+3,6)+
     $                          i2(m+lx+4)*i2(m+ly+4)*d2i(m+lz+4,6)+
     $                          i2(m+lx+5)*i2(m+ly+5)*d2i(m+lz+5,6)+
     $                          i2(m+lx+6)*i2(m+ly+6)*d2i(m+lz+6,6))
                           m=m+7
 1102                   continue
                     end if
                     go to 15
c
c
 12               continue
                     m=1
                     do 1200 prim=1,nv
                        energy=energy+dens(prim,n)*
     $                      (i2(m+lx  )*i2(m+ly  )*i2(m+lz  )+
     $                       i2(m+lx+1)*i2(m+ly+1)*i2(m+lz+1)+
     $                       i2(m+lx+2)*i2(m+ly+2)*i2(m+lz+2)+
     $                       i2(m+lx+3)*i2(m+ly+3)*i2(m+lz+3)+
     $                       i2(m+lx+4)*i2(m+ly+4)*i2(m+lz+4)+
     $                       i2(m+lx+5)*i2(m+ly+5)*i2(m+lz+5)+
     $                       i2(m+lx+6)*i2(m+ly+6)*i2(m+lz+6)+
     $                       i2(m+lx+7)*i2(m+ly+7)*i2(m+lz+7))
                        der(1,1)=der(1,1)+dens(prim,n)*
     $                      (di(m+lx  ,1)*i2(m+ly  )*i2(m+lz  )+
     $                       di(m+lx+1,1)*i2(m+ly+1)*i2(m+lz+1)+
     $                       di(m+lx+2,1)*i2(m+ly+2)*i2(m+lz+2)+
     $                       di(m+lx+3,1)*i2(m+ly+3)*i2(m+lz+3)+
     $                       di(m+lx+4,1)*i2(m+ly+4)*i2(m+lz+4)+
     $                       di(m+lx+5,1)*i2(m+ly+5)*i2(m+lz+5)+
     $                       di(m+lx+6,1)*i2(m+ly+6)*i2(m+lz+6)+
     $                       di(m+lx+7,1)*i2(m+ly+7)*i2(m+lz+7))
                        der(2,1)=der(2,1)+dens(prim,n)*
     $                      (i2(m+lx  )*di(m+ly  ,1)*i2(m+lz  )+
     $                       i2(m+lx+1)*di(m+ly+1,1)*i2(m+lz+1)+
     $                       i2(m+lx+2)*di(m+ly+2,1)*i2(m+lz+2)+
     $                       i2(m+lx+3)*di(m+ly+3,1)*i2(m+lz+3)+
     $                       i2(m+lx+4)*di(m+ly+4,1)*i2(m+lz+4)+
     $                       i2(m+lx+5)*di(m+ly+5,1)*i2(m+lz+5)+
     $                       i2(m+lx+6)*di(m+ly+6,1)*i2(m+lz+6)+
     $                       i2(m+lx+7)*di(m+ly+7,1)*i2(m+lz+7))
                        der(3,1)=der(3,1)+dens(prim,n)*
     $                      (i2(m+lx  )*i2(m+ly  )*di(m+lz  ,1)+
     $                       i2(m+lx+1)*i2(m+ly+1)*di(m+lz+1,1)+
     $                       i2(m+lx+2)*i2(m+ly+2)*di(m+lz+2,1)+
     $                       i2(m+lx+3)*i2(m+ly+3)*di(m+lz+3,1)+
     $                       i2(m+lx+4)*i2(m+ly+4)*di(m+lz+4,1)+
     $                       i2(m+lx+5)*i2(m+ly+5)*di(m+lz+5,1)+
     $                       i2(m+lx+6)*i2(m+ly+6)*di(m+lz+6,1)+
     $                       i2(m+lx+7)*i2(m+ly+7)*di(m+lz+7,1))
                        d2(1)=d2(1)+dens(prim,n)*
     $                      (d2i(m+lx  ,1)*i2(m+ly  )*i2(m+lz  )+
     $                       d2i(m+lx+1,1)*i2(m+ly+1)*i2(m+lz+1)+
     $                       d2i(m+lx+2,1)*i2(m+ly+2)*i2(m+lz+2)+
     $                       d2i(m+lx+3,1)*i2(m+ly+3)*i2(m+lz+3)+
     $                       d2i(m+lx+4,1)*i2(m+ly+4)*i2(m+lz+4)+
     $                       d2i(m+lx+5,1)*i2(m+ly+5)*i2(m+lz+5)+
     $                       d2i(m+lx+6,1)*i2(m+ly+6)*i2(m+lz+6)+
     $                       d2i(m+lx+7,1)*i2(m+ly+7)*i2(m+lz+7))
                        d2(2)=d2(2)+dens(prim,n)*
     $                      (di(m+lx  ,1)*di(m+ly  ,1)*i2(m+lz  )+
     $                       di(m+lx+1,1)*di(m+ly+1,1)*i2(m+lz+1)+
     $                       di(m+lx+2,1)*di(m+ly+2,1)*i2(m+lz+2)+
     $                       di(m+lx+3,1)*di(m+ly+3,1)*i2(m+lz+3)+
     $                       di(m+lx+4,1)*di(m+ly+4,1)*i2(m+lz+4)+
     $                       di(m+lx+5,1)*di(m+ly+5,1)*i2(m+lz+5)+
     $                       di(m+lx+6,1)*di(m+ly+6,1)*i2(m+lz+6)+
     $                       di(m+lx+7,1)*di(m+ly+7,1)*i2(m+lz+7))
                        d2(3)=d2(3)+dens(prim,n)*
     $                      (i2(m+lx  )*d2i(m+ly  ,1)*i2(m+lz  )+
     $                       i2(m+lx+1)*d2i(m+ly+1,1)*i2(m+lz+1)+
     $                       i2(m+lx+2)*d2i(m+ly+2,1)*i2(m+lz+2)+
     $                       i2(m+lx+3)*d2i(m+ly+3,1)*i2(m+lz+3)+
     $                       i2(m+lx+4)*d2i(m+ly+4,1)*i2(m+lz+4)+
     $                       i2(m+lx+5)*d2i(m+ly+5,1)*i2(m+lz+5)+
     $                       i2(m+lx+6)*d2i(m+ly+6,1)*i2(m+lz+6)+
     $                       i2(m+lx+7)*d2i(m+ly+7,1)*i2(m+lz+7))
                        d2(4)=d2(4)+dens(prim,n)*
     $                      (di(m+lx  ,1)*i2(m+ly  )*di(m+lz  ,1)+
     $                       di(m+lx+1,1)*i2(m+ly+1)*di(m+lz+1,1)+
     $                       di(m+lx+2,1)*i2(m+ly+2)*di(m+lz+2,1)+
     $                       di(m+lx+3,1)*i2(m+ly+3)*di(m+lz+3,1)+
     $                       di(m+lx+4,1)*i2(m+ly+4)*di(m+lz+4,1)+
     $                       di(m+lx+5,1)*i2(m+ly+5)*di(m+lz+5,1)+
     $                       di(m+lx+6,1)*i2(m+ly+6)*di(m+lz+6,1)+
     $                       di(m+lx+7,1)*i2(m+ly+7)*di(m+lz+7,1))
                        d2(5)=d2(5)+dens(prim,n)*
     $                      (i2(m+lx  )*di(m+ly  ,1)*di(m+lz  ,1)+
     $                       i2(m+lx+1)*di(m+ly+1,1)*di(m+lz+1,1)+
     $                       i2(m+lx+2)*di(m+ly+2,1)*di(m+lz+2,1)+
     $                       i2(m+lx+3)*di(m+ly+3,1)*di(m+lz+3,1)+
     $                       i2(m+lx+4)*di(m+ly+4,1)*di(m+lz+4,1)+
     $                       i2(m+lx+5)*di(m+ly+5,1)*di(m+lz+5,1)+
     $                       i2(m+lx+6)*di(m+ly+6,1)*di(m+lz+6,1)+
     $                       i2(m+lx+7)*di(m+ly+7,1)*di(m+lz+7,1))
                        d2(6)=d2(6)+dens(prim,n)*
     $                      (i2(m+lx  )*i2(m+ly  )*d2i(m+lz  ,1)+
     $                       i2(m+lx+1)*i2(m+ly+1)*d2i(m+lz+1,1)+
     $                       i2(m+lx+2)*i2(m+ly+2)*d2i(m+lz+2,1)+
     $                       i2(m+lx+3)*i2(m+ly+3)*d2i(m+lz+3,1)+
     $                       i2(m+lx+4)*i2(m+ly+4)*d2i(m+lz+4,1)+
     $                       i2(m+lx+5)*i2(m+ly+5)*d2i(m+lz+5,1)+
     $                       i2(m+lx+6)*i2(m+ly+6)*d2i(m+lz+6,1)+
     $                       i2(m+lx+7)*i2(m+ly+7)*d2i(m+lz+7,1))
                        m=m+8
 1200                continue
                     if (npass.lt.4) then
                        m=1
                        do 1201 prim=1,nv
                           der(1,2)=der(1,2)+dens(prim,n)*
     $                        (di(m+lx  ,2)*i2(m+ly  )*i2(m+lz  )+
     $                         di(m+lx+1,2)*i2(m+ly+1)*i2(m+lz+1)+
     $                         di(m+lx+2,2)*i2(m+ly+2)*i2(m+lz+2)+
     $                         di(m+lx+3,2)*i2(m+ly+3)*i2(m+lz+3)+
     $                         di(m+lx+4,2)*i2(m+ly+4)*i2(m+lz+4)+
     $                         di(m+lx+5,2)*i2(m+ly+5)*i2(m+lz+5)+
     $                         di(m+lx+6,2)*i2(m+ly+6)*i2(m+lz+6)+
     $                         di(m+lx+7,2)*i2(m+ly+7)*i2(m+lz+7))
                           der(2,2)=der(2,2)+dens(prim,n)*
     $                        (i2(m+lx  )*di(m+ly  ,2)*i2(m+lz  )+
     $                         i2(m+lx+1)*di(m+ly+1,2)*i2(m+lz+1)+
     $                         i2(m+lx+2)*di(m+ly+2,2)*i2(m+lz+2)+
     $                         i2(m+lx+3)*di(m+ly+3,2)*i2(m+lz+3)+
     $                         i2(m+lx+4)*di(m+ly+4,2)*i2(m+lz+4)+
     $                         i2(m+lx+5)*di(m+ly+5,2)*i2(m+lz+5)+
     $                         i2(m+lx+6)*di(m+ly+6,2)*i2(m+lz+6)+
     $                         i2(m+lx+7)*di(m+ly+7,2)*i2(m+lz+7))
                           der(3,2)=der(3,2)+dens(prim,n)*
     $                        (i2(m+lx  )*i2(m+ly  )*di(m+lz  ,2)+
     $                         i2(m+lx+1)*i2(m+ly+1)*di(m+lz+1,2)+
     $                         i2(m+lx+2)*i2(m+ly+2)*di(m+lz+2,2)+
     $                         i2(m+lx+3)*i2(m+ly+3)*di(m+lz+3,2)+
     $                         i2(m+lx+4)*i2(m+ly+4)*di(m+lz+4,2)+
     $                         i2(m+lx+5)*i2(m+ly+5)*di(m+lz+5,2)+
     $                         i2(m+lx+6)*i2(m+ly+6)*di(m+lz+6,2)+
     $                         i2(m+lx+7)*i2(m+ly+7)*di(m+lz+7,2))
                           d2(7)=d2(7)+dens(prim,n)*
     $                         (d2i(m+lx  ,2)*i2(m+ly  )*i2(m+lz  )+
     $                          d2i(m+lx+1,2)*i2(m+ly+1)*i2(m+lz+1)+
     $                          d2i(m+lx+2,2)*i2(m+ly+2)*i2(m+lz+2)+
     $                          d2i(m+lx+3,2)*i2(m+ly+3)*i2(m+lz+3)+
     $                          d2i(m+lx+4,2)*i2(m+ly+4)*i2(m+lz+4)+
     $                          d2i(m+lx+5,2)*i2(m+ly+5)*i2(m+lz+5)+
     $                          d2i(m+lx+6,2)*i2(m+ly+6)*i2(m+lz+6)+
     $                          d2i(m+lx+7,2)*i2(m+ly+7)*i2(m+lz+7))
                           d2(8)=d2(8)+dens(prim,n)*
     $                         (di(m+lx  ,2)*di(m+ly  ,1)*i2(m+lz  )+
     $                          di(m+lx+1,2)*di(m+ly+1,1)*i2(m+lz+1)+
     $                          di(m+lx+2,2)*di(m+ly+2,1)*i2(m+lz+2)+
     $                          di(m+lx+3,2)*di(m+ly+3,1)*i2(m+lz+3)+
     $                          di(m+lx+4,2)*di(m+ly+4,1)*i2(m+lz+4)+
     $                          di(m+lx+5,2)*di(m+ly+5,1)*i2(m+lz+5)+
     $                          di(m+lx+6,2)*di(m+ly+6,1)*i2(m+lz+6)+
     $                          di(m+lx+7,2)*di(m+ly+7,1)*i2(m+lz+7))
                           d2(9)=d2(9)+dens(prim,n)*
     $                         (di(m+lx  ,2)*i2(m+ly  )*di(m+lz  ,1)+
     $                          di(m+lx+1,2)*i2(m+ly+1)*di(m+lz+1,1)+
     $                          di(m+lx+2,2)*i2(m+ly+2)*di(m+lz+2,1)+
     $                          di(m+lx+3,2)*i2(m+ly+3)*di(m+lz+3,1)+
     $                          di(m+lx+4,2)*i2(m+ly+4)*di(m+lz+4,1)+
     $                          di(m+lx+5,2)*i2(m+ly+5)*di(m+lz+5,1)+
     $                          di(m+lx+6,2)*i2(m+ly+6)*di(m+lz+6,1)+
     $                          di(m+lx+7,2)*i2(m+ly+7)*di(m+lz+7,1))
                           d2(10)=d2(10)+dens(prim,n)*
     $                         (d2i(m+lx  ,3)*i2(m+ly  )*i2(m+lz  )+
     $                          d2i(m+lx+1,3)*i2(m+ly+1)*i2(m+lz+1)+
     $                          d2i(m+lx+2,3)*i2(m+ly+2)*i2(m+lz+2)+
     $                          d2i(m+lx+3,3)*i2(m+ly+3)*i2(m+lz+3)+
     $                          d2i(m+lx+4,3)*i2(m+ly+4)*i2(m+lz+4)+
     $                          d2i(m+lx+5,3)*i2(m+ly+5)*i2(m+lz+5)+
     $                          d2i(m+lx+6,3)*i2(m+ly+6)*i2(m+lz+6)+
     $                          d2i(m+lx+7,3)*i2(m+ly+7)*i2(m+lz+7))
                           d2(11)=d2(11)+dens(prim,n)*
     $                         (di(m+lx  ,1)*di(m+ly  ,2)*i2(m+lz  )+
     $                          di(m+lx+1,1)*di(m+ly+1,2)*i2(m+lz+1)+
     $                          di(m+lx+2,1)*di(m+ly+2,2)*i2(m+lz+2)+
     $                          di(m+lx+3,1)*di(m+ly+3,2)*i2(m+lz+3)+
     $                          di(m+lx+4,1)*di(m+ly+4,2)*i2(m+lz+4)+
     $                          di(m+lx+5,1)*di(m+ly+5,2)*i2(m+lz+5)+
     $                          di(m+lx+6,1)*di(m+ly+6,2)*i2(m+lz+6)+
     $                          di(m+lx+7,1)*di(m+ly+7,2)*i2(m+lz+7))
                           d2(12)=d2(12)+dens(prim,n)*
     $                         (i2(m+lx  )*d2i(m+ly  ,2)*i2(m+lz  )+
     $                          i2(m+lx+1)*d2i(m+ly+1,2)*i2(m+lz+1)+
     $                          i2(m+lx+2)*d2i(m+ly+2,2)*i2(m+lz+2)+
     $                          i2(m+lx+3)*d2i(m+ly+3,2)*i2(m+lz+3)+
     $                          i2(m+lx+4)*d2i(m+ly+4,2)*i2(m+lz+4)+
     $                          i2(m+lx+5)*d2i(m+ly+5,2)*i2(m+lz+5)+
     $                          i2(m+lx+6)*d2i(m+ly+6,2)*i2(m+lz+6)+
     $                          i2(m+lx+7)*d2i(m+ly+7,2)*i2(m+lz+7))
                           d2(13)=d2(13)+dens(prim,n)*
     $                         (i2(m+lx  )*di(m+ly  ,2)*di(m+lz  ,1)+
     $                          i2(m+lx+1)*di(m+ly+1,2)*di(m+lz+1,1)+
     $                          i2(m+lx+2)*di(m+ly+2,2)*di(m+lz+2,1)+
     $                          i2(m+lx+3)*di(m+ly+3,2)*di(m+lz+3,1)+
     $                          i2(m+lx+4)*di(m+ly+4,2)*di(m+lz+4,1)+
     $                          i2(m+lx+5)*di(m+ly+5,2)*di(m+lz+5,1)+
     $                          i2(m+lx+6)*di(m+ly+6,2)*di(m+lz+6,1)+
     $                          i2(m+lx+7)*di(m+ly+7,2)*di(m+lz+7,1))
                           d2(14)=d2(14)+dens(prim,n)*
     $                         (di(m+lx  ,2)*di(m+ly  ,2)*i2(m+lz  )+
     $                          di(m+lx+1,2)*di(m+ly+1,2)*i2(m+lz+1)+
     $                          di(m+lx+2,2)*di(m+ly+2,2)*i2(m+lz+2)+
     $                          di(m+lx+3,2)*di(m+ly+3,2)*i2(m+lz+3)+
     $                          di(m+lx+4,2)*di(m+ly+4,2)*i2(m+lz+4)+
     $                          di(m+lx+5,2)*di(m+ly+5,2)*i2(m+lz+5)+
     $                          di(m+lx+6,2)*di(m+ly+6,2)*i2(m+lz+6)+
     $                          di(m+lx+7,2)*di(m+ly+7,2)*i2(m+lz+7))
                           d2(15)=d2(15)+dens(prim,n)*
     $                         (i2(m+lx  )*d2i(m+ly  ,3)*i2(m+lz  )+
     $                          i2(m+lx+1)*d2i(m+ly+1,3)*i2(m+lz+1)+
     $                          i2(m+lx+2)*d2i(m+ly+2,3)*i2(m+lz+2)+
     $                          i2(m+lx+3)*d2i(m+ly+3,3)*i2(m+lz+3)+
     $                          i2(m+lx+4)*d2i(m+ly+4,3)*i2(m+lz+4)+
     $                          i2(m+lx+5)*d2i(m+ly+5,3)*i2(m+lz+5)+
     $                          i2(m+lx+6)*d2i(m+ly+6,3)*i2(m+lz+6)+
     $                          i2(m+lx+7)*d2i(m+ly+7,3)*i2(m+lz+7))
                           d2(16)=d2(16)+dens(prim,n)*
     $                         (di(m+lx  ,1)*i2(m+ly  )*di(m+lz  ,2)+
     $                          di(m+lx+1,1)*i2(m+ly+1)*di(m+lz+1,2)+
     $                          di(m+lx+2,1)*i2(m+ly+2)*di(m+lz+2,2)+
     $                          di(m+lx+3,1)*i2(m+ly+3)*di(m+lz+3,2)+
     $                          di(m+lx+4,1)*i2(m+ly+4)*di(m+lz+4,2)+
     $                          di(m+lx+5,1)*i2(m+ly+5)*di(m+lz+5,2)+
     $                          di(m+lx+6,1)*i2(m+ly+6)*di(m+lz+6,2)+
     $                          di(m+lx+7,1)*i2(m+ly+7)*di(m+lz+7,2))
                           d2(17)=d2(17)+dens(prim,n)*
     $                         (i2(m+lx  )*di(m+ly  ,1)*di(m+lz  ,2)+
     $                          i2(m+lx+1)*di(m+ly+1,1)*di(m+lz+1,2)+
     $                          i2(m+lx+2)*di(m+ly+2,1)*di(m+lz+2,2)+
     $                          i2(m+lx+3)*di(m+ly+3,1)*di(m+lz+3,2)+
     $                          i2(m+lx+4)*di(m+ly+4,1)*di(m+lz+4,2)+
     $                          i2(m+lx+5)*di(m+ly+5,1)*di(m+lz+5,2)+
     $                          i2(m+lx+6)*di(m+ly+6,1)*di(m+lz+6,2)+
     $                          i2(m+lx+7)*di(m+ly+7,1)*di(m+lz+7,2))
                           d2(18)=d2(18)+dens(prim,n)*
     $                         (i2(m+lx  )*i2(m+ly  )*d2i(m+lz  ,2)+
     $                          i2(m+lx+1)*i2(m+ly+1)*d2i(m+lz+1,2)+
     $                          i2(m+lx+2)*i2(m+ly+2)*d2i(m+lz+2,2)+
     $                          i2(m+lx+3)*i2(m+ly+3)*d2i(m+lz+3,2)+
     $                          i2(m+lx+4)*i2(m+ly+4)*d2i(m+lz+4,2)+
     $                          i2(m+lx+5)*i2(m+ly+5)*d2i(m+lz+5,2)+
     $                          i2(m+lx+6)*i2(m+ly+6)*d2i(m+lz+6,2)+
     $                          i2(m+lx+7)*i2(m+ly+7)*d2i(m+lz+7,2))
                           d2(19)=d2(19)+dens(prim,n)*
     $                         (di(m+lx  ,2)*i2(m+ly  )*di(m+lz  ,2)+
     $                          di(m+lx+1,2)*i2(m+ly+1)*di(m+lz+1,2)+
     $                          di(m+lx+2,2)*i2(m+ly+2)*di(m+lz+2,2)+
     $                          di(m+lx+3,2)*i2(m+ly+3)*di(m+lz+3,2)+
     $                          di(m+lx+4,2)*i2(m+ly+4)*di(m+lz+4,2)+
     $                          di(m+lx+5,2)*i2(m+ly+5)*di(m+lz+5,2)+
     $                          di(m+lx+6,2)*i2(m+ly+6)*di(m+lz+6,2)+
     $                          di(m+lx+7,2)*i2(m+ly+7)*di(m+lz+7,2))
                           d2(20)=d2(20)+dens(prim,n)*
     $                         (i2(m+lx  )*di(m+ly  ,2)*di(m+lz  ,2)+
     $                          i2(m+lx+1)*di(m+ly+1,2)*di(m+lz+1,2)+
     $                          i2(m+lx+2)*di(m+ly+2,2)*di(m+lz+2,2)+
     $                          i2(m+lx+3)*di(m+ly+3,2)*di(m+lz+3,2)+
     $                          i2(m+lx+4)*di(m+ly+4,2)*di(m+lz+4,2)+
     $                          i2(m+lx+5)*di(m+ly+5,2)*di(m+lz+5,2)+
     $                          i2(m+lx+6)*di(m+ly+6,2)*di(m+lz+6,2)+
     $                          i2(m+lx+7)*di(m+ly+7,2)*di(m+lz+7,2))
                           d2(21)=d2(21)+dens(prim,n)*
     $                         (i2(m+lx  )*i2(m+ly  )*d2i(m+lz  ,3)+
     $                          i2(m+lx+1)*i2(m+ly+1)*d2i(m+lz+1,3)+
     $                          i2(m+lx+2)*i2(m+ly+2)*d2i(m+lz+2,3)+
     $                          i2(m+lx+3)*i2(m+ly+3)*d2i(m+lz+3,3)+
     $                          i2(m+lx+4)*i2(m+ly+4)*d2i(m+lz+4,3)+
     $                          i2(m+lx+5)*i2(m+ly+5)*d2i(m+lz+5,3)+
     $                          i2(m+lx+6)*i2(m+ly+6)*d2i(m+lz+6,3)+
     $                          i2(m+lx+7)*i2(m+ly+7)*d2i(m+lz+7,3))
                           m=m+8
 1201                   continue
                     end if
                     if (npass.eq.1) then
                        m=1
                        do 1202 prim=1,nv
                           der(1,3)=der(1,3)+dens(prim,n)*
     $                        (di(m+lx  ,3)*i2(m+ly  )*i2(m+lz  )+
     $                         di(m+lx+1,3)*i2(m+ly+1)*i2(m+lz+1)+
     $                         di(m+lx+2,3)*i2(m+ly+2)*i2(m+lz+2)+
     $                         di(m+lx+3,3)*i2(m+ly+3)*i2(m+lz+3)+
     $                         di(m+lx+4,3)*i2(m+ly+4)*i2(m+lz+4)+
     $                         di(m+lx+5,3)*i2(m+ly+5)*i2(m+lz+5)+
     $                         di(m+lx+6,3)*i2(m+ly+6)*i2(m+lz+6)+
     $                         di(m+lx+7,3)*i2(m+ly+7)*i2(m+lz+7))
                           der(2,3)=der(2,3)+dens(prim,n)*
     $                        (i2(m+lx  )*di(m+ly  ,3)*i2(m+lz  )+
     $                         i2(m+lx+1)*di(m+ly+1,3)*i2(m+lz+1)+
     $                         i2(m+lx+2)*di(m+ly+2,3)*i2(m+lz+2)+
     $                         i2(m+lx+3)*di(m+ly+3,3)*i2(m+lz+3)+
     $                         i2(m+lx+4)*di(m+ly+4,3)*i2(m+lz+4)+
     $                         i2(m+lx+5)*di(m+ly+5,3)*i2(m+lz+5)+
     $                         i2(m+lx+6)*di(m+ly+6,3)*i2(m+lz+6)+
     $                         i2(m+lx+7)*di(m+ly+7,3)*i2(m+lz+7))
                           der(3,3)=der(3,3)+dens(prim,n)*
     $                        (i2(m+lx  )*i2(m+ly  )*di(m+lz  ,3)+
     $                         i2(m+lx+1)*i2(m+ly+1)*di(m+lz+1,3)+
     $                         i2(m+lx+2)*i2(m+ly+2)*di(m+lz+2,3)+
     $                         i2(m+lx+3)*i2(m+ly+3)*di(m+lz+3,3)+
     $                         i2(m+lx+4)*i2(m+ly+4)*di(m+lz+4,3)+
     $                         i2(m+lx+5)*i2(m+ly+5)*di(m+lz+5,3)+
     $                         i2(m+lx+6)*i2(m+ly+6)*di(m+lz+6,3)+
     $                         i2(m+lx+7)*i2(m+ly+7)*di(m+lz+7,3))
                           d2(22)=d2(22)+dens(prim,n)*
     $                         (d2i(m+lx  ,4)*i2(m+ly  )*i2(m+lz  )+
     $                          d2i(m+lx+1,4)*i2(m+ly+1)*i2(m+lz+1)+
     $                          d2i(m+lx+2,4)*i2(m+ly+2)*i2(m+lz+2)+
     $                          d2i(m+lx+3,4)*i2(m+ly+3)*i2(m+lz+3)+
     $                          d2i(m+lx+4,4)*i2(m+ly+4)*i2(m+lz+4)+
     $                          d2i(m+lx+5,4)*i2(m+ly+5)*i2(m+lz+5)+
     $                          d2i(m+lx+6,4)*i2(m+ly+6)*i2(m+lz+6)+
     $                          d2i(m+lx+7,4)*i2(m+ly+7)*i2(m+lz+7))
                           d2(23)=d2(23)+dens(prim,n)*
     $                         (di(m+lx  ,3)*di(m+ly  ,1)*i2(m+lz  )+
     $                          di(m+lx+1,3)*di(m+ly+1,1)*i2(m+lz+1)+
     $                          di(m+lx+2,3)*di(m+ly+2,1)*i2(m+lz+2)+
     $                          di(m+lx+3,3)*di(m+ly+3,1)*i2(m+lz+3)+
     $                          di(m+lx+4,3)*di(m+ly+4,1)*i2(m+lz+4)+
     $                          di(m+lx+5,3)*di(m+ly+5,1)*i2(m+lz+5)+
     $                          di(m+lx+6,3)*di(m+ly+6,1)*i2(m+lz+6)+
     $                          di(m+lx+7,3)*di(m+ly+7,1)*i2(m+lz+7))
                           d2(24)=d2(24)+dens(prim,n)*
     $                         (di(m+lx  ,3)*i2(m+ly  )*di(m+lz  ,1)+
     $                          di(m+lx+1,3)*i2(m+ly+1)*di(m+lz+1,1)+
     $                          di(m+lx+2,3)*i2(m+ly+2)*di(m+lz+2,1)+
     $                          di(m+lx+3,3)*i2(m+ly+3)*di(m+lz+3,1)+
     $                          di(m+lx+4,3)*i2(m+ly+4)*di(m+lz+4,1)+
     $                          di(m+lx+5,3)*i2(m+ly+5)*di(m+lz+5,1)+
     $                          di(m+lx+6,3)*i2(m+ly+6)*di(m+lz+6,1)+
     $                          di(m+lx+7,3)*i2(m+ly+7)*di(m+lz+7,1))
                           d2(25)=d2(25)+dens(prim,n)*
     $                         (d2i(m+lx  ,5)*i2(m+ly  )*i2(m+lz  )+
     $                          d2i(m+lx+1,5)*i2(m+ly+1)*i2(m+lz+1)+
     $                          d2i(m+lx+2,5)*i2(m+ly+2)*i2(m+lz+2)+
     $                          d2i(m+lx+3,5)*i2(m+ly+3)*i2(m+lz+3)+
     $                          d2i(m+lx+4,5)*i2(m+ly+4)*i2(m+lz+4)+
     $                          d2i(m+lx+5,5)*i2(m+ly+5)*i2(m+lz+5)+
     $                          d2i(m+lx+6,5)*i2(m+ly+6)*i2(m+lz+6)+
     $                          d2i(m+lx+7,5)*i2(m+ly+7)*i2(m+lz+7))
                           d2(26)=d2(26)+dens(prim,n)*
     $                         (di(m+lx  ,3)*di(m+ly  ,2)*i2(m+lz  )+
     $                          di(m+lx+1,3)*di(m+ly+1,2)*i2(m+lz+1)+
     $                          di(m+lx+2,3)*di(m+ly+2,2)*i2(m+lz+2)+
     $                          di(m+lx+3,3)*di(m+ly+3,2)*i2(m+lz+3)+
     $                          di(m+lx+4,3)*di(m+ly+4,2)*i2(m+lz+4)+
     $                          di(m+lx+5,3)*di(m+ly+5,2)*i2(m+lz+5)+
     $                          di(m+lx+6,3)*di(m+ly+6,2)*i2(m+lz+6)+
     $                          di(m+lx+7,3)*di(m+ly+7,2)*i2(m+lz+7))
                           d2(27)=d2(27)+dens(prim,n)*
     $                         (di(m+lx  ,3)*i2(m+ly  )*di(m+lz  ,2)+
     $                          di(m+lx+1,3)*i2(m+ly+1)*di(m+lz+1,2)+
     $                          di(m+lx+2,3)*i2(m+ly+2)*di(m+lz+2,2)+
     $                          di(m+lx+3,3)*i2(m+ly+3)*di(m+lz+3,2)+
     $                          di(m+lx+4,3)*i2(m+ly+4)*di(m+lz+4,2)+
     $                          di(m+lx+5,3)*i2(m+ly+5)*di(m+lz+5,2)+
     $                          di(m+lx+6,3)*i2(m+ly+6)*di(m+lz+6,2)+
     $                          di(m+lx+7,3)*i2(m+ly+7)*di(m+lz+7,2))
                           d2(28)=d2(28)+dens(prim,n)*
     $                         (d2i(m+lx  ,6)*i2(m+ly  )*i2(m+lz  )+
     $                          d2i(m+lx+1,6)*i2(m+ly+1)*i2(m+lz+1)+
     $                          d2i(m+lx+2,6)*i2(m+ly+2)*i2(m+lz+2)+
     $                          d2i(m+lx+3,6)*i2(m+ly+3)*i2(m+lz+3)+
     $                          d2i(m+lx+4,6)*i2(m+ly+4)*i2(m+lz+4)+
     $                          d2i(m+lx+5,6)*i2(m+ly+5)*i2(m+lz+5)+
     $                          d2i(m+lx+6,6)*i2(m+ly+6)*i2(m+lz+6)+
     $                          d2i(m+lx+7,6)*i2(m+ly+7)*i2(m+lz+7))
                           d2(29)=d2(29)+dens(prim,n)*
     $                         (di(m+lx  ,1)*di(m+ly  ,3)*i2(m+lz  )+
     $                          di(m+lx+1,1)*di(m+ly+1,3)*i2(m+lz+1)+
     $                          di(m+lx+2,1)*di(m+ly+2,3)*i2(m+lz+2)+
     $                          di(m+lx+3,1)*di(m+ly+3,3)*i2(m+lz+3)+
     $                          di(m+lx+4,1)*di(m+ly+4,3)*i2(m+lz+4)+
     $                          di(m+lx+5,1)*di(m+ly+5,3)*i2(m+lz+5)+
     $                          di(m+lx+6,1)*di(m+ly+6,3)*i2(m+lz+6)+
     $                          di(m+lx+7,1)*di(m+ly+7,3)*i2(m+lz+7))
                           d2(30)=d2(30)+dens(prim,n)*
     $                         (i2(m+lx  )*d2i(m+ly  ,4)*i2(m+lz  )+
     $                          i2(m+lx+1)*d2i(m+ly+1,4)*i2(m+lz+1)+
     $                          i2(m+lx+2)*d2i(m+ly+2,4)*i2(m+lz+2)+
     $                          i2(m+lx+3)*d2i(m+ly+3,4)*i2(m+lz+3)+
     $                          i2(m+lx+4)*d2i(m+ly+4,4)*i2(m+lz+4)+
     $                          i2(m+lx+5)*d2i(m+ly+5,4)*i2(m+lz+5)+
     $                          i2(m+lx+6)*d2i(m+ly+6,4)*i2(m+lz+6)+
     $                          i2(m+lx+7)*d2i(m+ly+7,4)*i2(m+lz+7))
                           d2(31)=d2(31)+dens(prim,n)*
     $                         (i2(m+lx  )*di(m+ly  ,3)*di(m+lz  ,1)+
     $                          i2(m+lx+1)*di(m+ly+1,3)*di(m+lz+1,1)+
     $                          i2(m+lx+2)*di(m+ly+2,3)*di(m+lz+2,1)+
     $                          i2(m+lx+3)*di(m+ly+3,3)*di(m+lz+3,1)+
     $                          i2(m+lx+4)*di(m+ly+4,3)*di(m+lz+4,1)+
     $                          i2(m+lx+5)*di(m+ly+5,3)*di(m+lz+5,1)+
     $                          i2(m+lx+6)*di(m+ly+6,3)*di(m+lz+6,1)+
     $                          i2(m+lx+7)*di(m+ly+7,3)*di(m+lz+7,1))
                           d2(32)=d2(32)+dens(prim,n)*
     $                         (di(m+lx  ,2)*di(m+ly  ,3)*i2(m+lz  )+
     $                          di(m+lx+1,2)*di(m+ly+1,3)*i2(m+lz+1)+
     $                          di(m+lx+2,2)*di(m+ly+2,3)*i2(m+lz+2)+
     $                          di(m+lx+3,2)*di(m+ly+3,3)*i2(m+lz+3)+
     $                          di(m+lx+4,2)*di(m+ly+4,3)*i2(m+lz+4)+
     $                          di(m+lx+5,2)*di(m+ly+5,3)*i2(m+lz+5)+
     $                          di(m+lx+6,2)*di(m+ly+6,3)*i2(m+lz+6)+
     $                          di(m+lx+7,2)*di(m+ly+7,3)*i2(m+lz+7))
                           d2(33)=d2(33)+dens(prim,n)*
     $                         (i2(m+lx  )*d2i(m+ly  ,5)*i2(m+lz  )+
     $                          i2(m+lx+1)*d2i(m+ly+1,5)*i2(m+lz+1)+
     $                          i2(m+lx+2)*d2i(m+ly+2,5)*i2(m+lz+2)+
     $                          i2(m+lx+3)*d2i(m+ly+3,5)*i2(m+lz+3)+
     $                          i2(m+lx+4)*d2i(m+ly+4,5)*i2(m+lz+4)+
     $                          i2(m+lx+5)*d2i(m+ly+5,5)*i2(m+lz+5)+
     $                          i2(m+lx+6)*d2i(m+ly+6,5)*i2(m+lz+6)+
     $                          i2(m+lx+7)*d2i(m+ly+7,5)*i2(m+lz+7))
                           d2(34)=d2(34)+dens(prim,n)*
     $                         (i2(m+lx  )*di(m+ly  ,3)*di(m+lz  ,2)+
     $                          i2(m+lx+1)*di(m+ly+1,3)*di(m+lz+1,2)+
     $                          i2(m+lx+2)*di(m+ly+2,3)*di(m+lz+2,2)+
     $                          i2(m+lx+3)*di(m+ly+3,3)*di(m+lz+3,2)+
     $                          i2(m+lx+4)*di(m+ly+4,3)*di(m+lz+4,2)+
     $                          i2(m+lx+5)*di(m+ly+5,3)*di(m+lz+5,2)+
     $                          i2(m+lx+6)*di(m+ly+6,3)*di(m+lz+6,2)+
     $                          i2(m+lx+7)*di(m+ly+7,3)*di(m+lz+7,2))
                           d2(35)=d2(35)+dens(prim,n)*
     $                         (di(m+lx  ,3)*di(m+ly  ,3)*i2(m+lz  )+
     $                          di(m+lx+1,3)*di(m+ly+1,3)*i2(m+lz+1)+
     $                          di(m+lx+2,3)*di(m+ly+2,3)*i2(m+lz+2)+
     $                          di(m+lx+3,3)*di(m+ly+3,3)*i2(m+lz+3)+
     $                          di(m+lx+4,3)*di(m+ly+4,3)*i2(m+lz+4)+
     $                          di(m+lx+5,3)*di(m+ly+5,3)*i2(m+lz+5)+
     $                          di(m+lx+6,3)*di(m+ly+6,3)*i2(m+lz+6)+
     $                          di(m+lx+7,3)*di(m+ly+7,3)*i2(m+lz+7))
                           d2(36)=d2(36)+dens(prim,n)*
     $                         (i2(m+lx  )*d2i(m+ly  ,6)*i2(m+lz  )+
     $                          i2(m+lx+1)*d2i(m+ly+1,6)*i2(m+lz+1)+
     $                          i2(m+lx+2)*d2i(m+ly+2,6)*i2(m+lz+2)+
     $                          i2(m+lx+3)*d2i(m+ly+3,6)*i2(m+lz+3)+
     $                          i2(m+lx+4)*d2i(m+ly+4,6)*i2(m+lz+4)+
     $                          i2(m+lx+5)*d2i(m+ly+5,6)*i2(m+lz+5)+
     $                          i2(m+lx+6)*d2i(m+ly+6,6)*i2(m+lz+6)+
     $                          i2(m+lx+7)*d2i(m+ly+7,6)*i2(m+lz+7))
                           d2(37)=d2(37)+dens(prim,n)*
     $                         (di(m+lx  ,1)*i2(m+ly  )*di(m+lz  ,3)+
     $                          di(m+lx+1,1)*i2(m+ly+1)*di(m+lz+1,3)+
     $                          di(m+lx+2,1)*i2(m+ly+2)*di(m+lz+2,3)+
     $                          di(m+lx+3,1)*i2(m+ly+3)*di(m+lz+3,3)+
     $                          di(m+lx+4,1)*i2(m+ly+4)*di(m+lz+4,3)+
     $                          di(m+lx+5,1)*i2(m+ly+5)*di(m+lz+5,3)+
     $                          di(m+lx+6,1)*i2(m+ly+6)*di(m+lz+6,3)+
     $                          di(m+lx+7,1)*i2(m+ly+7)*di(m+lz+7,3))
                           d2(38)=d2(38)+dens(prim,n)*
     $                         (i2(m+lx  )*di(m+ly  ,1)*di(m+lz  ,3)+
     $                          i2(m+lx+1)*di(m+ly+1,1)*di(m+lz+1,3)+
     $                          i2(m+lx+2)*di(m+ly+2,1)*di(m+lz+2,3)+
     $                          i2(m+lx+3)*di(m+ly+3,1)*di(m+lz+3,3)+
     $                          i2(m+lx+4)*di(m+ly+4,1)*di(m+lz+4,3)+
     $                          i2(m+lx+5)*di(m+ly+5,1)*di(m+lz+5,3)+
     $                          i2(m+lx+6)*di(m+ly+6,1)*di(m+lz+6,3)+
     $                          i2(m+lx+7)*di(m+ly+7,1)*di(m+lz+7,3))
                           d2(39)=d2(39)+dens(prim,n)*
     $                         (i2(m+lx  )*i2(m+ly  )*d2i(m+lz  ,4)+
     $                          i2(m+lx+1)*i2(m+ly+1)*d2i(m+lz+1,4)+
     $                          i2(m+lx+2)*i2(m+ly+2)*d2i(m+lz+2,4)+
     $                          i2(m+lx+3)*i2(m+ly+3)*d2i(m+lz+3,4)+
     $                          i2(m+lx+4)*i2(m+ly+4)*d2i(m+lz+4,4)+
     $                          i2(m+lx+5)*i2(m+ly+5)*d2i(m+lz+5,4)+
     $                          i2(m+lx+6)*i2(m+ly+6)*d2i(m+lz+6,4)+
     $                          i2(m+lx+7)*i2(m+ly+7)*d2i(m+lz+7,4))
                           d2(40)=d2(40)+dens(prim,n)*
     $                         (di(m+lx  ,2)*i2(m+ly  )*di(m+lz  ,3)+
     $                          di(m+lx+1,2)*i2(m+ly+1)*di(m+lz+1,3)+
     $                          di(m+lx+2,2)*i2(m+ly+2)*di(m+lz+2,3)+
     $                          di(m+lx+3,2)*i2(m+ly+3)*di(m+lz+3,3)+
     $                          di(m+lx+4,2)*i2(m+ly+4)*di(m+lz+4,3)+
     $                          di(m+lx+5,2)*i2(m+ly+5)*di(m+lz+5,3)+
     $                          di(m+lx+6,2)*i2(m+ly+6)*di(m+lz+6,3)+
     $                          di(m+lx+7,2)*i2(m+ly+7)*di(m+lz+7,3))
                           d2(41)=d2(41)+dens(prim,n)*
     $                         (i2(m+lx  )*di(m+ly  ,2)*di(m+lz  ,3)+
     $                          i2(m+lx+1)*di(m+ly+1,2)*di(m+lz+1,3)+
     $                          i2(m+lx+2)*di(m+ly+2,2)*di(m+lz+2,3)+
     $                          i2(m+lx+3)*di(m+ly+3,2)*di(m+lz+3,3)+
     $                          i2(m+lx+4)*di(m+ly+4,2)*di(m+lz+4,3)+
     $                          i2(m+lx+5)*di(m+ly+5,2)*di(m+lz+5,3)+
     $                          i2(m+lx+6)*di(m+ly+6,2)*di(m+lz+6,3)+
     $                          i2(m+lx+7)*di(m+ly+7,2)*di(m+lz+7,3))
                           d2(42)=d2(42)+dens(prim,n)*
     $                         (i2(m+lx  )*i2(m+ly  )*d2i(m+lz  ,5)+
     $                          i2(m+lx+1)*i2(m+ly+1)*d2i(m+lz+1,5)+
     $                          i2(m+lx+2)*i2(m+ly+2)*d2i(m+lz+2,5)+
     $                          i2(m+lx+3)*i2(m+ly+3)*d2i(m+lz+3,5)+
     $                          i2(m+lx+4)*i2(m+ly+4)*d2i(m+lz+4,5)+
     $                          i2(m+lx+5)*i2(m+ly+5)*d2i(m+lz+5,5)+
     $                          i2(m+lx+6)*i2(m+ly+6)*d2i(m+lz+6,5)+
     $                          i2(m+lx+7)*i2(m+ly+7)*d2i(m+lz+7,5))
                           d2(43)=d2(43)+dens(prim,n)*
     $                         (di(m+lx  ,3)*i2(m+ly  )*di(m+lz  ,3)+
     $                          di(m+lx+1,3)*i2(m+ly+1)*di(m+lz+1,3)+
     $                          di(m+lx+2,3)*i2(m+ly+2)*di(m+lz+2,3)+
     $                          di(m+lx+3,3)*i2(m+ly+3)*di(m+lz+3,3)+
     $                          di(m+lx+4,3)*i2(m+ly+4)*di(m+lz+4,3)+
     $                          di(m+lx+5,3)*i2(m+ly+5)*di(m+lz+5,3)+
     $                          di(m+lx+6,3)*i2(m+ly+6)*di(m+lz+6,3)+
     $                          di(m+lx+7,3)*i2(m+ly+7)*di(m+lz+7,3))
                           d2(44)=d2(44)+dens(prim,n)*
     $                         (i2(m+lx  )*di(m+ly  ,3)*di(m+lz  ,3)+
     $                          i2(m+lx+1)*di(m+ly+1,3)*di(m+lz+1,3)+
     $                          i2(m+lx+2)*di(m+ly+2,3)*di(m+lz+2,3)+
     $                          i2(m+lx+3)*di(m+ly+3,3)*di(m+lz+3,3)+
     $                          i2(m+lx+4)*di(m+ly+4,3)*di(m+lz+4,3)+
     $                          i2(m+lx+5)*di(m+ly+5,3)*di(m+lz+5,3)+
     $                          i2(m+lx+6)*di(m+ly+6,3)*di(m+lz+6,3)+
     $                          i2(m+lx+7)*di(m+ly+7,3)*di(m+lz+7,3))
                           d2(45)=d2(45)+dens(prim,n)*
     $                         (i2(m+lx  )*i2(m+ly  )*d2i(m+lz  ,6)+
     $                          i2(m+lx+1)*i2(m+ly+1)*d2i(m+lz+1,6)+
     $                          i2(m+lx+2)*i2(m+ly+2)*d2i(m+lz+2,6)+
     $                          i2(m+lx+3)*i2(m+ly+3)*d2i(m+lz+3,6)+
     $                          i2(m+lx+4)*i2(m+ly+4)*d2i(m+lz+4,6)+
     $                          i2(m+lx+5)*i2(m+ly+5)*d2i(m+lz+5,6)+
     $                          i2(m+lx+6)*i2(m+ly+6)*d2i(m+lz+6,6)+
     $                          i2(m+lx+7)*i2(m+ly+7)*d2i(m+lz+7,6))
                           m=m+8
 1202                   continue
                     end if
                     go to 15
c
c
 13               continue
                     m=1
                     do 1300 prim=1,nv
                        energy=energy+dens(prim,n)*
     $                      (i2(m+lx  )*i2(m+ly  )*i2(m+lz  )+
     $                       i2(m+lx+1)*i2(m+ly+1)*i2(m+lz+1)+
     $                       i2(m+lx+2)*i2(m+ly+2)*i2(m+lz+2)+
     $                       i2(m+lx+3)*i2(m+ly+3)*i2(m+lz+3)+
     $                       i2(m+lx+4)*i2(m+ly+4)*i2(m+lz+4)+
     $                       i2(m+lx+5)*i2(m+ly+5)*i2(m+lz+5)+
     $                       i2(m+lx+6)*i2(m+ly+6)*i2(m+lz+6)+
     $                       i2(m+lx+7)*i2(m+ly+7)*i2(m+lz+7)+
     $                       i2(m+lx+8)*i2(m+ly+8)*i2(m+lz+8))
                        der(1,1)=der(1,1)+dens(prim,n)*
     $                      (di(m+lx  ,1)*i2(m+ly  )*i2(m+lz  )+
     $                       di(m+lx+1,1)*i2(m+ly+1)*i2(m+lz+1)+
     $                       di(m+lx+2,1)*i2(m+ly+2)*i2(m+lz+2)+
     $                       di(m+lx+3,1)*i2(m+ly+3)*i2(m+lz+3)+
     $                       di(m+lx+4,1)*i2(m+ly+4)*i2(m+lz+4)+
     $                       di(m+lx+5,1)*i2(m+ly+5)*i2(m+lz+5)+
     $                       di(m+lx+6,1)*i2(m+ly+6)*i2(m+lz+6)+
     $                       di(m+lx+7,1)*i2(m+ly+7)*i2(m+lz+7)+
     $                       di(m+lx+8,1)*i2(m+ly+8)*i2(m+lz+8))
                        der(2,1)=der(2,1)+dens(prim,n)*
     $                      (i2(m+lx  )*di(m+ly  ,1)*i2(m+lz  )+
     $                       i2(m+lx+1)*di(m+ly+1,1)*i2(m+lz+1)+
     $                       i2(m+lx+2)*di(m+ly+2,1)*i2(m+lz+2)+
     $                       i2(m+lx+3)*di(m+ly+3,1)*i2(m+lz+3)+
     $                       i2(m+lx+4)*di(m+ly+4,1)*i2(m+lz+4)+
     $                       i2(m+lx+5)*di(m+ly+5,1)*i2(m+lz+5)+
     $                       i2(m+lx+6)*di(m+ly+6,1)*i2(m+lz+6)+
     $                       i2(m+lx+7)*di(m+ly+7,1)*i2(m+lz+7)+
     $                       i2(m+lx+8)*di(m+ly+8,1)*i2(m+lz+8))
                        der(3,1)=der(3,1)+dens(prim,n)*
     $                      (i2(m+lx  )*i2(m+ly  )*di(m+lz  ,1)+
     $                       i2(m+lx+1)*i2(m+ly+1)*di(m+lz+1,1)+
     $                       i2(m+lx+2)*i2(m+ly+2)*di(m+lz+2,1)+
     $                       i2(m+lx+3)*i2(m+ly+3)*di(m+lz+3,1)+
     $                       i2(m+lx+4)*i2(m+ly+4)*di(m+lz+4,1)+
     $                       i2(m+lx+5)*i2(m+ly+5)*di(m+lz+5,1)+
     $                       i2(m+lx+6)*i2(m+ly+6)*di(m+lz+6,1)+
     $                       i2(m+lx+7)*i2(m+ly+7)*di(m+lz+7,1)+
     $                       i2(m+lx+8)*i2(m+ly+8)*di(m+lz+8,1))
                        d2(1)=d2(1)+dens(prim,n)*
     $                      (d2i(m+lx  ,1)*i2(m+ly  )*i2(m+lz  )+
     $                       d2i(m+lx+1,1)*i2(m+ly+1)*i2(m+lz+1)+
     $                       d2i(m+lx+2,1)*i2(m+ly+2)*i2(m+lz+2)+
     $                       d2i(m+lx+3,1)*i2(m+ly+3)*i2(m+lz+3)+
     $                       d2i(m+lx+4,1)*i2(m+ly+4)*i2(m+lz+4)+
     $                       d2i(m+lx+5,1)*i2(m+ly+5)*i2(m+lz+5)+
     $                       d2i(m+lx+6,1)*i2(m+ly+6)*i2(m+lz+6)+
     $                       d2i(m+lx+7,1)*i2(m+ly+7)*i2(m+lz+7)+
     $                       d2i(m+lx+8,1)*i2(m+ly+8)*i2(m+lz+8))
                        d2(2)=d2(2)+dens(prim,n)*
     $                      (di(m+lx  ,1)*di(m+ly  ,1)*i2(m+lz  )+
     $                       di(m+lx+1,1)*di(m+ly+1,1)*i2(m+lz+1)+
     $                       di(m+lx+2,1)*di(m+ly+2,1)*i2(m+lz+2)+
     $                       di(m+lx+3,1)*di(m+ly+3,1)*i2(m+lz+3)+
     $                       di(m+lx+4,1)*di(m+ly+4,1)*i2(m+lz+4)+
     $                       di(m+lx+5,1)*di(m+ly+5,1)*i2(m+lz+5)+
     $                       di(m+lx+6,1)*di(m+ly+6,1)*i2(m+lz+6)+
     $                       di(m+lx+7,1)*di(m+ly+7,1)*i2(m+lz+7)+
     $                       di(m+lx+8,1)*di(m+ly+8,1)*i2(m+lz+8))
                        d2(3)=d2(3)+dens(prim,n)*
     $                      (i2(m+lx  )*d2i(m+ly  ,1)*i2(m+lz  )+
     $                       i2(m+lx+1)*d2i(m+ly+1,1)*i2(m+lz+1)+
     $                       i2(m+lx+2)*d2i(m+ly+2,1)*i2(m+lz+2)+
     $                       i2(m+lx+3)*d2i(m+ly+3,1)*i2(m+lz+3)+
     $                       i2(m+lx+4)*d2i(m+ly+4,1)*i2(m+lz+4)+
     $                       i2(m+lx+5)*d2i(m+ly+5,1)*i2(m+lz+5)+
     $                       i2(m+lx+6)*d2i(m+ly+6,1)*i2(m+lz+6)+
     $                       i2(m+lx+7)*d2i(m+ly+7,1)*i2(m+lz+7)+
     $                       i2(m+lx+8)*d2i(m+ly+8,1)*i2(m+lz+8))
                        d2(4)=d2(4)+dens(prim,n)*
     $                      (di(m+lx  ,1)*i2(m+ly  )*di(m+lz  ,1)+
     $                       di(m+lx+1,1)*i2(m+ly+1)*di(m+lz+1,1)+
     $                       di(m+lx+2,1)*i2(m+ly+2)*di(m+lz+2,1)+
     $                       di(m+lx+3,1)*i2(m+ly+3)*di(m+lz+3,1)+
     $                       di(m+lx+4,1)*i2(m+ly+4)*di(m+lz+4,1)+
     $                       di(m+lx+5,1)*i2(m+ly+5)*di(m+lz+5,1)+
     $                       di(m+lx+6,1)*i2(m+ly+6)*di(m+lz+6,1)+
     $                       di(m+lx+7,1)*i2(m+ly+7)*di(m+lz+7,1)+
     $                       di(m+lx+8,1)*i2(m+ly+8)*di(m+lz+8,1))
                        d2(5)=d2(5)+dens(prim,n)*
     $                      (i2(m+lx  )*di(m+ly  ,1)*di(m+lz  ,1)+
     $                       i2(m+lx+1)*di(m+ly+1,1)*di(m+lz+1,1)+
     $                       i2(m+lx+2)*di(m+ly+2,1)*di(m+lz+2,1)+
     $                       i2(m+lx+3)*di(m+ly+3,1)*di(m+lz+3,1)+
     $                       i2(m+lx+4)*di(m+ly+4,1)*di(m+lz+4,1)+
     $                       i2(m+lx+5)*di(m+ly+5,1)*di(m+lz+5,1)+
     $                       i2(m+lx+6)*di(m+ly+6,1)*di(m+lz+6,1)+
     $                       i2(m+lx+7)*di(m+ly+7,1)*di(m+lz+7,1)+
     $                       i2(m+lx+8)*di(m+ly+8,1)*di(m+lz+8,1))
                        d2(6)=d2(6)+dens(prim,n)*
     $                      (i2(m+lx  )*i2(m+ly  )*d2i(m+lz  ,1)+
     $                       i2(m+lx+1)*i2(m+ly+1)*d2i(m+lz+1,1)+
     $                       i2(m+lx+2)*i2(m+ly+2)*d2i(m+lz+2,1)+
     $                       i2(m+lx+3)*i2(m+ly+3)*d2i(m+lz+3,1)+
     $                       i2(m+lx+4)*i2(m+ly+4)*d2i(m+lz+4,1)+
     $                       i2(m+lx+5)*i2(m+ly+5)*d2i(m+lz+5,1)+
     $                       i2(m+lx+6)*i2(m+ly+6)*d2i(m+lz+6,1)+
     $                       i2(m+lx+7)*i2(m+ly+7)*d2i(m+lz+7,1)+
     $                       i2(m+lx+8)*i2(m+ly+8)*d2i(m+lz+8,1))
                        m=m+9
 1300                continue
                     if (npass.lt.4) then
                        m=1
                        do 1301 prim=1,nv
                           der(1,2)=der(1,2)+dens(prim,n)*
     $                        (di(m+lx  ,2)*i2(m+ly  )*i2(m+lz  )+
     $                         di(m+lx+1,2)*i2(m+ly+1)*i2(m+lz+1)+
     $                         di(m+lx+2,2)*i2(m+ly+2)*i2(m+lz+2)+
     $                         di(m+lx+3,2)*i2(m+ly+3)*i2(m+lz+3)+
     $                         di(m+lx+4,2)*i2(m+ly+4)*i2(m+lz+4)+
     $                         di(m+lx+5,2)*i2(m+ly+5)*i2(m+lz+5)+
     $                         di(m+lx+6,2)*i2(m+ly+6)*i2(m+lz+6)+
     $                         di(m+lx+7,2)*i2(m+ly+7)*i2(m+lz+7)+
     $                         di(m+lx+8,2)*i2(m+ly+8)*i2(m+lz+8))
                           der(2,2)=der(2,2)+dens(prim,n)*
     $                        (i2(m+lx  )*di(m+ly  ,2)*i2(m+lz  )+
     $                         i2(m+lx+1)*di(m+ly+1,2)*i2(m+lz+1)+
     $                         i2(m+lx+2)*di(m+ly+2,2)*i2(m+lz+2)+
     $                         i2(m+lx+3)*di(m+ly+3,2)*i2(m+lz+3)+
     $                         i2(m+lx+4)*di(m+ly+4,2)*i2(m+lz+4)+
     $                         i2(m+lx+5)*di(m+ly+5,2)*i2(m+lz+5)+
     $                         i2(m+lx+6)*di(m+ly+6,2)*i2(m+lz+6)+
     $                         i2(m+lx+7)*di(m+ly+7,2)*i2(m+lz+7)+
     $                         i2(m+lx+8)*di(m+ly+8,2)*i2(m+lz+8))
                           der(3,2)=der(3,2)+dens(prim,n)*
     $                        (i2(m+lx  )*i2(m+ly  )*di(m+lz  ,2)+
     $                         i2(m+lx+1)*i2(m+ly+1)*di(m+lz+1,2)+
     $                         i2(m+lx+2)*i2(m+ly+2)*di(m+lz+2,2)+
     $                         i2(m+lx+3)*i2(m+ly+3)*di(m+lz+3,2)+
     $                         i2(m+lx+4)*i2(m+ly+4)*di(m+lz+4,2)+
     $                         i2(m+lx+5)*i2(m+ly+5)*di(m+lz+5,2)+
     $                         i2(m+lx+6)*i2(m+ly+6)*di(m+lz+6,2)+
     $                         i2(m+lx+7)*i2(m+ly+7)*di(m+lz+7,2)+
     $                         i2(m+lx+8)*i2(m+ly+8)*di(m+lz+8,2))
                           d2(7)=d2(7)+dens(prim,n)*
     $                         (d2i(m+lx  ,2)*i2(m+ly  )*i2(m+lz  )+
     $                          d2i(m+lx+1,2)*i2(m+ly+1)*i2(m+lz+1)+
     $                          d2i(m+lx+2,2)*i2(m+ly+2)*i2(m+lz+2)+
     $                          d2i(m+lx+3,2)*i2(m+ly+3)*i2(m+lz+3)+
     $                          d2i(m+lx+4,2)*i2(m+ly+4)*i2(m+lz+4)+
     $                          d2i(m+lx+5,2)*i2(m+ly+5)*i2(m+lz+5)+
     $                          d2i(m+lx+6,2)*i2(m+ly+6)*i2(m+lz+6)+
     $                          d2i(m+lx+7,2)*i2(m+ly+7)*i2(m+lz+7)+
     $                          d2i(m+lx+8,2)*i2(m+ly+8)*i2(m+lz+8))
                           d2(8)=d2(8)+dens(prim,n)*
     $                         (di(m+lx  ,2)*di(m+ly  ,1)*i2(m+lz  )+
     $                          di(m+lx+1,2)*di(m+ly+1,1)*i2(m+lz+1)+
     $                          di(m+lx+2,2)*di(m+ly+2,1)*i2(m+lz+2)+
     $                          di(m+lx+3,2)*di(m+ly+3,1)*i2(m+lz+3)+
     $                          di(m+lx+4,2)*di(m+ly+4,1)*i2(m+lz+4)+
     $                          di(m+lx+5,2)*di(m+ly+5,1)*i2(m+lz+5)+
     $                          di(m+lx+6,2)*di(m+ly+6,1)*i2(m+lz+6)+
     $                          di(m+lx+7,2)*di(m+ly+7,1)*i2(m+lz+7)+
     $                          di(m+lx+8,2)*di(m+ly+8,1)*i2(m+lz+8))
                           d2(9)=d2(9)+dens(prim,n)*
     $                         (di(m+lx  ,2)*i2(m+ly  )*di(m+lz  ,1)+
     $                          di(m+lx+1,2)*i2(m+ly+1)*di(m+lz+1,1)+
     $                          di(m+lx+2,2)*i2(m+ly+2)*di(m+lz+2,1)+
     $                          di(m+lx+3,2)*i2(m+ly+3)*di(m+lz+3,1)+
     $                          di(m+lx+4,2)*i2(m+ly+4)*di(m+lz+4,1)+
     $                          di(m+lx+5,2)*i2(m+ly+5)*di(m+lz+5,1)+
     $                          di(m+lx+6,2)*i2(m+ly+6)*di(m+lz+6,1)+
     $                          di(m+lx+7,2)*i2(m+ly+7)*di(m+lz+7,1)+
     $                          di(m+lx+8,2)*i2(m+ly+8)*di(m+lz+8,1))
                           d2(10)=d2(10)+dens(prim,n)*
     $                         (d2i(m+lx  ,3)*i2(m+ly  )*i2(m+lz  )+
     $                          d2i(m+lx+1,3)*i2(m+ly+1)*i2(m+lz+1)+
     $                          d2i(m+lx+2,3)*i2(m+ly+2)*i2(m+lz+2)+
     $                          d2i(m+lx+3,3)*i2(m+ly+3)*i2(m+lz+3)+
     $                          d2i(m+lx+4,3)*i2(m+ly+4)*i2(m+lz+4)+
     $                          d2i(m+lx+5,3)*i2(m+ly+5)*i2(m+lz+5)+
     $                          d2i(m+lx+6,3)*i2(m+ly+6)*i2(m+lz+6)+
     $                          d2i(m+lx+7,3)*i2(m+ly+7)*i2(m+lz+7)+
     $                          d2i(m+lx+8,3)*i2(m+ly+8)*i2(m+lz+8))
                           d2(11)=d2(11)+dens(prim,n)*
     $                         (di(m+lx  ,1)*di(m+ly  ,2)*i2(m+lz  )+
     $                          di(m+lx+1,1)*di(m+ly+1,2)*i2(m+lz+1)+
     $                          di(m+lx+2,1)*di(m+ly+2,2)*i2(m+lz+2)+
     $                          di(m+lx+3,1)*di(m+ly+3,2)*i2(m+lz+3)+
     $                          di(m+lx+4,1)*di(m+ly+4,2)*i2(m+lz+4)+
     $                          di(m+lx+5,1)*di(m+ly+5,2)*i2(m+lz+5)+
     $                          di(m+lx+6,1)*di(m+ly+6,2)*i2(m+lz+6)+
     $                          di(m+lx+7,1)*di(m+ly+7,2)*i2(m+lz+7)+
     $                          di(m+lx+8,1)*di(m+ly+8,2)*i2(m+lz+8))
                           d2(12)=d2(12)+dens(prim,n)*
     $                         (i2(m+lx  )*d2i(m+ly  ,2)*i2(m+lz  )+
     $                          i2(m+lx+1)*d2i(m+ly+1,2)*i2(m+lz+1)+
     $                          i2(m+lx+2)*d2i(m+ly+2,2)*i2(m+lz+2)+
     $                          i2(m+lx+3)*d2i(m+ly+3,2)*i2(m+lz+3)+
     $                          i2(m+lx+4)*d2i(m+ly+4,2)*i2(m+lz+4)+
     $                          i2(m+lx+5)*d2i(m+ly+5,2)*i2(m+lz+5)+
     $                          i2(m+lx+6)*d2i(m+ly+6,2)*i2(m+lz+6)+
     $                          i2(m+lx+7)*d2i(m+ly+7,2)*i2(m+lz+7)+
     $                          i2(m+lx+8)*d2i(m+ly+8,2)*i2(m+lz+8))
                           d2(13)=d2(13)+dens(prim,n)*
     $                         (i2(m+lx  )*di(m+ly  ,2)*di(m+lz  ,1)+
     $                          i2(m+lx+1)*di(m+ly+1,2)*di(m+lz+1,1)+
     $                          i2(m+lx+2)*di(m+ly+2,2)*di(m+lz+2,1)+
     $                          i2(m+lx+3)*di(m+ly+3,2)*di(m+lz+3,1)+
     $                          i2(m+lx+4)*di(m+ly+4,2)*di(m+lz+4,1)+
     $                          i2(m+lx+5)*di(m+ly+5,2)*di(m+lz+5,1)+
     $                          i2(m+lx+6)*di(m+ly+6,2)*di(m+lz+6,1)+
     $                          i2(m+lx+7)*di(m+ly+7,2)*di(m+lz+7,1)+
     $                          i2(m+lx+8)*di(m+ly+8,2)*di(m+lz+8,1))
                           d2(14)=d2(14)+dens(prim,n)*
     $                         (di(m+lx  ,2)*di(m+ly  ,2)*i2(m+lz  )+
     $                          di(m+lx+1,2)*di(m+ly+1,2)*i2(m+lz+1)+
     $                          di(m+lx+2,2)*di(m+ly+2,2)*i2(m+lz+2)+
     $                          di(m+lx+3,2)*di(m+ly+3,2)*i2(m+lz+3)+
     $                          di(m+lx+4,2)*di(m+ly+4,2)*i2(m+lz+4)+
     $                          di(m+lx+5,2)*di(m+ly+5,2)*i2(m+lz+5)+
     $                          di(m+lx+6,2)*di(m+ly+6,2)*i2(m+lz+6)+
     $                          di(m+lx+7,2)*di(m+ly+7,2)*i2(m+lz+7)+
     $                          di(m+lx+8,2)*di(m+ly+8,2)*i2(m+lz+8))
                           d2(15)=d2(15)+dens(prim,n)*
     $                         (i2(m+lx  )*d2i(m+ly  ,3)*i2(m+lz  )+
     $                          i2(m+lx+1)*d2i(m+ly+1,3)*i2(m+lz+1)+
     $                          i2(m+lx+2)*d2i(m+ly+2,3)*i2(m+lz+2)+
     $                          i2(m+lx+3)*d2i(m+ly+3,3)*i2(m+lz+3)+
     $                          i2(m+lx+4)*d2i(m+ly+4,3)*i2(m+lz+4)+
     $                          i2(m+lx+5)*d2i(m+ly+5,3)*i2(m+lz+5)+
     $                          i2(m+lx+6)*d2i(m+ly+6,3)*i2(m+lz+6)+
     $                          i2(m+lx+7)*d2i(m+ly+7,3)*i2(m+lz+7)+
     $                          i2(m+lx+8)*d2i(m+ly+8,3)*i2(m+lz+8))
                           d2(16)=d2(16)+dens(prim,n)*
     $                         (di(m+lx  ,1)*i2(m+ly  )*di(m+lz  ,2)+
     $                          di(m+lx+1,1)*i2(m+ly+1)*di(m+lz+1,2)+
     $                          di(m+lx+2,1)*i2(m+ly+2)*di(m+lz+2,2)+
     $                          di(m+lx+3,1)*i2(m+ly+3)*di(m+lz+3,2)+
     $                          di(m+lx+4,1)*i2(m+ly+4)*di(m+lz+4,2)+
     $                          di(m+lx+5,1)*i2(m+ly+5)*di(m+lz+5,2)+
     $                          di(m+lx+6,1)*i2(m+ly+6)*di(m+lz+6,2)+
     $                          di(m+lx+7,1)*i2(m+ly+7)*di(m+lz+7,2)+
     $                          di(m+lx+8,1)*i2(m+ly+8)*di(m+lz+8,2))
                           d2(17)=d2(17)+dens(prim,n)*
     $                         (i2(m+lx  )*di(m+ly  ,1)*di(m+lz  ,2)+
     $                          i2(m+lx+1)*di(m+ly+1,1)*di(m+lz+1,2)+
     $                          i2(m+lx+2)*di(m+ly+2,1)*di(m+lz+2,2)+
     $                          i2(m+lx+3)*di(m+ly+3,1)*di(m+lz+3,2)+
     $                          i2(m+lx+4)*di(m+ly+4,1)*di(m+lz+4,2)+
     $                          i2(m+lx+5)*di(m+ly+5,1)*di(m+lz+5,2)+
     $                          i2(m+lx+6)*di(m+ly+6,1)*di(m+lz+6,2)+
     $                          i2(m+lx+7)*di(m+ly+7,1)*di(m+lz+7,2)+
     $                          i2(m+lx+8)*di(m+ly+8,1)*di(m+lz+8,2))
                           d2(18)=d2(18)+dens(prim,n)*
     $                         (i2(m+lx  )*i2(m+ly  )*d2i(m+lz  ,2)+
     $                          i2(m+lx+1)*i2(m+ly+1)*d2i(m+lz+1,2)+
     $                          i2(m+lx+2)*i2(m+ly+2)*d2i(m+lz+2,2)+
     $                          i2(m+lx+3)*i2(m+ly+3)*d2i(m+lz+3,2)+
     $                          i2(m+lx+4)*i2(m+ly+4)*d2i(m+lz+4,2)+
     $                          i2(m+lx+5)*i2(m+ly+5)*d2i(m+lz+5,2)+
     $                          i2(m+lx+6)*i2(m+ly+6)*d2i(m+lz+6,2)+
     $                          i2(m+lx+7)*i2(m+ly+7)*d2i(m+lz+7,2)+
     $                          i2(m+lx+8)*i2(m+ly+8)*d2i(m+lz+8,2))
                           d2(19)=d2(19)+dens(prim,n)*
     $                         (di(m+lx  ,2)*i2(m+ly  )*di(m+lz  ,2)+
     $                          di(m+lx+1,2)*i2(m+ly+1)*di(m+lz+1,2)+
     $                          di(m+lx+2,2)*i2(m+ly+2)*di(m+lz+2,2)+
     $                          di(m+lx+3,2)*i2(m+ly+3)*di(m+lz+3,2)+
     $                          di(m+lx+4,2)*i2(m+ly+4)*di(m+lz+4,2)+
     $                          di(m+lx+5,2)*i2(m+ly+5)*di(m+lz+5,2)+
     $                          di(m+lx+6,2)*i2(m+ly+6)*di(m+lz+6,2)+
     $                          di(m+lx+7,2)*i2(m+ly+7)*di(m+lz+7,2)+
     $                          di(m+lx+8,2)*i2(m+ly+8)*di(m+lz+8,2))
                           d2(20)=d2(20)+dens(prim,n)*
     $                         (i2(m+lx  )*di(m+ly  ,2)*di(m+lz  ,2)+
     $                          i2(m+lx+1)*di(m+ly+1,2)*di(m+lz+1,2)+
     $                          i2(m+lx+2)*di(m+ly+2,2)*di(m+lz+2,2)+
     $                          i2(m+lx+3)*di(m+ly+3,2)*di(m+lz+3,2)+
     $                          i2(m+lx+4)*di(m+ly+4,2)*di(m+lz+4,2)+
     $                          i2(m+lx+5)*di(m+ly+5,2)*di(m+lz+5,2)+
     $                          i2(m+lx+6)*di(m+ly+6,2)*di(m+lz+6,2)+
     $                          i2(m+lx+7)*di(m+ly+7,2)*di(m+lz+7,2)+
     $                          i2(m+lx+8)*di(m+ly+8,2)*di(m+lz+8,2))
                           d2(21)=d2(21)+dens(prim,n)*
     $                         (i2(m+lx  )*i2(m+ly  )*d2i(m+lz  ,3)+
     $                          i2(m+lx+1)*i2(m+ly+1)*d2i(m+lz+1,3)+
     $                          i2(m+lx+2)*i2(m+ly+2)*d2i(m+lz+2,3)+
     $                          i2(m+lx+3)*i2(m+ly+3)*d2i(m+lz+3,3)+
     $                          i2(m+lx+4)*i2(m+ly+4)*d2i(m+lz+4,3)+
     $                          i2(m+lx+5)*i2(m+ly+5)*d2i(m+lz+5,3)+
     $                          i2(m+lx+6)*i2(m+ly+6)*d2i(m+lz+6,3)+
     $                          i2(m+lx+7)*i2(m+ly+7)*d2i(m+lz+7,3)+
     $                          i2(m+lx+8)*i2(m+ly+8)*d2i(m+lz+8,3))
                           m=m+9
 1301                   continue
                     end if
                     if (npass.eq.1) then
                        m=1
                        do 1302 prim=1,nv
                           der(1,3)=der(1,3)+dens(prim,n)*
     $                        (di(m+lx  ,3)*i2(m+ly  )*i2(m+lz  )+
     $                         di(m+lx+1,3)*i2(m+ly+1)*i2(m+lz+1)+
     $                         di(m+lx+2,3)*i2(m+ly+2)*i2(m+lz+2)+
     $                         di(m+lx+3,3)*i2(m+ly+3)*i2(m+lz+3)+
     $                         di(m+lx+4,3)*i2(m+ly+4)*i2(m+lz+4)+
     $                         di(m+lx+5,3)*i2(m+ly+5)*i2(m+lz+5)+
     $                         di(m+lx+6,3)*i2(m+ly+6)*i2(m+lz+6)+
     $                         di(m+lx+7,3)*i2(m+ly+7)*i2(m+lz+7)+
     $                         di(m+lx+8,3)*i2(m+ly+8)*i2(m+lz+8))
                           der(2,3)=der(2,3)+dens(prim,n)*
     $                        (i2(m+lx  )*di(m+ly  ,3)*i2(m+lz  )+
     $                         i2(m+lx+1)*di(m+ly+1,3)*i2(m+lz+1)+
     $                         i2(m+lx+2)*di(m+ly+2,3)*i2(m+lz+2)+
     $                         i2(m+lx+3)*di(m+ly+3,3)*i2(m+lz+3)+
     $                         i2(m+lx+4)*di(m+ly+4,3)*i2(m+lz+4)+
     $                         i2(m+lx+5)*di(m+ly+5,3)*i2(m+lz+5)+
     $                         i2(m+lx+6)*di(m+ly+6,3)*i2(m+lz+6)+
     $                         i2(m+lx+7)*di(m+ly+7,3)*i2(m+lz+7)+
     $                         i2(m+lx+8)*di(m+ly+8,3)*i2(m+lz+8))
                           der(3,3)=der(3,3)+dens(prim,n)*
     $                        (i2(m+lx  )*i2(m+ly  )*di(m+lz  ,3)+
     $                         i2(m+lx+1)*i2(m+ly+1)*di(m+lz+1,3)+
     $                         i2(m+lx+2)*i2(m+ly+2)*di(m+lz+2,3)+
     $                         i2(m+lx+3)*i2(m+ly+3)*di(m+lz+3,3)+
     $                         i2(m+lx+4)*i2(m+ly+4)*di(m+lz+4,3)+
     $                         i2(m+lx+5)*i2(m+ly+5)*di(m+lz+5,3)+
     $                         i2(m+lx+6)*i2(m+ly+6)*di(m+lz+6,3)+
     $                         i2(m+lx+7)*i2(m+ly+7)*di(m+lz+7,3)+
     $                         i2(m+lx+8)*i2(m+ly+8)*di(m+lz+8,3))
                           d2(22)=d2(22)+dens(prim,n)*
     $                         (d2i(m+lx  ,4)*i2(m+ly  )*i2(m+lz  )+
     $                          d2i(m+lx+1,4)*i2(m+ly+1)*i2(m+lz+1)+
     $                          d2i(m+lx+2,4)*i2(m+ly+2)*i2(m+lz+2)+
     $                          d2i(m+lx+3,4)*i2(m+ly+3)*i2(m+lz+3)+
     $                          d2i(m+lx+4,4)*i2(m+ly+4)*i2(m+lz+4)+
     $                          d2i(m+lx+5,4)*i2(m+ly+5)*i2(m+lz+5)+
     $                          d2i(m+lx+6,4)*i2(m+ly+6)*i2(m+lz+6)+
     $                          d2i(m+lx+7,4)*i2(m+ly+7)*i2(m+lz+7)+
     $                          d2i(m+lx+8,4)*i2(m+ly+8)*i2(m+lz+8))
                           d2(23)=d2(23)+dens(prim,n)*
     $                         (di(m+lx  ,3)*di(m+ly  ,1)*i2(m+lz  )+
     $                          di(m+lx+1,3)*di(m+ly+1,1)*i2(m+lz+1)+
     $                          di(m+lx+2,3)*di(m+ly+2,1)*i2(m+lz+2)+
     $                          di(m+lx+3,3)*di(m+ly+3,1)*i2(m+lz+3)+
     $                          di(m+lx+4,3)*di(m+ly+4,1)*i2(m+lz+4)+
     $                          di(m+lx+5,3)*di(m+ly+5,1)*i2(m+lz+5)+
     $                          di(m+lx+6,3)*di(m+ly+6,1)*i2(m+lz+6)+
     $                          di(m+lx+7,3)*di(m+ly+7,1)*i2(m+lz+7)+
     $                          di(m+lx+8,3)*di(m+ly+8,1)*i2(m+lz+8))
                           d2(24)=d2(24)+dens(prim,n)*
     $                         (di(m+lx  ,3)*i2(m+ly  )*di(m+lz  ,1)+
     $                          di(m+lx+1,3)*i2(m+ly+1)*di(m+lz+1,1)+
     $                          di(m+lx+2,3)*i2(m+ly+2)*di(m+lz+2,1)+
     $                          di(m+lx+3,3)*i2(m+ly+3)*di(m+lz+3,1)+
     $                          di(m+lx+4,3)*i2(m+ly+4)*di(m+lz+4,1)+
     $                          di(m+lx+5,3)*i2(m+ly+5)*di(m+lz+5,1)+
     $                          di(m+lx+6,3)*i2(m+ly+6)*di(m+lz+6,1)+
     $                          di(m+lx+7,3)*i2(m+ly+7)*di(m+lz+7,1)+
     $                          di(m+lx+8,3)*i2(m+ly+8)*di(m+lz+8,1))
                           d2(25)=d2(25)+dens(prim,n)*
     $                         (d2i(m+lx  ,5)*i2(m+ly  )*i2(m+lz  )+
     $                          d2i(m+lx+1,5)*i2(m+ly+1)*i2(m+lz+1)+
     $                          d2i(m+lx+2,5)*i2(m+ly+2)*i2(m+lz+2)+
     $                          d2i(m+lx+3,5)*i2(m+ly+3)*i2(m+lz+3)+
     $                          d2i(m+lx+4,5)*i2(m+ly+4)*i2(m+lz+4)+
     $                          d2i(m+lx+5,5)*i2(m+ly+5)*i2(m+lz+5)+
     $                          d2i(m+lx+6,5)*i2(m+ly+6)*i2(m+lz+6)+
     $                          d2i(m+lx+7,5)*i2(m+ly+7)*i2(m+lz+7)+
     $                          d2i(m+lx+8,5)*i2(m+ly+8)*i2(m+lz+8))
                           d2(26)=d2(26)+dens(prim,n)*
     $                         (di(m+lx  ,3)*di(m+ly  ,2)*i2(m+lz  )+
     $                          di(m+lx+1,3)*di(m+ly+1,2)*i2(m+lz+1)+
     $                          di(m+lx+2,3)*di(m+ly+2,2)*i2(m+lz+2)+
     $                          di(m+lx+3,3)*di(m+ly+3,2)*i2(m+lz+3)+
     $                          di(m+lx+4,3)*di(m+ly+4,2)*i2(m+lz+4)+
     $                          di(m+lx+5,3)*di(m+ly+5,2)*i2(m+lz+5)+
     $                          di(m+lx+6,3)*di(m+ly+6,2)*i2(m+lz+6)+
     $                          di(m+lx+7,3)*di(m+ly+7,2)*i2(m+lz+7)+
     $                          di(m+lx+8,3)*di(m+ly+8,2)*i2(m+lz+8))
                           d2(27)=d2(27)+dens(prim,n)*
     $                         (di(m+lx  ,3)*i2(m+ly  )*di(m+lz  ,2)+
     $                          di(m+lx+1,3)*i2(m+ly+1)*di(m+lz+1,2)+
     $                          di(m+lx+2,3)*i2(m+ly+2)*di(m+lz+2,2)+
     $                          di(m+lx+3,3)*i2(m+ly+3)*di(m+lz+3,2)+
     $                          di(m+lx+4,3)*i2(m+ly+4)*di(m+lz+4,2)+
     $                          di(m+lx+5,3)*i2(m+ly+5)*di(m+lz+5,2)+
     $                          di(m+lx+6,3)*i2(m+ly+6)*di(m+lz+6,2)+
     $                          di(m+lx+7,3)*i2(m+ly+7)*di(m+lz+7,2)+
     $                          di(m+lx+8,3)*i2(m+ly+8)*di(m+lz+8,2))
                           d2(28)=d2(28)+dens(prim,n)*
     $                         (d2i(m+lx  ,6)*i2(m+ly  )*i2(m+lz  )+
     $                          d2i(m+lx+1,6)*i2(m+ly+1)*i2(m+lz+1)+
     $                          d2i(m+lx+2,6)*i2(m+ly+2)*i2(m+lz+2)+
     $                          d2i(m+lx+3,6)*i2(m+ly+3)*i2(m+lz+3)+
     $                          d2i(m+lx+4,6)*i2(m+ly+4)*i2(m+lz+4)+
     $                          d2i(m+lx+5,6)*i2(m+ly+5)*i2(m+lz+5)+
     $                          d2i(m+lx+6,6)*i2(m+ly+6)*i2(m+lz+6)+
     $                          d2i(m+lx+7,6)*i2(m+ly+7)*i2(m+lz+7)+
     $                          d2i(m+lx+8,6)*i2(m+ly+8)*i2(m+lz+8))
                           d2(29)=d2(29)+dens(prim,n)*
     $                         (di(m+lx  ,1)*di(m+ly  ,3)*i2(m+lz  )+
     $                          di(m+lx+1,1)*di(m+ly+1,3)*i2(m+lz+1)+
     $                          di(m+lx+2,1)*di(m+ly+2,3)*i2(m+lz+2)+
     $                          di(m+lx+3,1)*di(m+ly+3,3)*i2(m+lz+3)+
     $                          di(m+lx+4,1)*di(m+ly+4,3)*i2(m+lz+4)+
     $                          di(m+lx+5,1)*di(m+ly+5,3)*i2(m+lz+5)+
     $                          di(m+lx+6,1)*di(m+ly+6,3)*i2(m+lz+6)+
     $                          di(m+lx+7,1)*di(m+ly+7,3)*i2(m+lz+7)+
     $                          di(m+lx+8,1)*di(m+ly+8,3)*i2(m+lz+8))
                           d2(30)=d2(30)+dens(prim,n)*
     $                         (i2(m+lx  )*d2i(m+ly  ,4)*i2(m+lz  )+
     $                          i2(m+lx+1)*d2i(m+ly+1,4)*i2(m+lz+1)+
     $                          i2(m+lx+2)*d2i(m+ly+2,4)*i2(m+lz+2)+
     $                          i2(m+lx+3)*d2i(m+ly+3,4)*i2(m+lz+3)+
     $                          i2(m+lx+4)*d2i(m+ly+4,4)*i2(m+lz+4)+
     $                          i2(m+lx+5)*d2i(m+ly+5,4)*i2(m+lz+5)+
     $                          i2(m+lx+6)*d2i(m+ly+6,4)*i2(m+lz+6)+
     $                          i2(m+lx+7)*d2i(m+ly+7,4)*i2(m+lz+7)+
     $                          i2(m+lx+8)*d2i(m+ly+8,4)*i2(m+lz+8))
                           d2(31)=d2(31)+dens(prim,n)*
     $                         (i2(m+lx  )*di(m+ly  ,3)*di(m+lz  ,1)+
     $                          i2(m+lx+1)*di(m+ly+1,3)*di(m+lz+1,1)+
     $                          i2(m+lx+2)*di(m+ly+2,3)*di(m+lz+2,1)+
     $                          i2(m+lx+3)*di(m+ly+3,3)*di(m+lz+3,1)+
     $                          i2(m+lx+4)*di(m+ly+4,3)*di(m+lz+4,1)+
     $                          i2(m+lx+5)*di(m+ly+5,3)*di(m+lz+5,1)+
     $                          i2(m+lx+6)*di(m+ly+6,3)*di(m+lz+6,1)+
     $                          i2(m+lx+7)*di(m+ly+7,3)*di(m+lz+7,1)+
     $                          i2(m+lx+8)*di(m+ly+8,3)*di(m+lz+8,1))
                           d2(32)=d2(32)+dens(prim,n)*
     $                         (di(m+lx  ,2)*di(m+ly  ,3)*i2(m+lz  )+
     $                          di(m+lx+1,2)*di(m+ly+1,3)*i2(m+lz+1)+
     $                          di(m+lx+2,2)*di(m+ly+2,3)*i2(m+lz+2)+
     $                          di(m+lx+3,2)*di(m+ly+3,3)*i2(m+lz+3)+
     $                          di(m+lx+4,2)*di(m+ly+4,3)*i2(m+lz+4)+
     $                          di(m+lx+5,2)*di(m+ly+5,3)*i2(m+lz+5)+
     $                          di(m+lx+6,2)*di(m+ly+6,3)*i2(m+lz+6)+
     $                          di(m+lx+7,2)*di(m+ly+7,3)*i2(m+lz+7)+
     $                          di(m+lx+8,2)*di(m+ly+8,3)*i2(m+lz+8))
                           d2(33)=d2(33)+dens(prim,n)*
     $                         (i2(m+lx  )*d2i(m+ly  ,5)*i2(m+lz  )+
     $                          i2(m+lx+1)*d2i(m+ly+1,5)*i2(m+lz+1)+
     $                          i2(m+lx+2)*d2i(m+ly+2,5)*i2(m+lz+2)+
     $                          i2(m+lx+3)*d2i(m+ly+3,5)*i2(m+lz+3)+
     $                          i2(m+lx+4)*d2i(m+ly+4,5)*i2(m+lz+4)+
     $                          i2(m+lx+5)*d2i(m+ly+5,5)*i2(m+lz+5)+
     $                          i2(m+lx+6)*d2i(m+ly+6,5)*i2(m+lz+6)+
     $                          i2(m+lx+7)*d2i(m+ly+7,5)*i2(m+lz+7)+
     $                          i2(m+lx+8)*d2i(m+ly+8,5)*i2(m+lz+8))
                           d2(34)=d2(34)+dens(prim,n)*
     $                         (i2(m+lx  )*di(m+ly  ,3)*di(m+lz  ,2)+
     $                          i2(m+lx+1)*di(m+ly+1,3)*di(m+lz+1,2)+
     $                          i2(m+lx+2)*di(m+ly+2,3)*di(m+lz+2,2)+
     $                          i2(m+lx+3)*di(m+ly+3,3)*di(m+lz+3,2)+
     $                          i2(m+lx+4)*di(m+ly+4,3)*di(m+lz+4,2)+
     $                          i2(m+lx+5)*di(m+ly+5,3)*di(m+lz+5,2)+
     $                          i2(m+lx+6)*di(m+ly+6,3)*di(m+lz+6,2)+
     $                          i2(m+lx+7)*di(m+ly+7,3)*di(m+lz+7,2)+
     $                          i2(m+lx+8)*di(m+ly+8,3)*di(m+lz+8,2))
                           d2(35)=d2(35)+dens(prim,n)*
     $                         (di(m+lx  ,3)*di(m+ly  ,3)*i2(m+lz  )+
     $                          di(m+lx+1,3)*di(m+ly+1,3)*i2(m+lz+1)+
     $                          di(m+lx+2,3)*di(m+ly+2,3)*i2(m+lz+2)+
     $                          di(m+lx+3,3)*di(m+ly+3,3)*i2(m+lz+3)+
     $                          di(m+lx+4,3)*di(m+ly+4,3)*i2(m+lz+4)+
     $                          di(m+lx+5,3)*di(m+ly+5,3)*i2(m+lz+5)+
     $                          di(m+lx+6,3)*di(m+ly+6,3)*i2(m+lz+6)+
     $                          di(m+lx+7,3)*di(m+ly+7,3)*i2(m+lz+7)+
     $                          di(m+lx+8,3)*di(m+ly+8,3)*i2(m+lz+8))
                           d2(36)=d2(36)+dens(prim,n)*
     $                         (i2(m+lx  )*d2i(m+ly  ,6)*i2(m+lz  )+
     $                          i2(m+lx+1)*d2i(m+ly+1,6)*i2(m+lz+1)+
     $                          i2(m+lx+2)*d2i(m+ly+2,6)*i2(m+lz+2)+
     $                          i2(m+lx+3)*d2i(m+ly+3,6)*i2(m+lz+3)+
     $                          i2(m+lx+4)*d2i(m+ly+4,6)*i2(m+lz+4)+
     $                          i2(m+lx+5)*d2i(m+ly+5,6)*i2(m+lz+5)+
     $                          i2(m+lx+6)*d2i(m+ly+6,6)*i2(m+lz+6)+
     $                          i2(m+lx+7)*d2i(m+ly+7,6)*i2(m+lz+7)+
     $                          i2(m+lx+8)*d2i(m+ly+8,6)*i2(m+lz+8))
                           d2(37)=d2(37)+dens(prim,n)*
     $                         (di(m+lx  ,1)*i2(m+ly  )*di(m+lz  ,3)+
     $                          di(m+lx+1,1)*i2(m+ly+1)*di(m+lz+1,3)+
     $                          di(m+lx+2,1)*i2(m+ly+2)*di(m+lz+2,3)+
     $                          di(m+lx+3,1)*i2(m+ly+3)*di(m+lz+3,3)+
     $                          di(m+lx+4,1)*i2(m+ly+4)*di(m+lz+4,3)+
     $                          di(m+lx+5,1)*i2(m+ly+5)*di(m+lz+5,3)+
     $                          di(m+lx+6,1)*i2(m+ly+6)*di(m+lz+6,3)+
     $                          di(m+lx+7,1)*i2(m+ly+7)*di(m+lz+7,3)+
     $                          di(m+lx+8,1)*i2(m+ly+8)*di(m+lz+8,3))
                           d2(38)=d2(38)+dens(prim,n)*
     $                         (i2(m+lx  )*di(m+ly  ,1)*di(m+lz  ,3)+
     $                          i2(m+lx+1)*di(m+ly+1,1)*di(m+lz+1,3)+
     $                          i2(m+lx+2)*di(m+ly+2,1)*di(m+lz+2,3)+
     $                          i2(m+lx+3)*di(m+ly+3,1)*di(m+lz+3,3)+
     $                          i2(m+lx+4)*di(m+ly+4,1)*di(m+lz+4,3)+
     $                          i2(m+lx+5)*di(m+ly+5,1)*di(m+lz+5,3)+
     $                          i2(m+lx+6)*di(m+ly+6,1)*di(m+lz+6,3)+
     $                          i2(m+lx+7)*di(m+ly+7,1)*di(m+lz+7,3)+
     $                          i2(m+lx+8)*di(m+ly+8,1)*di(m+lz+8,3))
                           d2(39)=d2(39)+dens(prim,n)*
     $                         (i2(m+lx  )*i2(m+ly  )*d2i(m+lz  ,4)+
     $                          i2(m+lx+1)*i2(m+ly+1)*d2i(m+lz+1,4)+
     $                          i2(m+lx+2)*i2(m+ly+2)*d2i(m+lz+2,4)+
     $                          i2(m+lx+3)*i2(m+ly+3)*d2i(m+lz+3,4)+
     $                          i2(m+lx+4)*i2(m+ly+4)*d2i(m+lz+4,4)+
     $                          i2(m+lx+5)*i2(m+ly+5)*d2i(m+lz+5,4)+
     $                          i2(m+lx+6)*i2(m+ly+6)*d2i(m+lz+6,4)+
     $                          i2(m+lx+7)*i2(m+ly+7)*d2i(m+lz+7,4)+
     $                          i2(m+lx+8)*i2(m+ly+8)*d2i(m+lz+8,4))
                           d2(40)=d2(40)+dens(prim,n)*
     $                         (di(m+lx  ,2)*i2(m+ly  )*di(m+lz  ,3)+
     $                          di(m+lx+1,2)*i2(m+ly+1)*di(m+lz+1,3)+
     $                          di(m+lx+2,2)*i2(m+ly+2)*di(m+lz+2,3)+
     $                          di(m+lx+3,2)*i2(m+ly+3)*di(m+lz+3,3)+
     $                          di(m+lx+4,2)*i2(m+ly+4)*di(m+lz+4,3)+
     $                          di(m+lx+5,2)*i2(m+ly+5)*di(m+lz+5,3)+
     $                          di(m+lx+6,2)*i2(m+ly+6)*di(m+lz+6,3)+
     $                          di(m+lx+7,2)*i2(m+ly+7)*di(m+lz+7,3)+
     $                          di(m+lx+8,2)*i2(m+ly+8)*di(m+lz+8,3))
                           d2(41)=d2(41)+dens(prim,n)*
     $                         (i2(m+lx  )*di(m+ly  ,2)*di(m+lz  ,3)+
     $                          i2(m+lx+1)*di(m+ly+1,2)*di(m+lz+1,3)+
     $                          i2(m+lx+2)*di(m+ly+2,2)*di(m+lz+2,3)+
     $                          i2(m+lx+3)*di(m+ly+3,2)*di(m+lz+3,3)+
     $                          i2(m+lx+4)*di(m+ly+4,2)*di(m+lz+4,3)+
     $                          i2(m+lx+5)*di(m+ly+5,2)*di(m+lz+5,3)+
     $                          i2(m+lx+6)*di(m+ly+6,2)*di(m+lz+6,3)+
     $                          i2(m+lx+7)*di(m+ly+7,2)*di(m+lz+7,3)+
     $                          i2(m+lx+8)*di(m+ly+8,2)*di(m+lz+8,3))
                           d2(42)=d2(42)+dens(prim,n)*
     $                         (i2(m+lx  )*i2(m+ly  )*d2i(m+lz  ,5)+
     $                          i2(m+lx+1)*i2(m+ly+1)*d2i(m+lz+1,5)+
     $                          i2(m+lx+2)*i2(m+ly+2)*d2i(m+lz+2,5)+
     $                          i2(m+lx+3)*i2(m+ly+3)*d2i(m+lz+3,5)+
     $                          i2(m+lx+4)*i2(m+ly+4)*d2i(m+lz+4,5)+
     $                          i2(m+lx+5)*i2(m+ly+5)*d2i(m+lz+5,5)+
     $                          i2(m+lx+6)*i2(m+ly+6)*d2i(m+lz+6,5)+
     $                          i2(m+lx+7)*i2(m+ly+7)*d2i(m+lz+7,5)+
     $                          i2(m+lx+8)*i2(m+ly+8)*d2i(m+lz+8,5))
                           d2(43)=d2(43)+dens(prim,n)*
     $                         (di(m+lx  ,3)*i2(m+ly  )*di(m+lz  ,3)+
     $                          di(m+lx+1,3)*i2(m+ly+1)*di(m+lz+1,3)+
     $                          di(m+lx+2,3)*i2(m+ly+2)*di(m+lz+2,3)+
     $                          di(m+lx+3,3)*i2(m+ly+3)*di(m+lz+3,3)+
     $                          di(m+lx+4,3)*i2(m+ly+4)*di(m+lz+4,3)+
     $                          di(m+lx+5,3)*i2(m+ly+5)*di(m+lz+5,3)+
     $                          di(m+lx+6,3)*i2(m+ly+6)*di(m+lz+6,3)+
     $                          di(m+lx+7,3)*i2(m+ly+7)*di(m+lz+7,3)+
     $                          di(m+lx+8,3)*i2(m+ly+8)*di(m+lz+8,3))
                           d2(44)=d2(44)+dens(prim,n)*
     $                         (i2(m+lx  )*di(m+ly  ,3)*di(m+lz  ,3)+
     $                          i2(m+lx+1)*di(m+ly+1,3)*di(m+lz+1,3)+
     $                          i2(m+lx+2)*di(m+ly+2,3)*di(m+lz+2,3)+
     $                          i2(m+lx+3)*di(m+ly+3,3)*di(m+lz+3,3)+
     $                          i2(m+lx+4)*di(m+ly+4,3)*di(m+lz+4,3)+
     $                          i2(m+lx+5)*di(m+ly+5,3)*di(m+lz+5,3)+
     $                          i2(m+lx+6)*di(m+ly+6,3)*di(m+lz+6,3)+
     $                          i2(m+lx+7)*di(m+ly+7,3)*di(m+lz+7,3)+
     $                          i2(m+lx+8)*di(m+ly+8,3)*di(m+lz+8,3))
                           d2(45)=d2(45)+dens(prim,n)*
     $                         (i2(m+lx  )*i2(m+ly  )*d2i(m+lz  ,6)+
     $                          i2(m+lx+1)*i2(m+ly+1)*d2i(m+lz+1,6)+
     $                          i2(m+lx+2)*i2(m+ly+2)*d2i(m+lz+2,6)+
     $                          i2(m+lx+3)*i2(m+ly+3)*d2i(m+lz+3,6)+
     $                          i2(m+lx+4)*i2(m+ly+4)*d2i(m+lz+4,6)+
     $                          i2(m+lx+5)*i2(m+ly+5)*d2i(m+lz+5,6)+
     $                          i2(m+lx+6)*i2(m+ly+6)*d2i(m+lz+6,6)+
     $                          i2(m+lx+7)*i2(m+ly+7)*d2i(m+lz+7,6)+
     $                          i2(m+lx+8)*i2(m+ly+8)*d2i(m+lz+8,6))
                           m=m+9
 1302                   continue
                     end if
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
