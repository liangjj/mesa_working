*deck @(#)fd2bnt.f	5.1  11/6/94
      subroutine fd2bnt(i2,dens,lenblk,rv,angmom,imax,jmax,
     #     mmax,lmax,nroots,px,py,pz,lenxyz,mintyp,
     #     nocart,nbtype,nv,energy,di,ndcen,der,leni,d2i,nd2,d2,
     $     npass)
c
c***begin prologue     fd2int
c***date written       840713   (yymmdd)
c***revision date      851113   (yymmdd)
c
c     modification of vfmint from m312
c
c   11 november 1985     pws  at lanl
c          modifying original version from m312 to sum up energy by directly
c          multiplying primitive integrals by two-particle density matrix.
c***keywords
c***author             saxe, paul    (lanl)
c***source             @(#)fd2bnt.f	5.1   11/6/94
c***purpose
c***description
c***references         (none)
c
c***routines called    (none)
c
c***end prologue       fd2int
c
      implicit integer (a-h,o-z)
c
c     real*8 i2(rv,3,0:imax,0:jmax,0:mmax,0:lmax)
      real*8 i2(leni)
      real*8 dens(nv,lenblk)
      real*8 energy
      real*8 di(leni,ndcen)
      real*8 d2i(leni,nd2)
      real*8 der(3,ndcen)
      real*8 d2(78)
      integer angmom(4)
      integer px(lenxyz),py(lenxyz),pz(lenxyz)
      integer mintyp(nbtype),nocart(nbtype)
      parameter (nindex=56)
      integer ipx(nindex),ipy(nindex),ipz(nindex)
      integer jpx(nindex),jpy(nindex),jpz(nindex)
      integer kpx(nindex),kpy(nindex),kpz(nindex)
      integer lpx(nindex),lpy(nindex),lpz(nindex)
c
c     ----- set up minimum and maximum functions in each momentum group -----
c
      mini=mintyp(angmom(1))
      maxi=mini+nocart(angmom(1))-1
      minj=mintyp(angmom(2))
      maxj=minj+nocart(angmom(2))-1
      mink=mintyp(angmom(3))
      maxk=mink+nocart(angmom(3))-1
      minl=mintyp(angmom(4))
      maxl=minl+nocart(angmom(4))-1
c
c     ----- set up offsets -----
c
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
c
                  call lnkerr('bad number of rys roots')
c
    5             continue
                     do 500 m=1,nv
                        energy=energy+i2(m+lx)*i2(m+ly)*
     #                                i2(m+lz)*dens(m,n)
                        der(1,1)=der(1,1)+dens(m,n)*
     #                       di(m+lx,1)*i2(m+ly)*i2(m+lz)
                        der(2,1)=der(2,1)+dens(m,n)*
     #                       i2(m+lx)*di(m+ly,1)*i2(m+lz)
                        der(3,1)=der(3,1)+dens(m,n)*
     #                       i2(m+lx)*i2(m+ly)*di(m+lz,1)
                        d2(1)=d2(1)+dens(m,n)*
     $                       d2i(m+lx,1)*i2(m+ly)*i2(m+lz)
                        d2(2)=d2(2)+dens(m,n)*
     $                       di(m+lx,1)*di(m+ly,1)*i2(m+lz)
                        d2(3)=d2(3)+dens(m,n)*
     $                       i2(m+lx)*d2i(m+ly,1)*i2(m+lz)
                        d2(4)=d2(4)+dens(m,n)*
     $                       di(m+lx,1)*i2(m+ly)*di(m+lz,1)
                        d2(5)=d2(5)+dens(m,n)*
     $                       i2(m+lx)*di(m+ly,1)*di(m+lz,1)
                        d2(6)=d2(6)+dens(m,n)*
     $                       i2(m+lx)*i2(m+ly)*d2i(m+lz,1)
 500                 continue
                     if (npass.lt.4) then
                        do 501 m=1,nv
                           der(1,2)=der(1,2)+dens(m,n)*
     #                         di(m+lx,2)*i2(m+ly)*i2(m+lz)
                           der(2,2)=der(2,2)+dens(m,n)*
     #                         i2(m+lx)*di(m+ly,2)*i2(m+lz)
                           der(3,2)=der(3,2)+dens(m,n)*
     #                         i2(m+lx)*i2(m+ly)*di(m+lz,2)
                           d2(7)=d2(7)+dens(m,n)*
     $                          d2i(m+lx,2)*i2(m+ly)*i2(m+lz)
                           d2(8)=d2(8)+dens(m,n)*
     $                          di(m+lx,2)*di(m+ly,1)*i2(m+lz)
                           d2(9)=d2(9)+dens(m,n)*
     $                          di(m+lx,2)*i2(m+ly)*di(m+lz,1)
                           d2(10)=d2(10)+dens(m,n)*
     $                          d2i(m+lx,3)*i2(m+ly)*i2(m+lz)
                           d2(11)=d2(11)+dens(m,n)*
     $                          di(m+lx,1)*di(m+ly,2)*i2(m+lz)
                           d2(12)=d2(12)+dens(m,n)*
     $                          i2(m+lx)*d2i(m+ly,2)*i2(m+lz)
                           d2(13)=d2(13)+dens(m,n)*
     $                          i2(m+lx)*di(m+ly,2)*di(m+lz,1)
                           d2(14)=d2(14)+dens(m,n)*
     $                          di(m+lx,2)*di(m+ly,2)*i2(m+lz)
                           d2(15)=d2(15)+dens(m,n)*
     $                          i2(m+lx)*d2i(m+ly,3)*i2(m+lz)
                           d2(16)=d2(16)+dens(m,n)*
     $                          di(m+lx,1)*i2(m+ly)*di(m+lz,2)
                           d2(17)=d2(17)+dens(m,n)*
     $                          i2(m+lx)*di(m+ly,1)*di(m+lz,2)
                           d2(18)=d2(18)+dens(m,n)*
     $                          i2(m+lx)*i2(m+ly)*d2i(m+lz,2)
                           d2(19)=d2(19)+dens(m,n)*
     $                          di(m+lx,2)*i2(m+ly)*di(m+lz,2)
                           d2(20)=d2(20)+dens(m,n)*
     $                          i2(m+lx)*di(m+ly,2)*di(m+lz,2)
                           d2(21)=d2(21)+dens(m,n)*
     $                          i2(m+lx)*i2(m+ly)*d2i(m+lz,3)
 501                    continue
                     end if
                     if (npass.eq.1) then
                        do 502 m=1,nv
                           der(1,3)=der(1,3)+dens(m,n)*
     #                         di(m+lx,3)*i2(m+ly)*i2(m+lz)
                           der(2,3)=der(2,3)+dens(m,n)*
     #                         i2(m+lx)*di(m+ly,3)*i2(m+lz)
                           der(3,3)=der(3,3)+dens(m,n)*
     #                         i2(m+lx)*i2(m+ly)*di(m+lz,3)
                           d2(22)=d2(22)+dens(m,n)*
     $                          d2i(m+lx,4)*i2(m+ly)*i2(m+lz)
                           d2(23)=d2(23)+dens(m,n)*
     $                          di(m+lx,3)*di(m+ly,1)*i2(m+lz)
                           d2(24)=d2(24)+dens(m,n)*
     $                          di(m+lx,3)*i2(m+ly)*di(m+lz,1)
                           d2(25)=d2(25)+dens(m,n)*
     $                          d2i(m+lx,5)*i2(m+ly)*i2(m+lz)
                           d2(26)=d2(26)+dens(m,n)*
     $                          di(m+lx,3)*di(m+ly,2)*i2(m+lz)
                           d2(27)=d2(27)+dens(m,n)*
     $                          di(m+lx,3)*i2(m+ly)*di(m+lz,2)
                           d2(28)=d2(28)+dens(m,n)*
     $                          d2i(m+lx,6)*i2(m+ly)*i2(m+lz)
                           d2(29)=d2(29)+dens(m,n)*
     $                          di(m+lx,1)*di(m+ly,3)*i2(m+lz)
                           d2(30)=d2(30)+dens(m,n)*
     $                          i2(m+lx)*d2i(m+ly,4)*i2(m+lz)
                           d2(31)=d2(31)+dens(m,n)*
     $                          i2(m+lx)*di(m+ly,3)*di(m+lz,1)
                           d2(32)=d2(32)+dens(m,n)*
     $                          di(m+lx,2)*di(m+ly,3)*i2(m+lz)
                           d2(33)=d2(33)+dens(m,n)*
     $                          i2(m+lx)*d2i(m+ly,5)*i2(m+lz)
                           d2(34)=d2(34)+dens(m,n)*
     $                          i2(m+lx)*di(m+ly,3)*di(m+lz,2)
                           d2(35)=d2(35)+dens(m,n)*
     $                          di(m+lx,3)*di(m+ly,3)*i2(m+lz)
                           d2(36)=d2(36)+dens(m,n)*
     $                          i2(m+lx)*d2i(m+ly,6)*i2(m+lz)
                           d2(37)=d2(37)+dens(m,n)*
     $                          di(m+lx,1)*i2(m+ly)*di(m+lz,3)
                           d2(38)=d2(38)+dens(m,n)*
     $                          i2(m+lx)*di(m+ly,1)*di(m+lz,3)
                           d2(39)=d2(39)+dens(m,n)*
     $                          i2(m+lx)*i2(m+ly)*d2i(m+lz,4)
                           d2(40)=d2(40)+dens(m,n)*
     $                          di(m+lx,2)*i2(m+ly)*di(m+lz,3)
                           d2(41)=d2(41)+dens(m,n)*
     $                          i2(m+lx)*di(m+ly,2)*di(m+lz,3)
                           d2(42)=d2(42)+dens(m,n)*
     $                          i2(m+lx)*i2(m+ly)*d2i(m+lz,5)
                           d2(43)=d2(43)+dens(m,n)*
     $                          di(m+lx,3)*i2(m+ly)*di(m+lz,3)
                           d2(44)=d2(44)+dens(m,n)*
     $                          i2(m+lx)*di(m+ly,3)*di(m+lz,3)
                           d2(45)=d2(45)+dens(m,n)*
     $                          i2(m+lx)*i2(m+ly)*d2i(m+lz,6)
 502                    continue
                     end if
                     go to 15
    6             continue
                     m=1
                     do 600 prim=1,nv
                        energy=energy+dens(prim,n)*
     $                      (i2(m+lx  )*i2(m+ly  )*i2(m+lz  )+
     $                       i2(m+lx+1)*i2(m+ly+1)*i2(m+lz+1))
                        der(1,1)=der(1,1)+dens(prim,n)*
     #                      (di(m+lx  ,1)*i2(m+ly  )*i2(m+lz  )+
     #                       di(m+lx+1,1)*i2(m+ly+1)*i2(m+lz+1))
                        der(2,1)=der(2,1)+dens(prim,n)*
     #                      (i2(m+lx  )*di(m+ly  ,1)*i2(m+lz  )+
     #                       i2(m+lx+1)*di(m+ly+1,1)*i2(m+lz+1))
                        der(3,1)=der(3,1)+dens(prim,n)*
     #                      (i2(m+lx  )*i2(m+ly  )*di(m+lz  ,1)+
     #                       i2(m+lx+1)*i2(m+ly+1)*di(m+lz+1,1))
                        d2(1)=d2(1)+dens(prim,n)*
     $                      (d2i(m+lx  ,1)*i2(m+ly  )*i2(m+lz  )+
     $                       d2i(m+lx+1,1)*i2(m+ly+1)*i2(m+lz+1))
                        d2(2)=d2(2)+dens(prim,n)*
     $                      (di(m+lx  ,1)*di(m+ly  ,1)*i2(m+lz  )+
     $                       di(m+lx+1,1)*di(m+ly+1,1)*i2(m+lz+1))
                        d2(3)=d2(3)+dens(prim,n)*
     $                      (i2(m+lx  )*d2i(m+ly  ,1)*i2(m+lz  )+
     $                       i2(m+lx+1)*d2i(m+ly+1,1)*i2(m+lz+1))
                        d2(4)=d2(4)+dens(prim,n)*
     $                      (di(m+lx  ,1)*i2(m+ly  )*di(m+lz  ,1)+
     $                       di(m+lx+1,1)*i2(m+ly+1)*di(m+lz+1,1))
                        d2(5)=d2(5)+dens(prim,n)*
     $                      (i2(m+lx  )*di(m+ly  ,1)*di(m+lz  ,1)+
     $                       i2(m+lx+1)*di(m+ly+1,1)*di(m+lz+1,1))
                        d2(6)=d2(6)+dens(prim,n)*
     $                      (i2(m+lx  )*i2(m+ly  )*d2i(m+lz  ,1)+
     $                       i2(m+lx+1)*i2(m+ly+1)*d2i(m+lz+1,1))
                        m=m+2
 600                 continue
                     if (npass.lt.4) then
                        m=1
                        do 601 prim=1,nv
                           der(1,2)=der(1,2)+dens(prim,n)*
     #                        (di(m+lx  ,2)*i2(m+ly  )*i2(m+lz  )+
     #                         di(m+lx+1,2)*i2(m+ly+1)*i2(m+lz+1))
                           der(2,2)=der(2,2)+dens(prim,n)*
     #                        (i2(m+lx  )*di(m+ly  ,2)*i2(m+lz  )+
     #                         i2(m+lx+1)*di(m+ly+1,2)*i2(m+lz+1))
                           der(3,2)=der(3,2)+dens(prim,n)*
     #                        (i2(m+lx  )*i2(m+ly  )*di(m+lz  ,2)+
     #                         i2(m+lx+1)*i2(m+ly+1)*di(m+lz+1,2))
                           d2(7)=d2(7)+dens(prim,n)*
     $                         (d2i(m+lx  ,2)*i2(m+ly  )*i2(m+lz  )+
     $                          d2i(m+lx+1,2)*i2(m+ly+1)*i2(m+lz+1))
                           d2(8)=d2(8)+dens(prim,n)*
     $                         (di(m+lx  ,2)*di(m+ly  ,1)*i2(m+lz  )+
     $                          di(m+lx+1,2)*di(m+ly+1,1)*i2(m+lz+1))
                           d2(9)=d2(9)+dens(prim,n)*
     $                         (di(m+lx  ,2)*i2(m+ly  )*di(m+lz  ,1)+
     $                          di(m+lx+1,2)*i2(m+ly+1)*di(m+lz+1,1))
                           d2(10)=d2(10)+dens(prim,n)*
     $                         (d2i(m+lx  ,3)*i2(m+ly  )*i2(m+lz  )+
     $                          d2i(m+lx+1,3)*i2(m+ly+1)*i2(m+lz+1))
                           d2(11)=d2(11)+dens(prim,n)*
     $                         (di(m+lx  ,1)*di(m+ly  ,2)*i2(m+lz  )+
     $                          di(m+lx+1,1)*di(m+ly+1,2)*i2(m+lz+1))
                           d2(12)=d2(12)+dens(prim,n)*
     $                         (i2(m+lx  )*d2i(m+ly  ,2)*i2(m+lz  )+
     $                          i2(m+lx+1)*d2i(m+ly+1,2)*i2(m+lz+1))
                           d2(13)=d2(13)+dens(prim,n)*
     $                         (i2(m+lx  )*di(m+ly  ,2)*di(m+lz  ,1)+
     $                          i2(m+lx+1)*di(m+ly+1,2)*di(m+lz+1,1))
                           d2(14)=d2(14)+dens(prim,n)*
     $                         (di(m+lx  ,2)*di(m+ly  ,2)*i2(m+lz  )+
     $                          di(m+lx+1,2)*di(m+ly+1,2)*i2(m+lz+1))
                           d2(15)=d2(15)+dens(prim,n)*
     $                         (i2(m+lx  )*d2i(m+ly  ,3)*i2(m+lz  )+
     $                          i2(m+lx+1)*d2i(m+ly+1,3)*i2(m+lz+1))
                           d2(16)=d2(16)+dens(prim,n)*
     $                         (di(m+lx  ,1)*i2(m+ly  )*di(m+lz  ,2)+
     $                          di(m+lx+1,1)*i2(m+ly+1)*di(m+lz+1,2))
                           d2(17)=d2(17)+dens(prim,n)*
     $                         (i2(m+lx  )*di(m+ly  ,1)*di(m+lz  ,2)+
     $                          i2(m+lx+1)*di(m+ly+1,1)*di(m+lz+1,2))
                           d2(18)=d2(18)+dens(prim,n)*
     $                         (i2(m+lx  )*i2(m+ly  )*d2i(m+lz  ,2)+
     $                          i2(m+lx+1)*i2(m+ly+1)*d2i(m+lz+1,2))
                           d2(19)=d2(19)+dens(prim,n)*
     $                         (di(m+lx  ,2)*i2(m+ly  )*di(m+lz  ,2)+
     $                          di(m+lx+1,2)*i2(m+ly+1)*di(m+lz+1,2))
                           d2(20)=d2(20)+dens(prim,n)*
     $                         (i2(m+lx  )*di(m+ly  ,2)*di(m+lz  ,2)+
     $                          i2(m+lx+1)*di(m+ly+1,2)*di(m+lz+1,2))
                           d2(21)=d2(21)+dens(prim,n)*
     $                         (i2(m+lx  )*i2(m+ly  )*d2i(m+lz  ,3)+
     $                          i2(m+lx+1)*i2(m+ly+1)*d2i(m+lz+1,3))
                           m=m+2
 601                    continue
                     end if
                     if (npass.eq.1) then
                        m=1
                        do 602 prim=1,nv
                           der(1,3)=der(1,3)+dens(prim,n)*
     #                        (di(m+lx  ,3)*i2(m+ly  )*i2(m+lz  )+
     #                         di(m+lx+1,3)*i2(m+ly+1)*i2(m+lz+1))
                           der(2,3)=der(2,3)+dens(prim,n)*
     #                        (i2(m+lx  )*di(m+ly  ,3)*i2(m+lz  )+
     #                         i2(m+lx+1)*di(m+ly+1,3)*i2(m+lz+1))
                           der(3,3)=der(3,3)+dens(prim,n)*
     #                        (i2(m+lx  )*i2(m+ly  )*di(m+lz  ,3)+
     #                         i2(m+lx+1)*i2(m+ly+1)*di(m+lz+1,3))
                           d2(22)=d2(22)+dens(prim,n)*
     $                         (d2i(m+lx  ,4)*i2(m+ly  )*i2(m+lz  )+
     $                          d2i(m+lx+1,4)*i2(m+ly+1)*i2(m+lz+1))
                           d2(23)=d2(23)+dens(prim,n)*
     $                         (di(m+lx  ,3)*di(m+ly  ,1)*i2(m+lz  )+
     $                          di(m+lx+1,3)*di(m+ly+1,1)*i2(m+lz+1))
                           d2(24)=d2(24)+dens(prim,n)*
     $                         (di(m+lx  ,3)*i2(m+ly  )*di(m+lz  ,1)+
     $                          di(m+lx+1,3)*i2(m+ly+1)*di(m+lz+1,1))
                           d2(25)=d2(25)+dens(prim,n)*
     $                         (d2i(m+lx  ,5)*i2(m+ly  )*i2(m+lz  )+
     $                          d2i(m+lx+1,5)*i2(m+ly+1)*i2(m+lz+1))
                           d2(26)=d2(26)+dens(prim,n)*
     $                         (di(m+lx  ,3)*di(m+ly  ,2)*i2(m+lz  )+
     $                          di(m+lx+1,3)*di(m+ly+1,2)*i2(m+lz+1))
                           d2(27)=d2(27)+dens(prim,n)*
     $                         (di(m+lx  ,3)*i2(m+ly  )*di(m+lz  ,2)+
     $                          di(m+lx+1,3)*i2(m+ly+1)*di(m+lz+1,2))
                           d2(28)=d2(28)+dens(prim,n)*
     $                         (d2i(m+lx  ,6)*i2(m+ly  )*i2(m+lz  )+
     $                          d2i(m+lx+1,6)*i2(m+ly+1)*i2(m+lz+1))
                           d2(29)=d2(29)+dens(prim,n)*
     $                         (di(m+lx  ,1)*di(m+ly  ,3)*i2(m+lz  )+
     $                          di(m+lx+1,1)*di(m+ly+1,3)*i2(m+lz+1))
                           d2(30)=d2(30)+dens(prim,n)*
     $                         (i2(m+lx  )*d2i(m+ly  ,4)*i2(m+lz  )+
     $                          i2(m+lx+1)*d2i(m+ly+1,4)*i2(m+lz+1))
                           d2(31)=d2(31)+dens(prim,n)*
     $                         (i2(m+lx  )*di(m+ly  ,3)*di(m+lz  ,1)+
     $                          i2(m+lx+1)*di(m+ly+1,3)*di(m+lz+1,1))
                           d2(32)=d2(32)+dens(prim,n)*
     $                         (di(m+lx  ,2)*di(m+ly  ,3)*i2(m+lz  )+
     $                          di(m+lx+1,2)*di(m+ly+1,3)*i2(m+lz+1))
                           d2(33)=d2(33)+dens(prim,n)*
     $                         (i2(m+lx  )*d2i(m+ly  ,5)*i2(m+lz  )+
     $                          i2(m+lx+1)*d2i(m+ly+1,5)*i2(m+lz+1))
                           d2(34)=d2(34)+dens(prim,n)*
     $                         (i2(m+lx  )*di(m+ly  ,3)*di(m+lz  ,2)+
     $                          i2(m+lx+1)*di(m+ly+1,3)*di(m+lz+1,2))
                           d2(35)=d2(35)+dens(prim,n)*
     $                         (di(m+lx  ,3)*di(m+ly  ,3)*i2(m+lz  )+
     $                          di(m+lx+1,3)*di(m+ly+1,3)*i2(m+lz+1))
                           d2(36)=d2(36)+dens(prim,n)*
     $                         (i2(m+lx  )*d2i(m+ly  ,6)*i2(m+lz  )+
     $                          i2(m+lx+1)*d2i(m+ly+1,6)*i2(m+lz+1))
                           d2(37)=d2(37)+dens(prim,n)*
     $                         (di(m+lx  ,1)*i2(m+ly  )*di(m+lz  ,3)+
     $                          di(m+lx+1,1)*i2(m+ly+1)*di(m+lz+1,3))
                           d2(38)=d2(38)+dens(prim,n)*
     $                         (i2(m+lx  )*di(m+ly  ,1)*di(m+lz  ,3)+
     $                          i2(m+lx+1)*di(m+ly+1,1)*di(m+lz+1,3))
                           d2(39)=d2(39)+dens(prim,n)*
     $                         (i2(m+lx  )*i2(m+ly  )*d2i(m+lz  ,4)+
     $                          i2(m+lx+1)*i2(m+ly+1)*d2i(m+lz+1,4))
                           d2(40)=d2(40)+dens(prim,n)*
     $                         (di(m+lx  ,2)*i2(m+ly  )*di(m+lz  ,3)+
     $                          di(m+lx+1,2)*i2(m+ly+1)*di(m+lz+1,3))
                           d2(41)=d2(41)+dens(prim,n)*
     $                         (i2(m+lx  )*di(m+ly  ,2)*di(m+lz  ,3)+
     $                          i2(m+lx+1)*di(m+ly+1,2)*di(m+lz+1,3))
                           d2(42)=d2(42)+dens(prim,n)*
     $                         (i2(m+lx  )*i2(m+ly  )*d2i(m+lz  ,5)+
     $                          i2(m+lx+1)*i2(m+ly+1)*d2i(m+lz+1,5))
                           d2(43)=d2(43)+dens(prim,n)*
     $                         (di(m+lx  ,3)*i2(m+ly  )*di(m+lz  ,3)+
     $                          di(m+lx+1,3)*i2(m+ly+1)*di(m+lz+1,3))
                           d2(44)=d2(44)+dens(prim,n)*
     $                         (i2(m+lx  )*di(m+ly  ,3)*di(m+lz  ,3)+
     $                          i2(m+lx+1)*di(m+ly+1,3)*di(m+lz+1,3))
                           d2(45)=d2(45)+dens(prim,n)*
     $                         (i2(m+lx  )*i2(m+ly  )*d2i(m+lz  ,6)+
     $                          i2(m+lx+1)*i2(m+ly+1)*d2i(m+lz+1,6))
                           m=m+2
 602                    continue
                     end if
                     go to 15
    7             continue
                     m=1
                     do 700 prim=1,nv
                        energy=energy+dens(prim,n)*
     $                      (i2(m+lx  )*i2(m+ly  )*i2(m+lz  )+
     $                       i2(m+lx+1)*i2(m+ly+1)*i2(m+lz+1)+
     $                       i2(m+lx+2)*i2(m+ly+2)*i2(m+lz+2))
                        der(1,1)=der(1,1)+dens(prim,n)*
     #                      (di(m+lx  ,1)*i2(m+ly  )*i2(m+lz  )+
     #                       di(m+lx+1,1)*i2(m+ly+1)*i2(m+lz+1)+
     #                       di(m+lx+2,1)*i2(m+ly+2)*i2(m+lz+2))
                        der(2,1)=der(2,1)+dens(prim,n)*
     #                      (i2(m+lx  )*di(m+ly  ,1)*i2(m+lz  )+
     #                       i2(m+lx+1)*di(m+ly+1,1)*i2(m+lz+1)+
     #                       i2(m+lx+2)*di(m+ly+2,1)*i2(m+lz+2))
                        der(3,1)=der(3,1)+dens(prim,n)*
     #                      (i2(m+lx  )*i2(m+ly  )*di(m+lz  ,1)+
     #                       i2(m+lx+1)*i2(m+ly+1)*di(m+lz+1,1)+
     #                       i2(m+lx+2)*i2(m+ly+2)*di(m+lz+2,1))
                        d2(1)=d2(1)+dens(prim,n)*
     $                      (d2i(m+lx  ,1)*i2(m+ly  )*i2(m+lz  )+
     $                       d2i(m+lx+1,1)*i2(m+ly+1)*i2(m+lz+1)+
     $                       d2i(m+lx+2,1)*i2(m+ly+2)*i2(m+lz+2))
                        d2(2)=d2(2)+dens(prim,n)*
     $                      (di(m+lx  ,1)*di(m+ly  ,1)*i2(m+lz  )+
     $                       di(m+lx+1,1)*di(m+ly+1,1)*i2(m+lz+1)+
     $                       di(m+lx+2,1)*di(m+ly+2,1)*i2(m+lz+2))
                        d2(3)=d2(3)+dens(prim,n)*
     $                      (i2(m+lx  )*d2i(m+ly  ,1)*i2(m+lz  )+
     $                       i2(m+lx+1)*d2i(m+ly+1,1)*i2(m+lz+1)+
     $                       i2(m+lx+2)*d2i(m+ly+2,1)*i2(m+lz+2))
                        d2(4)=d2(4)+dens(prim,n)*
     $                      (di(m+lx  ,1)*i2(m+ly  )*di(m+lz  ,1)+
     $                       di(m+lx+1,1)*i2(m+ly+1)*di(m+lz+1,1)+
     $                       di(m+lx+2,1)*i2(m+ly+2)*di(m+lz+2,1))
                        d2(5)=d2(5)+dens(prim,n)*
     $                      (i2(m+lx  )*di(m+ly  ,1)*di(m+lz  ,1)+
     $                       i2(m+lx+1)*di(m+ly+1,1)*di(m+lz+1,1)+
     $                       i2(m+lx+2)*di(m+ly+2,1)*di(m+lz+2,1))
                        d2(6)=d2(6)+dens(prim,n)*
     $                      (i2(m+lx  )*i2(m+ly  )*d2i(m+lz  ,1)+
     $                       i2(m+lx+1)*i2(m+ly+1)*d2i(m+lz+1,1)+
     $                       i2(m+lx+2)*i2(m+ly+2)*d2i(m+lz+2,1))
                        m=m+3
 700                 continue
                     if (npass.lt.4) then
                        m=1
                        do 701 prim=1,nv
                           der(1,2)=der(1,2)+dens(prim,n)*
     #                        (di(m+lx  ,2)*i2(m+ly  )*i2(m+lz  )+
     #                         di(m+lx+1,2)*i2(m+ly+1)*i2(m+lz+1)+
     #                         di(m+lx+2,2)*i2(m+ly+2)*i2(m+lz+2))
                           der(2,2)=der(2,2)+dens(prim,n)*
     #                        (i2(m+lx  )*di(m+ly  ,2)*i2(m+lz  )+
     #                         i2(m+lx+1)*di(m+ly+1,2)*i2(m+lz+1)+
     #                         i2(m+lx+2)*di(m+ly+2,2)*i2(m+lz+2))
                           der(3,2)=der(3,2)+dens(prim,n)*
     #                        (i2(m+lx  )*i2(m+ly  )*di(m+lz  ,2)+
     #                         i2(m+lx+1)*i2(m+ly+1)*di(m+lz+1,2)+
     #                         i2(m+lx+2)*i2(m+ly+2)*di(m+lz+2,2))
                           d2(7)=d2(7)+dens(prim,n)*
     $                         (d2i(m+lx  ,2)*i2(m+ly  )*i2(m+lz  )+
     $                          d2i(m+lx+1,2)*i2(m+ly+1)*i2(m+lz+1)+
     $                          d2i(m+lx+2,2)*i2(m+ly+2)*i2(m+lz+2))
                           d2(8)=d2(8)+dens(prim,n)*
     $                         (di(m+lx  ,2)*di(m+ly  ,1)*i2(m+lz  )+
     $                          di(m+lx+1,2)*di(m+ly+1,1)*i2(m+lz+1)+
     $                          di(m+lx+2,2)*di(m+ly+2,1)*i2(m+lz+2))
                           d2(9)=d2(9)+dens(prim,n)*
     $                         (di(m+lx  ,2)*i2(m+ly  )*di(m+lz  ,1)+
     $                          di(m+lx+1,2)*i2(m+ly+1)*di(m+lz+1,1)+
     $                          di(m+lx+2,2)*i2(m+ly+2)*di(m+lz+2,1))
                           d2(10)=d2(10)+dens(prim,n)*
     $                         (d2i(m+lx  ,3)*i2(m+ly  )*i2(m+lz  )+
     $                          d2i(m+lx+1,3)*i2(m+ly+1)*i2(m+lz+1)+
     $                          d2i(m+lx+2,3)*i2(m+ly+2)*i2(m+lz+2))
                           d2(11)=d2(11)+dens(prim,n)*
     $                         (di(m+lx  ,1)*di(m+ly  ,2)*i2(m+lz  )+
     $                          di(m+lx+1,1)*di(m+ly+1,2)*i2(m+lz+1)+
     $                          di(m+lx+2,1)*di(m+ly+2,2)*i2(m+lz+2))
                           d2(12)=d2(12)+dens(prim,n)*
     $                         (i2(m+lx  )*d2i(m+ly  ,2)*i2(m+lz  )+
     $                          i2(m+lx+1)*d2i(m+ly+1,2)*i2(m+lz+1)+
     $                          i2(m+lx+2)*d2i(m+ly+2,2)*i2(m+lz+2))
                           d2(13)=d2(13)+dens(prim,n)*
     $                         (i2(m+lx  )*di(m+ly  ,2)*di(m+lz  ,1)+
     $                          i2(m+lx+1)*di(m+ly+1,2)*di(m+lz+1,1)+
     $                          i2(m+lx+2)*di(m+ly+2,2)*di(m+lz+2,1))
                           d2(14)=d2(14)+dens(prim,n)*
     $                         (di(m+lx  ,2)*di(m+ly  ,2)*i2(m+lz  )+
     $                          di(m+lx+1,2)*di(m+ly+1,2)*i2(m+lz+1)+
     $                          di(m+lx+2,2)*di(m+ly+2,2)*i2(m+lz+2))
                           d2(15)=d2(15)+dens(prim,n)*
     $                         (i2(m+lx  )*d2i(m+ly  ,3)*i2(m+lz  )+
     $                          i2(m+lx+1)*d2i(m+ly+1,3)*i2(m+lz+1)+
     $                          i2(m+lx+2)*d2i(m+ly+2,3)*i2(m+lz+2))
                           d2(16)=d2(16)+dens(prim,n)*
     $                         (di(m+lx  ,1)*i2(m+ly  )*di(m+lz  ,2)+
     $                          di(m+lx+1,1)*i2(m+ly+1)*di(m+lz+1,2)+
     $                          di(m+lx+2,1)*i2(m+ly+2)*di(m+lz+2,2))
                           d2(17)=d2(17)+dens(prim,n)*
     $                         (i2(m+lx  )*di(m+ly  ,1)*di(m+lz  ,2)+
     $                          i2(m+lx+1)*di(m+ly+1,1)*di(m+lz+1,2)+
     $                          i2(m+lx+2)*di(m+ly+2,1)*di(m+lz+2,2))
                           d2(18)=d2(18)+dens(prim,n)*
     $                         (i2(m+lx  )*i2(m+ly  )*d2i(m+lz  ,2)+
     $                          i2(m+lx+1)*i2(m+ly+1)*d2i(m+lz+1,2)+
     $                          i2(m+lx+2)*i2(m+ly+2)*d2i(m+lz+2,2))
                           d2(19)=d2(19)+dens(prim,n)*
     $                         (di(m+lx  ,2)*i2(m+ly  )*di(m+lz  ,2)+
     $                          di(m+lx+1,2)*i2(m+ly+1)*di(m+lz+1,2)+
     $                          di(m+lx+2,2)*i2(m+ly+2)*di(m+lz+2,2))
                           d2(20)=d2(20)+dens(prim,n)*
     $                         (i2(m+lx  )*di(m+ly  ,2)*di(m+lz  ,2)+
     $                          i2(m+lx+1)*di(m+ly+1,2)*di(m+lz+1,2)+
     $                          i2(m+lx+2)*di(m+ly+2,2)*di(m+lz+2,2))
                           d2(21)=d2(21)+dens(prim,n)*
     $                         (i2(m+lx  )*i2(m+ly  )*d2i(m+lz  ,3)+
     $                          i2(m+lx+1)*i2(m+ly+1)*d2i(m+lz+1,3)+
     $                          i2(m+lx+2)*i2(m+ly+2)*d2i(m+lz+2,3))
                           m=m+3
 701                    continue
                     end if
                     if (npass.eq.1) then
                        m=1
                        do 702 prim=1,nv
                           der(1,3)=der(1,3)+dens(prim,n)*
     #                        (di(m+lx  ,3)*i2(m+ly  )*i2(m+lz  )+
     #                         di(m+lx+1,3)*i2(m+ly+1)*i2(m+lz+1)+
     #                         di(m+lx+2,3)*i2(m+ly+2)*i2(m+lz+2))
                           der(2,3)=der(2,3)+dens(prim,n)*
     #                        (i2(m+lx  )*di(m+ly  ,3)*i2(m+lz  )+
     #                         i2(m+lx+1)*di(m+ly+1,3)*i2(m+lz+1)+
     #                         i2(m+lx+2)*di(m+ly+2,3)*i2(m+lz+2))
                           der(3,3)=der(3,3)+dens(prim,n)*
     #                        (i2(m+lx  )*i2(m+ly  )*di(m+lz  ,3)+
     #                         i2(m+lx+1)*i2(m+ly+1)*di(m+lz+1,3)+
     #                         i2(m+lx+2)*i2(m+ly+2)*di(m+lz+2,3))
                           d2(22)=d2(22)+dens(prim,n)*
     $                         (d2i(m+lx  ,4)*i2(m+ly  )*i2(m+lz  )+
     $                          d2i(m+lx+1,4)*i2(m+ly+1)*i2(m+lz+1)+
     $                          d2i(m+lx+2,4)*i2(m+ly+2)*i2(m+lz+2))
                           d2(23)=d2(23)+dens(prim,n)*
     $                         (di(m+lx  ,3)*di(m+ly  ,1)*i2(m+lz  )+
     $                          di(m+lx+1,3)*di(m+ly+1,1)*i2(m+lz+1)+
     $                          di(m+lx+2,3)*di(m+ly+2,1)*i2(m+lz+2))
                           d2(24)=d2(24)+dens(prim,n)*
     $                         (di(m+lx  ,3)*i2(m+ly  )*di(m+lz  ,1)+
     $                          di(m+lx+1,3)*i2(m+ly+1)*di(m+lz+1,1)+
     $                          di(m+lx+2,3)*i2(m+ly+2)*di(m+lz+2,1))
                           d2(25)=d2(25)+dens(prim,n)*
     $                         (d2i(m+lx  ,5)*i2(m+ly  )*i2(m+lz  )+
     $                          d2i(m+lx+1,5)*i2(m+ly+1)*i2(m+lz+1)+
     $                          d2i(m+lx+2,5)*i2(m+ly+2)*i2(m+lz+2))
                           d2(26)=d2(26)+dens(prim,n)*
     $                         (di(m+lx  ,3)*di(m+ly  ,2)*i2(m+lz  )+
     $                          di(m+lx+1,3)*di(m+ly+1,2)*i2(m+lz+1)+
     $                          di(m+lx+2,3)*di(m+ly+2,2)*i2(m+lz+2))
                           d2(27)=d2(27)+dens(prim,n)*
     $                         (di(m+lx  ,3)*i2(m+ly  )*di(m+lz  ,2)+
     $                          di(m+lx+1,3)*i2(m+ly+1)*di(m+lz+1,2)+
     $                          di(m+lx+2,3)*i2(m+ly+2)*di(m+lz+2,2))
                           d2(28)=d2(28)+dens(prim,n)*
     $                         (d2i(m+lx  ,6)*i2(m+ly  )*i2(m+lz  )+
     $                          d2i(m+lx+1,6)*i2(m+ly+1)*i2(m+lz+1)+
     $                          d2i(m+lx+2,6)*i2(m+ly+2)*i2(m+lz+2))
                           d2(29)=d2(29)+dens(prim,n)*
     $                         (di(m+lx  ,1)*di(m+ly  ,3)*i2(m+lz  )+
     $                          di(m+lx+1,1)*di(m+ly+1,3)*i2(m+lz+1)+
     $                          di(m+lx+2,1)*di(m+ly+2,3)*i2(m+lz+2))
                           d2(30)=d2(30)+dens(prim,n)*
     $                         (i2(m+lx  )*d2i(m+ly  ,4)*i2(m+lz  )+
     $                          i2(m+lx+1)*d2i(m+ly+1,4)*i2(m+lz+1)+
     $                          i2(m+lx+2)*d2i(m+ly+2,4)*i2(m+lz+2))
                           d2(31)=d2(31)+dens(prim,n)*
     $                         (i2(m+lx  )*di(m+ly  ,3)*di(m+lz  ,1)+
     $                          i2(m+lx+1)*di(m+ly+1,3)*di(m+lz+1,1)+
     $                          i2(m+lx+2)*di(m+ly+2,3)*di(m+lz+2,1))
                           d2(32)=d2(32)+dens(prim,n)*
     $                         (di(m+lx  ,2)*di(m+ly  ,3)*i2(m+lz  )+
     $                          di(m+lx+1,2)*di(m+ly+1,3)*i2(m+lz+1)+
     $                          di(m+lx+2,2)*di(m+ly+2,3)*i2(m+lz+2))
                           d2(33)=d2(33)+dens(prim,n)*
     $                         (i2(m+lx  )*d2i(m+ly  ,5)*i2(m+lz  )+
     $                          i2(m+lx+1)*d2i(m+ly+1,5)*i2(m+lz+1)+
     $                          i2(m+lx+2)*d2i(m+ly+2,5)*i2(m+lz+2))
                           d2(34)=d2(34)+dens(prim,n)*
     $                         (i2(m+lx  )*di(m+ly  ,3)*di(m+lz  ,2)+
     $                          i2(m+lx+1)*di(m+ly+1,3)*di(m+lz+1,2)+
     $                          i2(m+lx+2)*di(m+ly+2,3)*di(m+lz+2,2))
                           d2(35)=d2(35)+dens(prim,n)*
     $                         (di(m+lx  ,3)*di(m+ly  ,3)*i2(m+lz  )+
     $                          di(m+lx+1,3)*di(m+ly+1,3)*i2(m+lz+1)+
     $                          di(m+lx+2,3)*di(m+ly+2,3)*i2(m+lz+2))
                           d2(36)=d2(36)+dens(prim,n)*
     $                         (i2(m+lx  )*d2i(m+ly  ,6)*i2(m+lz  )+
     $                          i2(m+lx+1)*d2i(m+ly+1,6)*i2(m+lz+1)+
     $                          i2(m+lx+2)*d2i(m+ly+2,6)*i2(m+lz+2))
                           d2(37)=d2(37)+dens(prim,n)*
     $                         (di(m+lx  ,1)*i2(m+ly  )*di(m+lz  ,3)+
     $                          di(m+lx+1,1)*i2(m+ly+1)*di(m+lz+1,3)+
     $                          di(m+lx+2,1)*i2(m+ly+2)*di(m+lz+2,3))
                           d2(38)=d2(38)+dens(prim,n)*
     $                         (i2(m+lx  )*di(m+ly  ,1)*di(m+lz  ,3)+
     $                          i2(m+lx+1)*di(m+ly+1,1)*di(m+lz+1,3)+
     $                          i2(m+lx+2)*di(m+ly+2,1)*di(m+lz+2,3))
                           d2(39)=d2(39)+dens(prim,n)*
     $                         (i2(m+lx  )*i2(m+ly  )*d2i(m+lz  ,4)+
     $                          i2(m+lx+1)*i2(m+ly+1)*d2i(m+lz+1,4)+
     $                          i2(m+lx+2)*i2(m+ly+2)*d2i(m+lz+2,4))
                           d2(40)=d2(40)+dens(prim,n)*
     $                         (di(m+lx  ,2)*i2(m+ly  )*di(m+lz  ,3)+
     $                          di(m+lx+1,2)*i2(m+ly+1)*di(m+lz+1,3)+
     $                          di(m+lx+2,2)*i2(m+ly+2)*di(m+lz+2,3))
                           d2(41)=d2(41)+dens(prim,n)*
     $                         (i2(m+lx  )*di(m+ly  ,2)*di(m+lz  ,3)+
     $                          i2(m+lx+1)*di(m+ly+1,2)*di(m+lz+1,3)+
     $                          i2(m+lx+2)*di(m+ly+2,2)*di(m+lz+2,3))
                           d2(42)=d2(42)+dens(prim,n)*
     $                         (i2(m+lx  )*i2(m+ly  )*d2i(m+lz  ,5)+
     $                          i2(m+lx+1)*i2(m+ly+1)*d2i(m+lz+1,5)+
     $                          i2(m+lx+2)*i2(m+ly+2)*d2i(m+lz+2,5))
                           d2(43)=d2(43)+dens(prim,n)*
     $                         (di(m+lx  ,3)*i2(m+ly  )*di(m+lz  ,3)+
     $                          di(m+lx+1,3)*i2(m+ly+1)*di(m+lz+1,3)+
     $                          di(m+lx+2,3)*i2(m+ly+2)*di(m+lz+2,3))
                           d2(44)=d2(44)+dens(prim,n)*
     $                         (i2(m+lx  )*di(m+ly  ,3)*di(m+lz  ,3)+
     $                          i2(m+lx+1)*di(m+ly+1,3)*di(m+lz+1,3)+
     $                          i2(m+lx+2)*di(m+ly+2,3)*di(m+lz+2,3))
                           d2(45)=d2(45)+dens(prim,n)*
     $                         (i2(m+lx  )*i2(m+ly  )*d2i(m+lz  ,6)+
     $                          i2(m+lx+1)*i2(m+ly+1)*d2i(m+lz+1,6)+
     $                          i2(m+lx+2)*i2(m+ly+2)*d2i(m+lz+2,6))
                           m=m+3
 702                    continue
                     end if
                     go to 15
    8             continue
                     m=1
                     do 800 prim=1,nv
                        energy=energy+dens(prim,n)*
     $                      (i2(m+lx  )*i2(m+ly  )*i2(m+lz  )+
     $                       i2(m+lx+1)*i2(m+ly+1)*i2(m+lz+1)+
     $                       i2(m+lx+2)*i2(m+ly+2)*i2(m+lz+2)+
     $                       i2(m+lx+3)*i2(m+ly+3)*i2(m+lz+3))
                        der(1,1)=der(1,1)+dens(prim,n)*
     #                      (di(m+lx  ,1)*i2(m+ly  )*i2(m+lz  )+
     #                       di(m+lx+1,1)*i2(m+ly+1)*i2(m+lz+1)+
     #                       di(m+lx+2,1)*i2(m+ly+2)*i2(m+lz+2)+
     #                       di(m+lx+3,1)*i2(m+ly+3)*i2(m+lz+3))
                        der(2,1)=der(2,1)+dens(prim,n)*
     #                      (i2(m+lx  )*di(m+ly  ,1)*i2(m+lz  )+
     #                       i2(m+lx+1)*di(m+ly+1,1)*i2(m+lz+1)+
     #                       i2(m+lx+2)*di(m+ly+2,1)*i2(m+lz+2)+
     #                       i2(m+lx+3)*di(m+ly+3,1)*i2(m+lz+3))
                        der(3,1)=der(3,1)+dens(prim,n)*
     #                      (i2(m+lx  )*i2(m+ly  )*di(m+lz  ,1)+
     #                       i2(m+lx+1)*i2(m+ly+1)*di(m+lz+1,1)+
     #                       i2(m+lx+2)*i2(m+ly+2)*di(m+lz+2,1)+
     #                       i2(m+lx+3)*i2(m+ly+3)*di(m+lz+3,1))
                        d2(1)=d2(1)+dens(prim,n)*
     $                      (d2i(m+lx  ,1)*i2(m+ly  )*i2(m+lz  )+
     $                       d2i(m+lx+1,1)*i2(m+ly+1)*i2(m+lz+1)+
     $                       d2i(m+lx+2,1)*i2(m+ly+2)*i2(m+lz+2)+
     $                       d2i(m+lx+3,1)*i2(m+ly+3)*i2(m+lz+3))
                        d2(2)=d2(2)+dens(prim,n)*
     $                      (di(m+lx  ,1)*di(m+ly  ,1)*i2(m+lz  )+
     $                       di(m+lx+1,1)*di(m+ly+1,1)*i2(m+lz+1)+
     $                       di(m+lx+2,1)*di(m+ly+2,1)*i2(m+lz+2)+
     $                       di(m+lx+3,1)*di(m+ly+3,1)*i2(m+lz+3))
                        d2(3)=d2(3)+dens(prim,n)*
     $                      (i2(m+lx  )*d2i(m+ly  ,1)*i2(m+lz  )+
     $                       i2(m+lx+1)*d2i(m+ly+1,1)*i2(m+lz+1)+
     $                       i2(m+lx+2)*d2i(m+ly+2,1)*i2(m+lz+2)+
     $                       i2(m+lx+3)*d2i(m+ly+3,1)*i2(m+lz+3))
                        d2(4)=d2(4)+dens(prim,n)*
     $                      (di(m+lx  ,1)*i2(m+ly  )*di(m+lz  ,1)+
     $                       di(m+lx+1,1)*i2(m+ly+1)*di(m+lz+1,1)+
     $                       di(m+lx+2,1)*i2(m+ly+2)*di(m+lz+2,1)+
     $                       di(m+lx+3,1)*i2(m+ly+3)*di(m+lz+3,1))
                        d2(5)=d2(5)+dens(prim,n)*
     $                      (i2(m+lx  )*di(m+ly  ,1)*di(m+lz  ,1)+
     $                       i2(m+lx+1)*di(m+ly+1,1)*di(m+lz+1,1)+
     $                       i2(m+lx+2)*di(m+ly+2,1)*di(m+lz+2,1)+
     $                       i2(m+lx+3)*di(m+ly+3,1)*di(m+lz+3,1))
                        d2(6)=d2(6)+dens(prim,n)*
     $                      (i2(m+lx  )*i2(m+ly  )*d2i(m+lz  ,1)+
     $                       i2(m+lx+1)*i2(m+ly+1)*d2i(m+lz+1,1)+
     $                       i2(m+lx+2)*i2(m+ly+2)*d2i(m+lz+2,1)+
     $                       i2(m+lx+3)*i2(m+ly+3)*d2i(m+lz+3,1))
                        m=m+4
 800                 continue
                     if (npass.lt.4) then
                        m=1
                        do 801 prim=1,nv
                           der(1,2)=der(1,2)+dens(prim,n)*
     #                        (di(m+lx  ,2)*i2(m+ly  )*i2(m+lz  )+
     #                         di(m+lx+1,2)*i2(m+ly+1)*i2(m+lz+1)+
     #                         di(m+lx+2,2)*i2(m+ly+2)*i2(m+lz+2)+
     #                         di(m+lx+3,2)*i2(m+ly+3)*i2(m+lz+3))
                           der(2,2)=der(2,2)+dens(prim,n)*
     #                        (i2(m+lx  )*di(m+ly  ,2)*i2(m+lz  )+
     #                         i2(m+lx+1)*di(m+ly+1,2)*i2(m+lz+1)+
     #                         i2(m+lx+2)*di(m+ly+2,2)*i2(m+lz+2)+
     #                         i2(m+lx+3)*di(m+ly+3,2)*i2(m+lz+3))
                           der(3,2)=der(3,2)+dens(prim,n)*
     #                        (i2(m+lx  )*i2(m+ly  )*di(m+lz  ,2)+
     #                         i2(m+lx+1)*i2(m+ly+1)*di(m+lz+1,2)+
     #                         i2(m+lx+2)*i2(m+ly+2)*di(m+lz+2,2)+
     #                         i2(m+lx+3)*i2(m+ly+3)*di(m+lz+3,2))
                           d2(7)=d2(7)+dens(prim,n)*
     $                         (d2i(m+lx  ,2)*i2(m+ly  )*i2(m+lz  )+
     $                          d2i(m+lx+1,2)*i2(m+ly+1)*i2(m+lz+1)+
     $                          d2i(m+lx+2,2)*i2(m+ly+2)*i2(m+lz+2)+
     $                          d2i(m+lx+3,2)*i2(m+ly+3)*i2(m+lz+3))
                           d2(8)=d2(8)+dens(prim,n)*
     $                         (di(m+lx  ,2)*di(m+ly  ,1)*i2(m+lz  )+
     $                          di(m+lx+1,2)*di(m+ly+1,1)*i2(m+lz+1)+
     $                          di(m+lx+2,2)*di(m+ly+2,1)*i2(m+lz+2)+
     $                          di(m+lx+3,2)*di(m+ly+3,1)*i2(m+lz+3))
                           d2(9)=d2(9)+dens(prim,n)*
     $                         (di(m+lx  ,2)*i2(m+ly  )*di(m+lz  ,1)+
     $                          di(m+lx+1,2)*i2(m+ly+1)*di(m+lz+1,1)+
     $                          di(m+lx+2,2)*i2(m+ly+2)*di(m+lz+2,1)+
     $                          di(m+lx+3,2)*i2(m+ly+3)*di(m+lz+3,1))
                           d2(10)=d2(10)+dens(prim,n)*
     $                         (d2i(m+lx  ,3)*i2(m+ly  )*i2(m+lz  )+
     $                          d2i(m+lx+1,3)*i2(m+ly+1)*i2(m+lz+1)+
     $                          d2i(m+lx+2,3)*i2(m+ly+2)*i2(m+lz+2)+
     $                          d2i(m+lx+3,3)*i2(m+ly+3)*i2(m+lz+3))
                           d2(11)=d2(11)+dens(prim,n)*
     $                         (di(m+lx  ,1)*di(m+ly  ,2)*i2(m+lz  )+
     $                          di(m+lx+1,1)*di(m+ly+1,2)*i2(m+lz+1)+
     $                          di(m+lx+2,1)*di(m+ly+2,2)*i2(m+lz+2)+
     $                          di(m+lx+3,1)*di(m+ly+3,2)*i2(m+lz+3))
                           d2(12)=d2(12)+dens(prim,n)*
     $                         (i2(m+lx  )*d2i(m+ly  ,2)*i2(m+lz  )+
     $                          i2(m+lx+1)*d2i(m+ly+1,2)*i2(m+lz+1)+
     $                          i2(m+lx+2)*d2i(m+ly+2,2)*i2(m+lz+2)+
     $                          i2(m+lx+3)*d2i(m+ly+3,2)*i2(m+lz+3))
                           d2(13)=d2(13)+dens(prim,n)*
     $                         (i2(m+lx  )*di(m+ly  ,2)*di(m+lz  ,1)+
     $                          i2(m+lx+1)*di(m+ly+1,2)*di(m+lz+1,1)+
     $                          i2(m+lx+2)*di(m+ly+2,2)*di(m+lz+2,1)+
     $                          i2(m+lx+3)*di(m+ly+3,2)*di(m+lz+3,1))
                           d2(14)=d2(14)+dens(prim,n)*
     $                         (di(m+lx  ,2)*di(m+ly  ,2)*i2(m+lz  )+
     $                          di(m+lx+1,2)*di(m+ly+1,2)*i2(m+lz+1)+
     $                          di(m+lx+2,2)*di(m+ly+2,2)*i2(m+lz+2)+
     $                          di(m+lx+3,2)*di(m+ly+3,2)*i2(m+lz+3))
                           d2(15)=d2(15)+dens(prim,n)*
     $                         (i2(m+lx  )*d2i(m+ly  ,3)*i2(m+lz  )+
     $                          i2(m+lx+1)*d2i(m+ly+1,3)*i2(m+lz+1)+
     $                          i2(m+lx+2)*d2i(m+ly+2,3)*i2(m+lz+2)+
     $                          i2(m+lx+3)*d2i(m+ly+3,3)*i2(m+lz+3))
                           d2(16)=d2(16)+dens(prim,n)*
     $                         (di(m+lx  ,1)*i2(m+ly  )*di(m+lz  ,2)+
     $                          di(m+lx+1,1)*i2(m+ly+1)*di(m+lz+1,2)+
     $                          di(m+lx+2,1)*i2(m+ly+2)*di(m+lz+2,2)+
     $                          di(m+lx+3,1)*i2(m+ly+3)*di(m+lz+3,2))
                           d2(17)=d2(17)+dens(prim,n)*
     $                         (i2(m+lx  )*di(m+ly  ,1)*di(m+lz  ,2)+
     $                          i2(m+lx+1)*di(m+ly+1,1)*di(m+lz+1,2)+
     $                          i2(m+lx+2)*di(m+ly+2,1)*di(m+lz+2,2)+
     $                          i2(m+lx+3)*di(m+ly+3,1)*di(m+lz+3,2))
                           d2(18)=d2(18)+dens(prim,n)*
     $                         (i2(m+lx  )*i2(m+ly  )*d2i(m+lz  ,2)+
     $                          i2(m+lx+1)*i2(m+ly+1)*d2i(m+lz+1,2)+
     $                          i2(m+lx+2)*i2(m+ly+2)*d2i(m+lz+2,2)+
     $                          i2(m+lx+3)*i2(m+ly+3)*d2i(m+lz+3,2))
                           d2(19)=d2(19)+dens(prim,n)*
     $                         (di(m+lx  ,2)*i2(m+ly  )*di(m+lz  ,2)+
     $                          di(m+lx+1,2)*i2(m+ly+1)*di(m+lz+1,2)+
     $                          di(m+lx+2,2)*i2(m+ly+2)*di(m+lz+2,2)+
     $                          di(m+lx+3,2)*i2(m+ly+3)*di(m+lz+3,2))
                           d2(20)=d2(20)+dens(prim,n)*
     $                         (i2(m+lx  )*di(m+ly  ,2)*di(m+lz  ,2)+
     $                          i2(m+lx+1)*di(m+ly+1,2)*di(m+lz+1,2)+
     $                          i2(m+lx+2)*di(m+ly+2,2)*di(m+lz+2,2)+
     $                          i2(m+lx+3)*di(m+ly+3,2)*di(m+lz+3,2))
                           d2(21)=d2(21)+dens(prim,n)*
     $                         (i2(m+lx  )*i2(m+ly  )*d2i(m+lz  ,3)+
     $                          i2(m+lx+1)*i2(m+ly+1)*d2i(m+lz+1,3)+
     $                          i2(m+lx+2)*i2(m+ly+2)*d2i(m+lz+2,3)+
     $                          i2(m+lx+3)*i2(m+ly+3)*d2i(m+lz+3,3))
                           m=m+4
 801                    continue
                     end if
                     if (npass.eq.1) then
                        m=1
                        do 802 prim=1,nv
                           der(1,3)=der(1,3)+dens(prim,n)*
     #                        (di(m+lx  ,3)*i2(m+ly  )*i2(m+lz  )+
     #                         di(m+lx+1,3)*i2(m+ly+1)*i2(m+lz+1)+
     #                         di(m+lx+2,3)*i2(m+ly+2)*i2(m+lz+2)+
     #                         di(m+lx+3,3)*i2(m+ly+3)*i2(m+lz+3))
                           der(2,3)=der(2,3)+dens(prim,n)*
     #                        (i2(m+lx  )*di(m+ly  ,3)*i2(m+lz  )+
     #                         i2(m+lx+1)*di(m+ly+1,3)*i2(m+lz+1)+
     #                         i2(m+lx+2)*di(m+ly+2,3)*i2(m+lz+2)+
     #                         i2(m+lx+3)*di(m+ly+3,3)*i2(m+lz+3))
                           der(3,3)=der(3,3)+dens(prim,n)*
     #                        (i2(m+lx  )*i2(m+ly  )*di(m+lz  ,3)+
     #                         i2(m+lx+1)*i2(m+ly+1)*di(m+lz+1,3)+
     #                         i2(m+lx+2)*i2(m+ly+2)*di(m+lz+2,3)+
     #                         i2(m+lx+3)*i2(m+ly+3)*di(m+lz+3,3))
                           d2(22)=d2(22)+dens(prim,n)*
     $                         (d2i(m+lx  ,4)*i2(m+ly  )*i2(m+lz  )+
     $                          d2i(m+lx+1,4)*i2(m+ly+1)*i2(m+lz+1)+
     $                          d2i(m+lx+2,4)*i2(m+ly+2)*i2(m+lz+2)+
     $                          d2i(m+lx+3,4)*i2(m+ly+3)*i2(m+lz+3))
                           d2(23)=d2(23)+dens(prim,n)*
     $                         (di(m+lx  ,3)*di(m+ly  ,1)*i2(m+lz  )+
     $                          di(m+lx+1,3)*di(m+ly+1,1)*i2(m+lz+1)+
     $                          di(m+lx+2,3)*di(m+ly+2,1)*i2(m+lz+2)+
     $                          di(m+lx+3,3)*di(m+ly+3,1)*i2(m+lz+3))
                           d2(24)=d2(24)+dens(prim,n)*
     $                         (di(m+lx  ,3)*i2(m+ly  )*di(m+lz  ,1)+
     $                          di(m+lx+1,3)*i2(m+ly+1)*di(m+lz+1,1)+
     $                          di(m+lx+2,3)*i2(m+ly+2)*di(m+lz+2,1)+
     $                          di(m+lx+3,3)*i2(m+ly+3)*di(m+lz+3,1))
                           d2(25)=d2(25)+dens(prim,n)*
     $                         (d2i(m+lx  ,5)*i2(m+ly  )*i2(m+lz  )+
     $                          d2i(m+lx+1,5)*i2(m+ly+1)*i2(m+lz+1)+
     $                          d2i(m+lx+2,5)*i2(m+ly+2)*i2(m+lz+2)+
     $                          d2i(m+lx+3,5)*i2(m+ly+3)*i2(m+lz+3))
                           d2(26)=d2(26)+dens(prim,n)*
     $                         (di(m+lx  ,3)*di(m+ly  ,2)*i2(m+lz  )+
     $                          di(m+lx+1,3)*di(m+ly+1,2)*i2(m+lz+1)+
     $                          di(m+lx+2,3)*di(m+ly+2,2)*i2(m+lz+2)+
     $                          di(m+lx+3,3)*di(m+ly+3,2)*i2(m+lz+3))
                           d2(27)=d2(27)+dens(prim,n)*
     $                         (di(m+lx  ,3)*i2(m+ly  )*di(m+lz  ,2)+
     $                          di(m+lx+1,3)*i2(m+ly+1)*di(m+lz+1,2)+
     $                          di(m+lx+2,3)*i2(m+ly+2)*di(m+lz+2,2)+
     $                          di(m+lx+3,3)*i2(m+ly+3)*di(m+lz+3,2))
                           d2(28)=d2(28)+dens(prim,n)*
     $                         (d2i(m+lx  ,6)*i2(m+ly  )*i2(m+lz  )+
     $                          d2i(m+lx+1,6)*i2(m+ly+1)*i2(m+lz+1)+
     $                          d2i(m+lx+2,6)*i2(m+ly+2)*i2(m+lz+2)+
     $                          d2i(m+lx+3,6)*i2(m+ly+3)*i2(m+lz+3))
                           d2(29)=d2(29)+dens(prim,n)*
     $                         (di(m+lx  ,1)*di(m+ly  ,3)*i2(m+lz  )+
     $                          di(m+lx+1,1)*di(m+ly+1,3)*i2(m+lz+1)+
     $                          di(m+lx+2,1)*di(m+ly+2,3)*i2(m+lz+2)+
     $                          di(m+lx+3,1)*di(m+ly+3,3)*i2(m+lz+3))
                           d2(30)=d2(30)+dens(prim,n)*
     $                         (i2(m+lx  )*d2i(m+ly  ,4)*i2(m+lz  )+
     $                          i2(m+lx+1)*d2i(m+ly+1,4)*i2(m+lz+1)+
     $                          i2(m+lx+2)*d2i(m+ly+2,4)*i2(m+lz+2)+
     $                          i2(m+lx+3)*d2i(m+ly+3,4)*i2(m+lz+3))
                           d2(31)=d2(31)+dens(prim,n)*
     $                         (i2(m+lx  )*di(m+ly  ,3)*di(m+lz  ,1)+
     $                          i2(m+lx+1)*di(m+ly+1,3)*di(m+lz+1,1)+
     $                          i2(m+lx+2)*di(m+ly+2,3)*di(m+lz+2,1)+
     $                          i2(m+lx+3)*di(m+ly+3,3)*di(m+lz+3,1))
                           d2(32)=d2(32)+dens(prim,n)*
     $                         (di(m+lx  ,2)*di(m+ly  ,3)*i2(m+lz  )+
     $                          di(m+lx+1,2)*di(m+ly+1,3)*i2(m+lz+1)+
     $                          di(m+lx+2,2)*di(m+ly+2,3)*i2(m+lz+2)+
     $                          di(m+lx+3,2)*di(m+ly+3,3)*i2(m+lz+3))
                           d2(33)=d2(33)+dens(prim,n)*
     $                         (i2(m+lx  )*d2i(m+ly  ,5)*i2(m+lz  )+
     $                          i2(m+lx+1)*d2i(m+ly+1,5)*i2(m+lz+1)+
     $                          i2(m+lx+2)*d2i(m+ly+2,5)*i2(m+lz+2)+
     $                          i2(m+lx+3)*d2i(m+ly+3,5)*i2(m+lz+3))
                           d2(34)=d2(34)+dens(prim,n)*
     $                         (i2(m+lx  )*di(m+ly  ,3)*di(m+lz  ,2)+
     $                          i2(m+lx+1)*di(m+ly+1,3)*di(m+lz+1,2)+
     $                          i2(m+lx+2)*di(m+ly+2,3)*di(m+lz+2,2)+
     $                          i2(m+lx+3)*di(m+ly+3,3)*di(m+lz+3,2))
                           d2(35)=d2(35)+dens(prim,n)*
     $                         (di(m+lx  ,3)*di(m+ly  ,3)*i2(m+lz  )+
     $                          di(m+lx+1,3)*di(m+ly+1,3)*i2(m+lz+1)+
     $                          di(m+lx+2,3)*di(m+ly+2,3)*i2(m+lz+2)+
     $                          di(m+lx+3,3)*di(m+ly+3,3)*i2(m+lz+3))
                           d2(36)=d2(36)+dens(prim,n)*
     $                         (i2(m+lx  )*d2i(m+ly  ,6)*i2(m+lz  )+
     $                          i2(m+lx+1)*d2i(m+ly+1,6)*i2(m+lz+1)+
     $                          i2(m+lx+2)*d2i(m+ly+2,6)*i2(m+lz+2)+
     $                          i2(m+lx+3)*d2i(m+ly+3,6)*i2(m+lz+3))
                           d2(37)=d2(37)+dens(prim,n)*
     $                         (di(m+lx  ,1)*i2(m+ly  )*di(m+lz  ,3)+
     $                          di(m+lx+1,1)*i2(m+ly+1)*di(m+lz+1,3)+
     $                          di(m+lx+2,1)*i2(m+ly+2)*di(m+lz+2,3)+
     $                          di(m+lx+3,1)*i2(m+ly+3)*di(m+lz+3,3))
                           d2(38)=d2(38)+dens(prim,n)*
     $                         (i2(m+lx  )*di(m+ly  ,1)*di(m+lz  ,3)+
     $                          i2(m+lx+1)*di(m+ly+1,1)*di(m+lz+1,3)+
     $                          i2(m+lx+2)*di(m+ly+2,1)*di(m+lz+2,3)+
     $                          i2(m+lx+3)*di(m+ly+3,1)*di(m+lz+3,3))
                           d2(39)=d2(39)+dens(prim,n)*
     $                         (i2(m+lx  )*i2(m+ly  )*d2i(m+lz  ,4)+
     $                          i2(m+lx+1)*i2(m+ly+1)*d2i(m+lz+1,4)+
     $                          i2(m+lx+2)*i2(m+ly+2)*d2i(m+lz+2,4)+
     $                          i2(m+lx+3)*i2(m+ly+3)*d2i(m+lz+3,4))
                           d2(40)=d2(40)+dens(prim,n)*
     $                         (di(m+lx  ,2)*i2(m+ly  )*di(m+lz  ,3)+
     $                          di(m+lx+1,2)*i2(m+ly+1)*di(m+lz+1,3)+
     $                          di(m+lx+2,2)*i2(m+ly+2)*di(m+lz+2,3)+
     $                          di(m+lx+3,2)*i2(m+ly+3)*di(m+lz+3,3))
                           d2(41)=d2(41)+dens(prim,n)*
     $                         (i2(m+lx  )*di(m+ly  ,2)*di(m+lz  ,3)+
     $                          i2(m+lx+1)*di(m+ly+1,2)*di(m+lz+1,3)+
     $                          i2(m+lx+2)*di(m+ly+2,2)*di(m+lz+2,3)+
     $                          i2(m+lx+3)*di(m+ly+3,2)*di(m+lz+3,3))
                           d2(42)=d2(42)+dens(prim,n)*
     $                         (i2(m+lx  )*i2(m+ly  )*d2i(m+lz  ,5)+
     $                          i2(m+lx+1)*i2(m+ly+1)*d2i(m+lz+1,5)+
     $                          i2(m+lx+2)*i2(m+ly+2)*d2i(m+lz+2,5)+
     $                          i2(m+lx+3)*i2(m+ly+3)*d2i(m+lz+3,5))
                           d2(43)=d2(43)+dens(prim,n)*
     $                         (di(m+lx  ,3)*i2(m+ly  )*di(m+lz  ,3)+
     $                          di(m+lx+1,3)*i2(m+ly+1)*di(m+lz+1,3)+
     $                          di(m+lx+2,3)*i2(m+ly+2)*di(m+lz+2,3)+
     $                          di(m+lx+3,3)*i2(m+ly+3)*di(m+lz+3,3))
                           d2(44)=d2(44)+dens(prim,n)*
     $                         (i2(m+lx  )*di(m+ly  ,3)*di(m+lz  ,3)+
     $                          i2(m+lx+1)*di(m+ly+1,3)*di(m+lz+1,3)+
     $                          i2(m+lx+2)*di(m+ly+2,3)*di(m+lz+2,3)+
     $                          i2(m+lx+3)*di(m+ly+3,3)*di(m+lz+3,3))
                           d2(45)=d2(45)+dens(prim,n)*
     $                         (i2(m+lx  )*i2(m+ly  )*d2i(m+lz  ,6)+
     $                          i2(m+lx+1)*i2(m+ly+1)*d2i(m+lz+1,6)+
     $                          i2(m+lx+2)*i2(m+ly+2)*d2i(m+lz+2,6)+
     $                          i2(m+lx+3)*i2(m+ly+3)*d2i(m+lz+3,6))
                           m=m+4
 802                    continue
                     end if
                     go to 15
    9             continue
                     m=1
                     do 900 prim=1,nv
                        energy=energy+dens(prim,n)*
     $                      (i2(m+lx  )*i2(m+ly  )*i2(m+lz  )+
     $                       i2(m+lx+1)*i2(m+ly+1)*i2(m+lz+1)+
     $                       i2(m+lx+2)*i2(m+ly+2)*i2(m+lz+2)+
     $                       i2(m+lx+3)*i2(m+ly+3)*i2(m+lz+3)+
     $                       i2(m+lx+4)*i2(m+ly+4)*i2(m+lz+4))
                        der(1,1)=der(1,1)+dens(prim,n)*
     #                      (di(m+lx  ,1)*i2(m+ly  )*i2(m+lz  )+
     #                       di(m+lx+1,1)*i2(m+ly+1)*i2(m+lz+1)+
     #                       di(m+lx+2,1)*i2(m+ly+2)*i2(m+lz+2)+
     #                       di(m+lx+3,1)*i2(m+ly+3)*i2(m+lz+3)+
     #                       di(m+lx+4,1)*i2(m+ly+4)*i2(m+lz+4))
                        der(2,1)=der(2,1)+dens(prim,n)*
     #                      (i2(m+lx  )*di(m+ly  ,1)*i2(m+lz  )+
     #                       i2(m+lx+1)*di(m+ly+1,1)*i2(m+lz+1)+
     #                       i2(m+lx+2)*di(m+ly+2,1)*i2(m+lz+2)+
     #                       i2(m+lx+3)*di(m+ly+3,1)*i2(m+lz+3)+
     #                       i2(m+lx+4)*di(m+ly+4,1)*i2(m+lz+4))
                        der(3,1)=der(3,1)+dens(prim,n)*
     #                      (i2(m+lx  )*i2(m+ly  )*di(m+lz  ,1)+
     #                       i2(m+lx+1)*i2(m+ly+1)*di(m+lz+1,1)+
     #                       i2(m+lx+2)*i2(m+ly+2)*di(m+lz+2,1)+
     #                       i2(m+lx+3)*i2(m+ly+3)*di(m+lz+3,1)+
     #                       i2(m+lx+4)*i2(m+ly+4)*di(m+lz+4,1))
                        d2(1)=d2(1)+dens(prim,n)*
     $                      (d2i(m+lx  ,1)*i2(m+ly  )*i2(m+lz  )+
     $                       d2i(m+lx+1,1)*i2(m+ly+1)*i2(m+lz+1)+
     $                       d2i(m+lx+2,1)*i2(m+ly+2)*i2(m+lz+2)+
     $                       d2i(m+lx+3,1)*i2(m+ly+3)*i2(m+lz+3)+
     $                       d2i(m+lx+4,1)*i2(m+ly+4)*i2(m+lz+4))
                        d2(2)=d2(2)+dens(prim,n)*
     $                      (di(m+lx  ,1)*di(m+ly  ,1)*i2(m+lz  )+
     $                       di(m+lx+1,1)*di(m+ly+1,1)*i2(m+lz+1)+
     $                       di(m+lx+2,1)*di(m+ly+2,1)*i2(m+lz+2)+
     $                       di(m+lx+3,1)*di(m+ly+3,1)*i2(m+lz+3)+
     $                       di(m+lx+4,1)*di(m+ly+4,1)*i2(m+lz+4))
                        d2(3)=d2(3)+dens(prim,n)*
     $                      (i2(m+lx  )*d2i(m+ly  ,1)*i2(m+lz  )+
     $                       i2(m+lx+1)*d2i(m+ly+1,1)*i2(m+lz+1)+
     $                       i2(m+lx+2)*d2i(m+ly+2,1)*i2(m+lz+2)+
     $                       i2(m+lx+3)*d2i(m+ly+3,1)*i2(m+lz+3)+
     $                       i2(m+lx+4)*d2i(m+ly+4,1)*i2(m+lz+4))
                        d2(4)=d2(4)+dens(prim,n)*
     $                      (di(m+lx  ,1)*i2(m+ly  )*di(m+lz  ,1)+
     $                       di(m+lx+1,1)*i2(m+ly+1)*di(m+lz+1,1)+
     $                       di(m+lx+2,1)*i2(m+ly+2)*di(m+lz+2,1)+
     $                       di(m+lx+3,1)*i2(m+ly+3)*di(m+lz+3,1)+
     $                       di(m+lx+4,1)*i2(m+ly+4)*di(m+lz+4,1))
                        d2(5)=d2(5)+dens(prim,n)*
     $                      (i2(m+lx  )*di(m+ly  ,1)*di(m+lz  ,1)+
     $                       i2(m+lx+1)*di(m+ly+1,1)*di(m+lz+1,1)+
     $                       i2(m+lx+2)*di(m+ly+2,1)*di(m+lz+2,1)+
     $                       i2(m+lx+3)*di(m+ly+3,1)*di(m+lz+3,1)+
     $                       i2(m+lx+4)*di(m+ly+4,1)*di(m+lz+4,1))
                        d2(6)=d2(6)+dens(prim,n)*
     $                      (i2(m+lx  )*i2(m+ly  )*d2i(m+lz  ,1)+
     $                       i2(m+lx+1)*i2(m+ly+1)*d2i(m+lz+1,1)+
     $                       i2(m+lx+2)*i2(m+ly+2)*d2i(m+lz+2,1)+
     $                       i2(m+lx+3)*i2(m+ly+3)*d2i(m+lz+3,1)+
     $                       i2(m+lx+4)*i2(m+ly+4)*d2i(m+lz+4,1))
                        m=m+5
 900                 continue
                     if (npass.lt.4) then
                        m=1
                        do 901 prim=1,nv
                           der(1,2)=der(1,2)+dens(prim,n)*
     #                        (di(m+lx  ,2)*i2(m+ly  )*i2(m+lz  )+
     #                         di(m+lx+1,2)*i2(m+ly+1)*i2(m+lz+1)+
     #                         di(m+lx+2,2)*i2(m+ly+2)*i2(m+lz+2)+
     #                         di(m+lx+3,2)*i2(m+ly+3)*i2(m+lz+3)+
     #                         di(m+lx+4,2)*i2(m+ly+4)*i2(m+lz+4))
                           der(2,2)=der(2,2)+dens(prim,n)*
     #                        (i2(m+lx  )*di(m+ly  ,2)*i2(m+lz  )+
     #                         i2(m+lx+1)*di(m+ly+1,2)*i2(m+lz+1)+
     #                         i2(m+lx+2)*di(m+ly+2,2)*i2(m+lz+2)+
     #                         i2(m+lx+3)*di(m+ly+3,2)*i2(m+lz+3)+
     #                         i2(m+lx+4)*di(m+ly+4,2)*i2(m+lz+4))
                           der(3,2)=der(3,2)+dens(prim,n)*
     #                        (i2(m+lx  )*i2(m+ly  )*di(m+lz  ,2)+
     #                         i2(m+lx+1)*i2(m+ly+1)*di(m+lz+1,2)+
     #                         i2(m+lx+2)*i2(m+ly+2)*di(m+lz+2,2)+
     #                         i2(m+lx+3)*i2(m+ly+3)*di(m+lz+3,2)+
     #                         i2(m+lx+4)*i2(m+ly+4)*di(m+lz+4,2))
                           d2(7)=d2(7)+dens(prim,n)*
     $                         (d2i(m+lx  ,2)*i2(m+ly  )*i2(m+lz  )+
     $                          d2i(m+lx+1,2)*i2(m+ly+1)*i2(m+lz+1)+
     $                          d2i(m+lx+2,2)*i2(m+ly+2)*i2(m+lz+2)+
     $                          d2i(m+lx+3,2)*i2(m+ly+3)*i2(m+lz+3)+
     $                          d2i(m+lx+4,2)*i2(m+ly+4)*i2(m+lz+4))
                           d2(8)=d2(8)+dens(prim,n)*
     $                         (di(m+lx  ,2)*di(m+ly  ,1)*i2(m+lz  )+
     $                          di(m+lx+1,2)*di(m+ly+1,1)*i2(m+lz+1)+
     $                          di(m+lx+2,2)*di(m+ly+2,1)*i2(m+lz+2)+
     $                          di(m+lx+3,2)*di(m+ly+3,1)*i2(m+lz+3)+
     $                          di(m+lx+4,2)*di(m+ly+4,1)*i2(m+lz+4))
                           d2(9)=d2(9)+dens(prim,n)*
     $                         (di(m+lx  ,2)*i2(m+ly  )*di(m+lz  ,1)+
     $                          di(m+lx+1,2)*i2(m+ly+1)*di(m+lz+1,1)+
     $                          di(m+lx+2,2)*i2(m+ly+2)*di(m+lz+2,1)+
     $                          di(m+lx+3,2)*i2(m+ly+3)*di(m+lz+3,1)+
     $                          di(m+lx+4,2)*i2(m+ly+4)*di(m+lz+4,1))
                           d2(10)=d2(10)+dens(prim,n)*
     $                         (d2i(m+lx  ,3)*i2(m+ly  )*i2(m+lz  )+
     $                          d2i(m+lx+1,3)*i2(m+ly+1)*i2(m+lz+1)+
     $                          d2i(m+lx+2,3)*i2(m+ly+2)*i2(m+lz+2)+
     $                          d2i(m+lx+3,3)*i2(m+ly+3)*i2(m+lz+3)+
     $                          d2i(m+lx+4,3)*i2(m+ly+4)*i2(m+lz+4))
                           d2(11)=d2(11)+dens(prim,n)*
     $                         (di(m+lx  ,1)*di(m+ly  ,2)*i2(m+lz  )+
     $                          di(m+lx+1,1)*di(m+ly+1,2)*i2(m+lz+1)+
     $                          di(m+lx+2,1)*di(m+ly+2,2)*i2(m+lz+2)+
     $                          di(m+lx+3,1)*di(m+ly+3,2)*i2(m+lz+3)+
     $                          di(m+lx+4,1)*di(m+ly+4,2)*i2(m+lz+4))
                           d2(12)=d2(12)+dens(prim,n)*
     $                         (i2(m+lx  )*d2i(m+ly  ,2)*i2(m+lz  )+
     $                          i2(m+lx+1)*d2i(m+ly+1,2)*i2(m+lz+1)+
     $                          i2(m+lx+2)*d2i(m+ly+2,2)*i2(m+lz+2)+
     $                          i2(m+lx+3)*d2i(m+ly+3,2)*i2(m+lz+3)+
     $                          i2(m+lx+4)*d2i(m+ly+4,2)*i2(m+lz+4))
                           d2(13)=d2(13)+dens(prim,n)*
     $                         (i2(m+lx  )*di(m+ly  ,2)*di(m+lz  ,1)+
     $                          i2(m+lx+1)*di(m+ly+1,2)*di(m+lz+1,1)+
     $                          i2(m+lx+2)*di(m+ly+2,2)*di(m+lz+2,1)+
     $                          i2(m+lx+3)*di(m+ly+3,2)*di(m+lz+3,1)+
     $                          i2(m+lx+4)*di(m+ly+4,2)*di(m+lz+4,1))
                           d2(14)=d2(14)+dens(prim,n)*
     $                         (di(m+lx  ,2)*di(m+ly  ,2)*i2(m+lz  )+
     $                          di(m+lx+1,2)*di(m+ly+1,2)*i2(m+lz+1)+
     $                          di(m+lx+2,2)*di(m+ly+2,2)*i2(m+lz+2)+
     $                          di(m+lx+3,2)*di(m+ly+3,2)*i2(m+lz+3)+
     $                          di(m+lx+4,2)*di(m+ly+4,2)*i2(m+lz+4))
                           d2(15)=d2(15)+dens(prim,n)*
     $                         (i2(m+lx  )*d2i(m+ly  ,3)*i2(m+lz  )+
     $                          i2(m+lx+1)*d2i(m+ly+1,3)*i2(m+lz+1)+
     $                          i2(m+lx+2)*d2i(m+ly+2,3)*i2(m+lz+2)+
     $                          i2(m+lx+3)*d2i(m+ly+3,3)*i2(m+lz+3)+
     $                          i2(m+lx+4)*d2i(m+ly+4,3)*i2(m+lz+4))
                           d2(16)=d2(16)+dens(prim,n)*
     $                         (di(m+lx  ,1)*i2(m+ly  )*di(m+lz  ,2)+
     $                          di(m+lx+1,1)*i2(m+ly+1)*di(m+lz+1,2)+
     $                          di(m+lx+2,1)*i2(m+ly+2)*di(m+lz+2,2)+
     $                          di(m+lx+3,1)*i2(m+ly+3)*di(m+lz+3,2)+
     $                          di(m+lx+4,1)*i2(m+ly+4)*di(m+lz+4,2))
                           d2(17)=d2(17)+dens(prim,n)*
     $                         (i2(m+lx  )*di(m+ly  ,1)*di(m+lz  ,2)+
     $                          i2(m+lx+1)*di(m+ly+1,1)*di(m+lz+1,2)+
     $                          i2(m+lx+2)*di(m+ly+2,1)*di(m+lz+2,2)+
     $                          i2(m+lx+3)*di(m+ly+3,1)*di(m+lz+3,2)+
     $                          i2(m+lx+4)*di(m+ly+4,1)*di(m+lz+4,2))
                           d2(18)=d2(18)+dens(prim,n)*
     $                         (i2(m+lx  )*i2(m+ly  )*d2i(m+lz  ,2)+
     $                          i2(m+lx+1)*i2(m+ly+1)*d2i(m+lz+1,2)+
     $                          i2(m+lx+2)*i2(m+ly+2)*d2i(m+lz+2,2)+
     $                          i2(m+lx+3)*i2(m+ly+3)*d2i(m+lz+3,2)+
     $                          i2(m+lx+4)*i2(m+ly+4)*d2i(m+lz+4,2))
                           d2(19)=d2(19)+dens(prim,n)*
     $                         (di(m+lx  ,2)*i2(m+ly  )*di(m+lz  ,2)+
     $                          di(m+lx+1,2)*i2(m+ly+1)*di(m+lz+1,2)+
     $                          di(m+lx+2,2)*i2(m+ly+2)*di(m+lz+2,2)+
     $                          di(m+lx+3,2)*i2(m+ly+3)*di(m+lz+3,2)+
     $                          di(m+lx+4,2)*i2(m+ly+4)*di(m+lz+4,2))
                           d2(20)=d2(20)+dens(prim,n)*
     $                         (i2(m+lx  )*di(m+ly  ,2)*di(m+lz  ,2)+
     $                          i2(m+lx+1)*di(m+ly+1,2)*di(m+lz+1,2)+
     $                          i2(m+lx+2)*di(m+ly+2,2)*di(m+lz+2,2)+
     $                          i2(m+lx+3)*di(m+ly+3,2)*di(m+lz+3,2)+
     $                          i2(m+lx+4)*di(m+ly+4,2)*di(m+lz+4,2))
                           d2(21)=d2(21)+dens(prim,n)*
     $                         (i2(m+lx  )*i2(m+ly  )*d2i(m+lz  ,3)+
     $                          i2(m+lx+1)*i2(m+ly+1)*d2i(m+lz+1,3)+
     $                          i2(m+lx+2)*i2(m+ly+2)*d2i(m+lz+2,3)+
     $                          i2(m+lx+3)*i2(m+ly+3)*d2i(m+lz+3,3)+
     $                          i2(m+lx+4)*i2(m+ly+4)*d2i(m+lz+4,3))
                           m=m+5
 901                    continue
                     end if
                     if (npass.eq.1) then
                        m=1
                        do 902 prim=1,nv
                           der(1,3)=der(1,3)+dens(prim,n)*
     #                        (di(m+lx  ,3)*i2(m+ly  )*i2(m+lz  )+
     #                         di(m+lx+1,3)*i2(m+ly+1)*i2(m+lz+1)+
     #                         di(m+lx+2,3)*i2(m+ly+2)*i2(m+lz+2)+
     #                         di(m+lx+3,3)*i2(m+ly+3)*i2(m+lz+3)+
     #                         di(m+lx+4,3)*i2(m+ly+4)*i2(m+lz+4))
                           der(2,3)=der(2,3)+dens(prim,n)*
     #                        (i2(m+lx  )*di(m+ly  ,3)*i2(m+lz  )+
     #                         i2(m+lx+1)*di(m+ly+1,3)*i2(m+lz+1)+
     #                         i2(m+lx+2)*di(m+ly+2,3)*i2(m+lz+2)+
     #                         i2(m+lx+3)*di(m+ly+3,3)*i2(m+lz+3)+
     #                         i2(m+lx+4)*di(m+ly+4,3)*i2(m+lz+4))
                           der(3,3)=der(3,3)+dens(prim,n)*
     #                        (i2(m+lx  )*i2(m+ly  )*di(m+lz  ,3)+
     #                         i2(m+lx+1)*i2(m+ly+1)*di(m+lz+1,3)+
     #                         i2(m+lx+2)*i2(m+ly+2)*di(m+lz+2,3)+
     #                         i2(m+lx+3)*i2(m+ly+3)*di(m+lz+3,3)+
     #                         i2(m+lx+4)*i2(m+ly+4)*di(m+lz+4,3))
                           d2(22)=d2(22)+dens(prim,n)*
     $                         (d2i(m+lx  ,4)*i2(m+ly  )*i2(m+lz  )+
     $                          d2i(m+lx+1,4)*i2(m+ly+1)*i2(m+lz+1)+
     $                          d2i(m+lx+2,4)*i2(m+ly+2)*i2(m+lz+2)+
     $                          d2i(m+lx+3,4)*i2(m+ly+3)*i2(m+lz+3)+
     $                          d2i(m+lx+4,4)*i2(m+ly+4)*i2(m+lz+4))
                           d2(23)=d2(23)+dens(prim,n)*
     $                         (di(m+lx  ,3)*di(m+ly  ,1)*i2(m+lz  )+
     $                          di(m+lx+1,3)*di(m+ly+1,1)*i2(m+lz+1)+
     $                          di(m+lx+2,3)*di(m+ly+2,1)*i2(m+lz+2)+
     $                          di(m+lx+3,3)*di(m+ly+3,1)*i2(m+lz+3)+
     $                          di(m+lx+4,3)*di(m+ly+4,1)*i2(m+lz+4))
                           d2(24)=d2(24)+dens(prim,n)*
     $                         (di(m+lx  ,3)*i2(m+ly  )*di(m+lz  ,1)+
     $                          di(m+lx+1,3)*i2(m+ly+1)*di(m+lz+1,1)+
     $                          di(m+lx+2,3)*i2(m+ly+2)*di(m+lz+2,1)+
     $                          di(m+lx+3,3)*i2(m+ly+3)*di(m+lz+3,1)+
     $                          di(m+lx+4,3)*i2(m+ly+4)*di(m+lz+4,1))
                           d2(25)=d2(25)+dens(prim,n)*
     $                         (d2i(m+lx  ,5)*i2(m+ly  )*i2(m+lz  )+
     $                          d2i(m+lx+1,5)*i2(m+ly+1)*i2(m+lz+1)+
     $                          d2i(m+lx+2,5)*i2(m+ly+2)*i2(m+lz+2)+
     $                          d2i(m+lx+3,5)*i2(m+ly+3)*i2(m+lz+3)+
     $                          d2i(m+lx+4,5)*i2(m+ly+4)*i2(m+lz+4))
                           d2(26)=d2(26)+dens(prim,n)*
     $                         (di(m+lx  ,3)*di(m+ly  ,2)*i2(m+lz  )+
     $                          di(m+lx+1,3)*di(m+ly+1,2)*i2(m+lz+1)+
     $                          di(m+lx+2,3)*di(m+ly+2,2)*i2(m+lz+2)+
     $                          di(m+lx+3,3)*di(m+ly+3,2)*i2(m+lz+3)+
     $                          di(m+lx+4,3)*di(m+ly+4,2)*i2(m+lz+4))
                           d2(27)=d2(27)+dens(prim,n)*
     $                         (di(m+lx  ,3)*i2(m+ly  )*di(m+lz  ,2)+
     $                          di(m+lx+1,3)*i2(m+ly+1)*di(m+lz+1,2)+
     $                          di(m+lx+2,3)*i2(m+ly+2)*di(m+lz+2,2)+
     $                          di(m+lx+3,3)*i2(m+ly+3)*di(m+lz+3,2)+
     $                          di(m+lx+4,3)*i2(m+ly+4)*di(m+lz+4,2))
                           d2(28)=d2(28)+dens(prim,n)*
     $                         (d2i(m+lx  ,6)*i2(m+ly  )*i2(m+lz  )+
     $                          d2i(m+lx+1,6)*i2(m+ly+1)*i2(m+lz+1)+
     $                          d2i(m+lx+2,6)*i2(m+ly+2)*i2(m+lz+2)+
     $                          d2i(m+lx+3,6)*i2(m+ly+3)*i2(m+lz+3)+
     $                          d2i(m+lx+4,6)*i2(m+ly+4)*i2(m+lz+4))
                           d2(29)=d2(29)+dens(prim,n)*
     $                         (di(m+lx  ,1)*di(m+ly  ,3)*i2(m+lz  )+
     $                          di(m+lx+1,1)*di(m+ly+1,3)*i2(m+lz+1)+
     $                          di(m+lx+2,1)*di(m+ly+2,3)*i2(m+lz+2)+
     $                          di(m+lx+3,1)*di(m+ly+3,3)*i2(m+lz+3)+
     $                          di(m+lx+4,1)*di(m+ly+4,3)*i2(m+lz+4))
                           d2(30)=d2(30)+dens(prim,n)*
     $                         (i2(m+lx  )*d2i(m+ly  ,4)*i2(m+lz  )+
     $                          i2(m+lx+1)*d2i(m+ly+1,4)*i2(m+lz+1)+
     $                          i2(m+lx+2)*d2i(m+ly+2,4)*i2(m+lz+2)+
     $                          i2(m+lx+3)*d2i(m+ly+3,4)*i2(m+lz+3)+
     $                          i2(m+lx+4)*d2i(m+ly+4,4)*i2(m+lz+4))
                           d2(31)=d2(31)+dens(prim,n)*
     $                         (i2(m+lx  )*di(m+ly  ,3)*di(m+lz  ,1)+
     $                          i2(m+lx+1)*di(m+ly+1,3)*di(m+lz+1,1)+
     $                          i2(m+lx+2)*di(m+ly+2,3)*di(m+lz+2,1)+
     $                          i2(m+lx+3)*di(m+ly+3,3)*di(m+lz+3,1)+
     $                          i2(m+lx+4)*di(m+ly+4,3)*di(m+lz+4,1))
                           d2(32)=d2(32)+dens(prim,n)*
     $                         (di(m+lx  ,2)*di(m+ly  ,3)*i2(m+lz  )+
     $                          di(m+lx+1,2)*di(m+ly+1,3)*i2(m+lz+1)+
     $                          di(m+lx+2,2)*di(m+ly+2,3)*i2(m+lz+2)+
     $                          di(m+lx+3,2)*di(m+ly+3,3)*i2(m+lz+3)+
     $                          di(m+lx+4,2)*di(m+ly+4,3)*i2(m+lz+4))
                           d2(33)=d2(33)+dens(prim,n)*
     $                         (i2(m+lx  )*d2i(m+ly  ,5)*i2(m+lz  )+
     $                          i2(m+lx+1)*d2i(m+ly+1,5)*i2(m+lz+1)+
     $                          i2(m+lx+2)*d2i(m+ly+2,5)*i2(m+lz+2)+
     $                          i2(m+lx+3)*d2i(m+ly+3,5)*i2(m+lz+3)+
     $                          i2(m+lx+4)*d2i(m+ly+4,5)*i2(m+lz+4))
                           d2(34)=d2(34)+dens(prim,n)*
     $                         (i2(m+lx  )*di(m+ly  ,3)*di(m+lz  ,2)+
     $                          i2(m+lx+1)*di(m+ly+1,3)*di(m+lz+1,2)+
     $                          i2(m+lx+2)*di(m+ly+2,3)*di(m+lz+2,2)+
     $                          i2(m+lx+3)*di(m+ly+3,3)*di(m+lz+3,2)+
     $                          i2(m+lx+4)*di(m+ly+4,3)*di(m+lz+4,2))
                           d2(35)=d2(35)+dens(prim,n)*
     $                         (di(m+lx  ,3)*di(m+ly  ,3)*i2(m+lz  )+
     $                          di(m+lx+1,3)*di(m+ly+1,3)*i2(m+lz+1)+
     $                          di(m+lx+2,3)*di(m+ly+2,3)*i2(m+lz+2)+
     $                          di(m+lx+3,3)*di(m+ly+3,3)*i2(m+lz+3)+
     $                          di(m+lx+4,3)*di(m+ly+4,3)*i2(m+lz+4))
                           d2(36)=d2(36)+dens(prim,n)*
     $                         (i2(m+lx  )*d2i(m+ly  ,6)*i2(m+lz  )+
     $                          i2(m+lx+1)*d2i(m+ly+1,6)*i2(m+lz+1)+
     $                          i2(m+lx+2)*d2i(m+ly+2,6)*i2(m+lz+2)+
     $                          i2(m+lx+3)*d2i(m+ly+3,6)*i2(m+lz+3)+
     $                          i2(m+lx+4)*d2i(m+ly+4,6)*i2(m+lz+4))
                           d2(37)=d2(37)+dens(prim,n)*
     $                         (di(m+lx  ,1)*i2(m+ly  )*di(m+lz  ,3)+
     $                          di(m+lx+1,1)*i2(m+ly+1)*di(m+lz+1,3)+
     $                          di(m+lx+2,1)*i2(m+ly+2)*di(m+lz+2,3)+
     $                          di(m+lx+3,1)*i2(m+ly+3)*di(m+lz+3,3)+
     $                          di(m+lx+4,1)*i2(m+ly+4)*di(m+lz+4,3))
                           d2(38)=d2(38)+dens(prim,n)*
     $                         (i2(m+lx  )*di(m+ly  ,1)*di(m+lz  ,3)+
     $                          i2(m+lx+1)*di(m+ly+1,1)*di(m+lz+1,3)+
     $                          i2(m+lx+2)*di(m+ly+2,1)*di(m+lz+2,3)+
     $                          i2(m+lx+3)*di(m+ly+3,1)*di(m+lz+3,3)+
     $                          i2(m+lx+4)*di(m+ly+4,1)*di(m+lz+4,3))
                           d2(39)=d2(39)+dens(prim,n)*
     $                         (i2(m+lx  )*i2(m+ly  )*d2i(m+lz  ,4)+
     $                          i2(m+lx+1)*i2(m+ly+1)*d2i(m+lz+1,4)+
     $                          i2(m+lx+2)*i2(m+ly+2)*d2i(m+lz+2,4)+
     $                          i2(m+lx+3)*i2(m+ly+3)*d2i(m+lz+3,4)+
     $                          i2(m+lx+4)*i2(m+ly+4)*d2i(m+lz+4,4))
                           d2(40)=d2(40)+dens(prim,n)*
     $                         (di(m+lx  ,2)*i2(m+ly  )*di(m+lz  ,3)+
     $                          di(m+lx+1,2)*i2(m+ly+1)*di(m+lz+1,3)+
     $                          di(m+lx+2,2)*i2(m+ly+2)*di(m+lz+2,3)+
     $                          di(m+lx+3,2)*i2(m+ly+3)*di(m+lz+3,3)+
     $                          di(m+lx+4,2)*i2(m+ly+4)*di(m+lz+4,3))
                           d2(41)=d2(41)+dens(prim,n)*
     $                         (i2(m+lx  )*di(m+ly  ,2)*di(m+lz  ,3)+
     $                          i2(m+lx+1)*di(m+ly+1,2)*di(m+lz+1,3)+
     $                          i2(m+lx+2)*di(m+ly+2,2)*di(m+lz+2,3)+
     $                          i2(m+lx+3)*di(m+ly+3,2)*di(m+lz+3,3)+
     $                          i2(m+lx+4)*di(m+ly+4,2)*di(m+lz+4,3))
                           d2(42)=d2(42)+dens(prim,n)*
     $                         (i2(m+lx  )*i2(m+ly  )*d2i(m+lz  ,5)+
     $                          i2(m+lx+1)*i2(m+ly+1)*d2i(m+lz+1,5)+
     $                          i2(m+lx+2)*i2(m+ly+2)*d2i(m+lz+2,5)+
     $                          i2(m+lx+3)*i2(m+ly+3)*d2i(m+lz+3,5)+
     $                          i2(m+lx+4)*i2(m+ly+4)*d2i(m+lz+4,5))
                           d2(43)=d2(43)+dens(prim,n)*
     $                         (di(m+lx  ,3)*i2(m+ly  )*di(m+lz  ,3)+
     $                          di(m+lx+1,3)*i2(m+ly+1)*di(m+lz+1,3)+
     $                          di(m+lx+2,3)*i2(m+ly+2)*di(m+lz+2,3)+
     $                          di(m+lx+3,3)*i2(m+ly+3)*di(m+lz+3,3)+
     $                          di(m+lx+4,3)*i2(m+ly+4)*di(m+lz+4,3))
                           d2(44)=d2(44)+dens(prim,n)*
     $                         (i2(m+lx  )*di(m+ly  ,3)*di(m+lz  ,3)+
     $                          i2(m+lx+1)*di(m+ly+1,3)*di(m+lz+1,3)+
     $                          i2(m+lx+2)*di(m+ly+2,3)*di(m+lz+2,3)+
     $                          i2(m+lx+3)*di(m+ly+3,3)*di(m+lz+3,3)+
     $                          i2(m+lx+4)*di(m+ly+4,3)*di(m+lz+4,3))
                           d2(45)=d2(45)+dens(prim,n)*
     $                         (i2(m+lx  )*i2(m+ly  )*d2i(m+lz  ,6)+
     $                          i2(m+lx+1)*i2(m+ly+1)*d2i(m+lz+1,6)+
     $                          i2(m+lx+2)*i2(m+ly+2)*d2i(m+lz+2,6)+
     $                          i2(m+lx+3)*i2(m+ly+3)*d2i(m+lz+3,6)+
     $                          i2(m+lx+4)*i2(m+ly+4)*d2i(m+lz+4,6))
                           m=m+5
 902                    continue
                     end if
                     go to 15
   10             continue
                  call lnkerr('  fd2bnt: go to error ')
                     go to 15
   11             continue
                  call lnkerr('  fd2bnt: go to error ')
                     go to 15
 12               continue
                  call lnkerr('  fd2bnt: go to error ')
                     go to 15
 13               continue
                  call lnkerr('  fd2bnt: go to error ')
                     go to 15
   15             continue
c
  100          continue
  200       continue
  300    continue
  400 continue
c
c     ----- timing -----
c
c
c
      return
      end
