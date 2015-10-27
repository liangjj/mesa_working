*deck @(#)fmint.f	1.1  11/30/90
      subroutine fmint(x,y,z,lenblk,rv,angmom,imax,jmax,
     #                  mmax,lmax,nroots,px,py,pz,lenxyz,mintyp,
     #                  nocart,nbtype,nv,pint,lnpint,ci,npi,nci,
     #                  cj,npj,ncj,ck,npk,nck,cl,npl,ncl,index,len,
     #                  cint)
c
c***begin prologue     vfmint
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
c***source             @(#)fmint.f	1.1   11/30/90
c***purpose
c***description
c***references         (none)
c
c***routines called    (none)
c
c***end prologue       vfmint
c
      implicit integer (a-h,o-z)
c
c     real*8 i2(rv,0:imax,0:jmax,0:mmax,0:lmax,3)
      real*8 x(*),y(*),z(*)
      real*8 pint(nv,*),ci(npi,nci),cj(npj,ncj),ck(npk,nck),cl(npl,ncl)
      real*8 cint(lenblk,nci,ncj,nck,ncl)
      real*8 sdot,cij,cijk,cijkl
      integer index(len,6)
      integer angmom(4)
      integer px(lenxyz),py(lenxyz),pz(lenxyz)
      integer mintyp(nbtype),nocart(nbtype)
      integer ipx(35),ipy(35),ipz(35)
      integer jpx(35),jpy(35),jpz(35)
      integer kpx(35),kpy(35),kpz(35)
      integer lpx(35),lpy(35),lpz(35)
c
c     ----- timing -----
c
c
c     ----- set up minimum and maximum functions in each shell -----
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
      junk=rv*(imax+1)*(jmax+1)*(mmax+1)*(lmax+1)
      do 1 i=mini,maxi
         ipx(i)=rv*px(i)
         ipy(i)=rv*py(i)
         ipz(i)=rv*pz(i)
    1 continue
      junk=rv*(imax+1)
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
c     ----- determine how many integrals we can hold -----
c
      maxgrp=lnpint/nv
      if (maxgrp.lt.1) call lnkerr('not enough space for primitives')
      ngrp=0
      stgrp=0
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
                  ngrp=ngrp+1
                  n=n+1
c
                  go to (5,6,7,8,9,10,11),nroots
c
                  call lnkerr('bad number of rys roots')
c
    5             continue
                     do 50 prim=1,nv
                        pint(prim,ngrp)=x(prim+lx)*y(prim+ly)*
     #                                z(prim+lz)
   50                continue
                     go to 15
    6             continue
                     m=1
                     do 60 prim=1,nv
                        pint(prim,ngrp)=(x(m+lx)*y(m+ly)*z(m+lz)+
     #                               x(m+1+lx)*y(m+1+ly)*z(m+1+lz))
                        m=m+2
   60                continue
                     go to 15
    7             continue
                     m=1
                     do 70 prim=1,nv
                        pint(prim,ngrp)=(x(m+lx)*y(m+ly)*z(m+lz)+
     #                               x(m+1+lx)*y(m+1+ly)*z(m+1+lz)+
     #                               x(m+2+lx)*y(m+2+ly)*z(m+2+lz))
                        m=m+3
   70                continue
                     go to 15
    8             continue
                     m=1
                     do 80 prim=1,nv
                        pint(prim,ngrp)=(x(m+lx)*y(m+ly)*z(m+lz)+
     #                               x(m+1+lx)*y(m+1+ly)*z(m+1+lz)+
     #                               x(m+2+lx)*y(m+2+ly)*z(m+2+lz)+
     #                               x(m+3+lx)*y(m+3+ly)*z(m+3+lz))
                        m=m+4
   80                continue
                     go to 15
    9             continue
                     m=1
                     do 90 prim=1,nv
                        pint(prim,ngrp)=(x(m+lx)*y(m+ly)*z(m+lz)+
     #                               x(m+1+lx)*y(m+1+ly)*z(m+1+lz)+
     #                               x(m+2+lx)*y(m+2+ly)*z(m+2+lz)+
     #                               x(m+3+lx)*y(m+3+ly)*z(m+3+lz)+
     #                               x(m+4+lx)*y(m+4+ly)*z(m+4+lz))
                        m=m+5
   90                continue
                     go to 15
   10             continue
                     m=1
                     do 110 prim=1,nv
                        pint(prim,ngrp)=(x(m+lx)*y(m+ly)*z(m+lz)+
     #                               x(m+1+lx)*y(m+1+ly)*z(m+1+lz)+
     #                               x(m+2+lx)*y(m+2+ly)*z(m+2+lz)+
     #                               x(m+3+lx)*y(m+3+ly)*z(m+3+lz)+
     #                               x(m+4+lx)*y(m+4+ly)*z(m+4+lz)+
     #                               x(m+5+lx)*y(m+5+ly)*z(m+5+lz))
                        m=m+6
  110                continue
                     go to 15
   11             continue
                     m=1
                     do 120 prim=1,nv
                        pint(prim,ngrp)=(x(m+lx)*y(m+ly)*z(m+lz)+
     #                               x(m+1+lx)*y(m+1+ly)*z(m+1+lz)+
     #                               x(m+2+lx)*y(m+2+ly)*z(m+2+lz)+
     #                               x(m+3+lx)*y(m+3+ly)*z(m+3+lz)+
     #                               x(m+4+lx)*y(m+4+ly)*z(m+4+lz)+
     #                               x(m+5+lx)*y(m+5+ly)*z(m+5+lz)+
     #                               x(m+6+lx)*y(m+6+ly)*z(m+6+lz))
                        m=m+7
  120                continue
   15             continue
c
c                 ----- transform the primitive integrals we've accumulated
c
                  if (ngrp.eq.maxgrp.or.n.eq.lenblk) then
c
                     do 71 prim=1,nv
                        ip=index(prim,1)
                        jp=index(prim,2)
                        kp=index(prim,3)
                        lp=index(prim,4)
                        do 69 ic=1,nci
                           if (ci(ip,ic).eq.0.0d+00) go to 69
                           do 68 jc=1,ncj
                              cij=ci(ip,ic)*cj(jp,jc)
                              if (cij.eq.0.0d+00) go to 68
                              do 67 kc=1,nck
                                 cijk=cij*ck(kp,kc)
                                 if (cijk.eq.0.0d+00) go to 67
                                 do 66 lc=1,ncl
                                    cijkl=cijk*cl(lp,lc)
                                    if (cijkl.eq.0.0d+00) go to 66
                                    do 65 grp=1,ngrp
                                       cint(stgrp+grp,ic,jc,kc,lc)=
     #                                 cint(stgrp+grp,ic,jc,kc,lc)+
     #                                   cijkl*pint(prim,grp)
   65                               continue
   66                            continue
   67                         continue
   68                      continue
   69                   continue
   71                continue
c
                     stgrp=n
                     ngrp=0
                  end if
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
