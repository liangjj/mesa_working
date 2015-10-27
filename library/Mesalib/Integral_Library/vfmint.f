*deck @(#)vfmint.f	5.1   11/6/94
      subroutine vfmint(i2,half,temp,lenblk,nv,angmom,imax,jmax,
     #                  mmax,lmax,npint,nroots,px,py,pz,lenxyz,mintyp,
     #                  nocart,nbtype,minkl,nijkl,test,nbatch,
     #                  a,b,ia,fa,ib,fb,nkl,numkl,temp1,nij,lent)
c
c vectorised module to form primitive integerals
c
c paul saxe                 13 july 1984           lanl
c
      implicit integer (a-h,o-z)
c
c      real i2(nv,3,0:imax,0:jmax,0:mmax,0:lmax)
      real*8 i2(*)
c     real half(fa,fb,nkl,lenblk),temp(nijkl)
      real*8 half(*),temp(lent)
      real*8 a(ia,fa),b(ib,fb),temp1(lent)
      integer test(nijkl)
      integer angmom(4)
      integer px(lenxyz),py(lenxyz),pz(lenxyz)
      integer mintyp(nbtype),nocart(nbtype)
      parameter (nindex=56)
      integer ipx(nindex),ipy(nindex),ipz(nindex)
      integer jpx(nindex),jpy(nindex),jpz(nindex)
      integer kpx(nindex),kpy(nindex),kpz(nindex)
      integer lpx(nindex),lpy(nindex),lpz(nindex)
c     ----- work out how many integrals to accumulate -----
c
      numn=min(nocart(angmom(3))*nocart(angmom(4)),lent/nijkl)
      numn=min(numn,64)
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
      junk=3*nv
      do 1 i=mini,maxi
         ipx(i)=     junk*px(i)
         ipy(i)=  nv+junk*py(i)
         ipz(i)=2*nv+junk*pz(i)
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
            minn=n+1
            nkflf=0
            pt=0
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
                  nkflf=nkflf+1
c
                  go to (5,6,7,8,9,10,11,12,13),nroots
c
                  call lnkerr('too many roots needed')
c
    5             continue
                     do 50 prim=1,nbatch
                        temp(pt+prim)=i2(prim+lx)*i2(prim+ly)*
     #                                 i2(prim+lz)
   50                continue
                     go to 15
    6             continue
                     m=1
                     do 60 prim=pt+1,pt+nbatch
                        temp(prim)=i2(m+lx)*i2(m+ly)*i2(m+lz)+
     #                                i2(m+1+lx)*i2(m+1+ly)*i2(m+1+lz)
                        m=m+2
   60                continue
                     go to 15
    7             continue
                     m=1
                     do 70 prim=pt+1,pt+nbatch
                        temp(prim)=i2(m+lx)*i2(m+ly)*i2(m+lz)+
     #                                i2(m+1+lx)*i2(m+1+ly)*i2(m+1+lz)+
     #                                i2(m+2+lx)*i2(m+2+ly)*i2(m+2+lz)
                        m=m+3
   70                continue
                     go to 15
    8             continue
                     m=1
                     do 80 prim=pt+1,pt+nbatch
                        temp(prim)=i2(m+lx)*i2(m+ly)*i2(m+lz)+
     #                                i2(m+1+lx)*i2(m+1+ly)*i2(m+1+lz)+
     #                                i2(m+2+lx)*i2(m+2+ly)*i2(m+2+lz)+
     #                                i2(m+3+lx)*i2(m+3+ly)*i2(m+3+lz)
                        m=m+4
   80                continue
                     go to 15
    9             continue
                     m=1
                     do 90 prim=pt+1,pt+nbatch
                        temp(prim)=i2(m+lx)*i2(m+ly)*i2(m+lz)+
     #                                i2(m+1+lx)*i2(m+1+ly)*i2(m+1+lz)+
     #                                i2(m+2+lx)*i2(m+2+ly)*i2(m+2+lz)+
     #                                i2(m+3+lx)*i2(m+3+ly)*i2(m+3+lz)+
     #                                i2(m+4+lx)*i2(m+4+ly)*i2(m+4+lz)
                        m=m+5
   90                continue
                     go to 15
   10             continue
                     m=1
                     do 110 prim=pt+1,pt+nbatch
                        temp(prim)=i2(m+lx)*i2(m+ly)*i2(m+lz)+
     #                                i2(m+1+lx)*i2(m+1+ly)*i2(m+1+lz)+
     #                                i2(m+2+lx)*i2(m+2+ly)*i2(m+2+lz)+
     #                                i2(m+3+lx)*i2(m+3+ly)*i2(m+3+lz)+
     #                                i2(m+4+lx)*i2(m+4+ly)*i2(m+4+lz)+
     #                                i2(m+5+lx)*i2(m+5+ly)*i2(m+5+lz)
                        m=m+6
  110                continue
                     go to 15
   11             continue
                     m=1
                     do 120 prim=pt+1,pt+nbatch
                        temp(prim)=i2(m+lx)*i2(m+ly)*i2(m+lz)+
     #                                i2(m+1+lx)*i2(m+1+ly)*i2(m+1+lz)+
     #                                i2(m+2+lx)*i2(m+2+ly)*i2(m+2+lz)+
     #                                i2(m+3+lx)*i2(m+3+ly)*i2(m+3+lz)+
     #                                i2(m+4+lx)*i2(m+4+ly)*i2(m+4+lz)+
     #                                i2(m+5+lx)*i2(m+5+ly)*i2(m+5+lz)+
     #                                i2(m+6+lx)*i2(m+6+ly)*i2(m+6+lz)
                        m=m+7
  120                continue
                     go to 15
 12                  continue
                     m=1
                     do 130 prim=pt+1,pt+nbatch
                        temp(prim)=i2(m+lx)*i2(m+ly)*i2(m+lz)+
     #                                i2(m+1+lx)*i2(m+1+ly)*i2(m+1+lz)+
     #                                i2(m+2+lx)*i2(m+2+ly)*i2(m+2+lz)+
     #                                i2(m+3+lx)*i2(m+3+ly)*i2(m+3+lz)+
     #                                i2(m+4+lx)*i2(m+4+ly)*i2(m+4+lz)+
     #                                i2(m+5+lx)*i2(m+5+ly)*i2(m+5+lz)+
     #                                i2(m+6+lx)*i2(m+6+ly)*i2(m+6+lz)+
     #                                i2(m+7+lx)*i2(m+7+ly)*i2(m+7+lz)
                        m=m+8
  130                continue
                     go to 15
 13                  continue
                     m=1
                     do 140 prim=pt+1,pt+nbatch
                        temp(prim)=i2(m+lx)*i2(m+ly)*i2(m+lz)+
     #                                i2(m+1+lx)*i2(m+1+ly)*i2(m+1+lz)+
     #                                i2(m+2+lx)*i2(m+2+ly)*i2(m+2+lz)+
     #                                i2(m+3+lx)*i2(m+3+ly)*i2(m+3+lz)+
     #                                i2(m+4+lx)*i2(m+4+ly)*i2(m+4+lz)+
     #                                i2(m+5+lx)*i2(m+5+ly)*i2(m+5+lz)+
     #                                i2(m+6+lx)*i2(m+6+ly)*i2(m+6+lz)+
     #                                i2(m+7+lx)*i2(m+7+ly)*i2(m+7+lz)+
     #                                i2(m+8+lx)*i2(m+8+ly)*i2(m+8+lz)
                        m=m+9
  140                continue
   15             continue
c
                  if (nkflf.eq.numn.or.(k.eq.maxk.and.l.eq.maxl)) then
c
                     call rzero(temp1,nijkl*nkflf)
c
                     maxb=0
                     if (nkflf.gt.64) call lnkerr('ask paul about '//
     #                  'nkflf')
                     do 21 junk=1,(nkflf+63)/64
                        minb=maxb+1
                        maxb=min(nkflf,maxb+64)
                        nb=maxb-minb+1
                        junk1=0
                        ptx=0
                        pty=0
                        do 102 iqb=1,nb
                           do 101 ipt=1,nbatch
                              temp1(ptx+test(ipt))=temp(pty+ipt)
  101                      continue
                           ptx=ptx+nijkl
                           pty=pty+nbatch
  102                   continue
   21                continue
c
                     call rzero(temp,fa*ib*numkl*nkflf)
c
                     do 30 fapt=1,fa
                        do 29 iapt=1,ia
                           if (a(iapt,fapt).ne.0.0d+00) then
                              call saxpy(ib*numkl*nkflf,a(iapt,fapt),
     #                                   temp1(iapt),ia,
     #                                   temp(fapt),fa)
                           end if
   29                   continue
   30                continue
c
                     do 40 fbpt=1,fb
                        do 39 ibpt=1,ib
                           if (b(ibpt,fbpt).ne.0.0d+00) then
                              do 38 fapt=1,fa
                                 abkln=minn+lenblk*(fapt-1+fa*(fbpt-1+
     #                                            fb*(minkl-1)))
                                 abkl=fapt+fa*(ibpt-1)
c                                if (numkl.gt.nkflf) then
                                 if (numkl.lt.0) then
                                    do 36 kflf=1,nkflf
                                       call saxpy(numkl,b(ibpt,fbpt),
     #                                      temp(abkl),fa*ib,
     #                                      half(abkln),fa*fb*lenblk)
                                       abkl=abkl+fa*ib*numkl
                                       abkln=abkln+1
   36                               continue
                                 else
                                    do 37 kl=1,numkl
                                       call saxpy(nkflf,b(ibpt,fbpt),
     #                                        temp(abkl),fa*ib*numkl,
     #                                        half(abkln),1)
                                       abkl=abkl+fa*ib
                                       abkln=abkln+lenblk*fa*fb
   37                               continue
                                 end if
   38                         continue
                           end if
   39                   continue
   40                continue
c
                     pt=0
                     minn=n+1
                     nkflf=0
                  else
                     pt=pt+nbatch
                  end if
c
  100          continue
  200       continue
  300    continue
  400 continue
c
c
      return
      end
