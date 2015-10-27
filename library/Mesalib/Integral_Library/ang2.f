*deck @(#)ang2.f	5.1  11/6/94
      subroutine ang2 (na,la,ma,l,xk,yk,zk,ang)
c***begin prologue     ang2.f
c***date written       yymmdd  
c***revision date      930224      
c   february 24,1993   akr at csu
c      new common blocks to reflect enlarged ztab
c
c***keywords           
c***author             
c***source             @(#)ang2.f	5.1   11/6/94
c***purpose            
c     computes angular integral for type 2 integrals......
c***description
c     
c    
c
c***references
c
c***routines called
c                      (none)
c
c***end prologue       ang2.f
      implicit integer(a-z)
c
      real*8 ang(18,7,9),xk,yk,zk
      real*8 angt,pre,aint,xkp,ykp,zkp
      real*8 zero,one
      real*8 zlm
      real*8 dfac
      real*8 pi,twopi,fpi,pi3haf,pi5hf2,piquart,sqpi,sqpi2
c
c     common/ztabcm/lf(7),lmf(49),lml(49),lmx(130),lmy(130),lmz(130)
      common/ztabcm/lf(13),lmf(169),lml(169),lmx(1036),lmy(1036),
     $ lmz(1036)
c     common/ztabwp/zlm(130)
      common/ztabwp/zlm(1036)
      common/dfac/dfac(30)
      common/pifac/pi,twopi,fpi,pi3haf,pi5hf2,piquart,sqpi,sqpi2
c
      parameter (zero=0.0d+00,one=1.0d+00)
c
      loc1=lf(l+1)
      mhi=l+l+1
      iabc=0
      na1=na+1
      la1=la+1
      ma1=ma+1
      do 130 ia=1,na1
      do 130 ib=1,la1
      do 130 ic=1,ma1
         iabc=iabc+1
         lambhi=ia+ib+ic-3+l+1
         lamblo=max(l-ia-ib-ic+3,0)+1
         if (mod(lambhi-lamblo,2).ne.0) lamblo=lamblo+1
         do 120 m=1,mhi
            mstart=lmf(loc1+m-1)
            mend=lml(loc1+m-1)
            do 110 lamb=1,lambhi
               ang(iabc,m,lamb)=zero
               if (lamb.lt.lamblo) go to 100
c              if (((lamb-lamblo).and.1).ne.0) go to 100
               if (and(lamb-lamblo,1).ne.0) go to 100
               l2=lamb+lamb-1
               angt=zero
               loc2=lf(lamb)
               do 90 mu=1,l2
                  istart=lmf(loc2+mu-1)
                  if (and(ia-1+lmx(mstart)+lmx(istart),1).ne.0) go to 90
                  if (and(ib-1+lmy(mstart)+lmy(istart),1).ne.0) go to 90
                  if (and(ic-1+lmz(mstart)+lmz(istart),1).ne.0) go to 90
c                 if (((ia-1+lmx(mstart)+lmx(istart)).and.1).ne.0)
c    $               go to 90
c                 if (((ib-1+lmy(mstart)+lmy(istart)).and.1).ne.0)
c    $               go to 90
c                 if (((ic-1+lmz(mstart)+lmz(istart)).and.1).ne.0)
c    $               go to 90
                  pre=zero
                  iend=lml(loc2+mu-1)
                  aint=zero
                  do 80 i=istart,iend
                     indx=lmx(i)
                     indy=lmy(i)
                     indz=lmz(i)
                     if (indx.eq.0) then
                        xkp=one
                     else
                        xkp=xk**indx
                     endif
                     if (indy.eq.0) then
                        ykp=one
                     else
                        ykp=yk**indy
                     endif
                     if (indz.eq.0) then
                        zkp=one
                     else
                        zkp=zk**indz
                     endif
                     pre=pre+zlm(i)*xkp*ykp*zkp
                     do 70 j=mstart,mend
                        mndx=lmx(j)
                        mndy=lmy(j)
                        mndz=lmz(j)
                        aint=aint+zlm(i)*zlm(j)*dfac(ia+indx+mndx)
     $                    *dfac(ib+indy+mndy)*dfac(ic+indz+mndz)
     $                    /dfac(ia+ib+ic+indx+mndx+indy+mndy+indz+mndz)
   70                continue
   80             continue
                  angt=angt+pre*aint*sqrt(fpi**3)
   90          continue
               ang(iabc,m,lamb)=angt
  100          continue
  110       continue
  120    continue
  130 continue
c
c
      return
      end
