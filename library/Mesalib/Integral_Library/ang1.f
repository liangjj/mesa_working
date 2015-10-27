*deck @(#)ang1.f	5.1  11/6/94
      subroutine ang1 (n,l,m,xk,yk,zk,ang)
c***begin prologue     ang1.f
c***date written       yymmdd  
c***revision date      930224
c
c   february 24,1993   akr at csu
c      new common blocks to reflect enlarged ztab.
c   november 17,1992   rlm at lanl
c      cosmetology
c***keywords           
c***author             
c***source             @(#)ang1.f	5.1   11/6/94
c***purpose            
c
c     computes angular integral for type 1 integrals.
c***description
c
c   input:   n         the principal quantum number
c            l         the 
c            m
c            xk        
c            yk
c            zk
c   output:  ang       the angular integrals
c***references
c      l.e.mcmurchie and e.r.davidson, j.comp.phys.
c***routines called
c                      (none)
c
c***end prologue       ang1.f
      implicit  integer(a-z)
c
      real*8 xk,yk,zk,ang(9)
      real*8 angt,pre,aint,xkp,ykp,zkp
      real*8 zlm
      real*8 dfac
      real*8 pi,twopi,fpi,pi3haf,pi5hf2,piquart,sqpi,sqpi2
      real*8 one,zero
c
c     common/ztabcm/lf(7),lmf(49),lml(49),lmx(130),lmy(130),lmz(130)
      common/ztabcm/lf(13),lmf(169),lml(169),lmx(1036),lmy(1036),
     $ lmz(1036)
c     common/ztabwp/zlm(130)
      common/ztabwp/zlm(1036)
      common/dfac/dfac(30)
      common/pifac/pi,twopi,fpi,pi3haf,pi5hf2,piquart,sqpi,sqpi2
c
      parameter (zero=0.0d+00, one=1.0d+00)
c
      nlm=n+l+m+1
      do 90 lambda=1,nlm,2
         lmb=nlm-lambda
         l2=lmb+lmb+1
         angt=zero
         loc=lf(lmb+1)
         do 80 mu1=1,l2
            mu=mu1-1
            istart=lmf(loc+mu)
c           if ((n.and.1).ne.(lmx(istart).and.1)) go to 80
c           if ((l.and.1).ne.(lmy(istart).and.1)) go to 80
c           if ((m.and.1).ne.(lmz(istart).and.1)) go to 80
            if (and(n,1).ne.(and(lmx(istart),1))) go to 80
            if (and(l,1).ne.(and(lmy(istart),1))) go to 80
            if (and(m,1).ne.(and(lmz(istart),1))) go to 80
            pre=zero
            aint=zero
            iend=lml(loc+mu)
            do 70 i=istart,iend
               indx=lmx(i)
               indy=lmy(i)
               indz=lmz(i)
               if (indx.eq.0)  then
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
               aint=aint+zlm(i)*dfac(n+indx+1)*dfac(l+indy+1)
     $                  *dfac(m+indz+1)/dfac(nlm-1+indx+indy+indz+3)
   70       continue
            angt=angt+pre*aint*fpi
   80    continue
         ang(lmb+1)=angt
   90 continue
c
c
      return
      end
