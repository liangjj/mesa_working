*deck  %W% %G%
      subroutine ecp(iatom,jatom,lamax,lbmax,nprimi,nprimj,
     $                   nprim,ncont,nat,nbtype,ntypes,maxp2,maxi2,
     $                   i1,j1,c,ex,noprim,ptprim,ptcont,cont,
     $                   maxmom,mini,maxi,minj,maxj,nx,ny,nz,sp,
     $                   ntpse,nlp,zlp,clp,acoef1,len1,ptcf1,
     $                   acoef2,len2,ptcf2,ntop,ltop,lpmax,z,a)
      implicit real*8(a-h,o-z)
c     ----- arguments unchanged -----
      integer iatom,jatom,lamax,lbmax,nprimi,nprimj,nprim
      integer ncont,nat,nbtype,ntypes,maxp2,maxi2
      integer i1,j1 
      real*8 c(3,nat),ex(nprim)
      integer noprim(nat,ntypes)
      integer ptprim(nat,ntypes),ptcont(nat,ntypes)
      real*8 cont(ncont)
      integer maxmom(ntypes)
      integer mini,maxi,minj,maxj
      integer nx(*),ny(*),nz(*)
      integer nlp(*),ntpse(*)
      real*8 zlp(*),clp(*)
      complex*16 acoef1(len1)
      integer ptcf1(0:ltop,0:ltop,0:ltop)
      complex*16 acoef2(*)
      integer ptcf2(0:ltop,0:ltop,0:ltop,0:lpmax)
      integer ntop,ltop
c     ----- arguments returned -----
      real*8 sp(maxp2,maxi2)
c     ----- scratch space -----
      real*8 z(*)
      integer a(*)
c     ----- local variables -----
      logical hasecp
      real*8 ca(3),cb(3)
      real*8 half
      parameter (half=5.0d-01)
c
c     controls computation of pseudopotential integrals.
c     details of the method are presented in
c         l.e. mcmurchie and e.r. davidson, j. comp. phys.,44,289(1981).
c         r.l. martin, to be published.
c
      common/pifac/pi,twopi,fpi,pi3haf,pi5hf2,piquart,sqpi,sqpi2
c
c     the parameter ltop sets the limit on the total angular momentum
c     which may be handled; i.e., ltop=6 implies f with f. 
c     the parameter ntop deals with the value of n in the ecp projector.
c     currently this is n=0,1,2, so ntop=lmax+2
c     parameter (ntop=8,ltop=6)
c     real*8 q(0:ntop,0:ltop)
c     real*8 ang(0:ltop)
c     real*8 xab(0:ltop+ltop),yab(0:ltop+ltop),zab(ltop+ltop)
c
c     ----- initialize the integral output array -----
      call rzero(sp,maxp2*maxi2)
      do 100 katom=1,nat
c
         do 90 i=1,3
            ca(i)=c(i,katom)-c(i,iatom)
            cb(i)=c(i,katom)-c(i,jatom)
   90    continue
         nt=0
         kk=0
         lmax=0
         hasecp=.false.
         do 10 ktype=nbtype+1,ntypes
            if(noprim(katom,ktype).ne.0) then
               hasecp=.true.
               lmax=max(lmax,maxmom(ktype))
               nt=nt+1
               ntpse(nt)=noprim(katom,ktype)
               do 5 jj=1,noprim(katom,ktype)
                  kk=kk+1
                  nlp(kk)=int(cont(ptcont(katom,ktype)
     $                             +noprim(katom,ktype)+jj-1)+half)
                  clp(kk)=cont(ptcont(katom,ktype)+jj-1)
                  zlp(kk)=ex(ptprim(katom,ktype)+jj-1)
    5          continue
            endif
   10    continue
         if(hasecp) then
c
c           ----- allocate core for type1 integrals -----
            qrad=1
            angsum=qrad+maxp2*(lamax+lbmax+1)*(lamax+lbmax+1)
            thetak=angsum+maxp2*(lamax+lbmax+1)*(lamax+lbmax+1)
            phik=thetak+maxp2
            q=phik+maxp2
            ang=q+(ntop+1)*(ltop+1)
            xab=ang+ltop+ltop+1
            yab=xab+ltop+ltop+1
            zab=yab+ltop+ltop+1
            call type1 (ntpse,nlp,zlp,clp,
     $                  sp,ex,i1,i1+nprimi-1,j1,j1+nprimj-1,
     $                  maxi2,maxp2,nx,ny,nz,mini,maxi,minj,maxj,ca,cb,
     $                  lamax,lbmax,acoef1,len1,ptcf1,z(qrad),z(angsum),
     $                  z(thetak),z(phik),z(q),z(ang),z(xab),
     $                  z(yab),z(zab),ntop,ltop)
c
c           ----- allocate core for type2 integrals -----
            qq=1
            atheta=qq+maxp2*(ltop+1)*(ltop+1)*(ltop+1)
            aphi=atheta+maxp2
            btheta=aphi+maxp2
            bphi=btheta+maxp2
            xab=bphi+maxp2
            yab=xab+ltop+ltop+1
            zab=yab+ltop+ltop+1
            q2=zab+ltop+ltop+1
            anga=q2+maxp2
            angb=anga+2*(lmax+1)*(lmax+lamax+1)*nprima
            top=angb+2*(lmax+1)*(lmax+lbmax+1)*nprimb
            call type2 (lmax,ntpse(2),nlp(ntpse(1)+1),zlp(ntpse(1)+1),
     $                  clp(ntpse(1)+1),sp,ex,maxi2,nprimi,nprimj,
     $                  nx,ny,nz,ca,cb,lamax,lbmax,i1,i1+nprimi-1,
     $                  j1,j1+nprimj-1,mini,maxi,minj,maxj,z(qq),
     $                  z(atheta),z(aphi),z(btheta),z(bphi),z(xab),
     $                  z(yab),z(zab),z(q2),z(anga),z(angb),acoef2,len2,
     $                  ptcf2,ltop)
         endif
  100 continue
c
      call smul(sp,sp,fpi,maxi2*maxp2)
c
c
      return
      end
