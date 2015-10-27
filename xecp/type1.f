*deck %W%  %G%
      subroutine type1(nt,np,alpp,dp,sp,exx,igbegn,igend,jbegn,jgend,
     $                 maxi2,maxp2,nx,ny,nz,istart,iend,jstart,jend,
     $                 ca,cb,lamax,lbmax,acoef1,len1,ptcf1,
     $                 qrad,angsum,thetak,phik,q,ang,xab,yab,zab,
     $                 ntop,ltop)
      implicit integer(a-z)
c
c     computes type 1 pseudopotential integrals.
c
c     ----- common -----
      real*8 pi,twopi,fpi,pi3haf,pi5hf2,piquart,sqpi,sqpi2
      common/pifac/pi,twopi,fpi,pi3haf,pi5hf2,piquart,sqpi,sqpi2
c
c     ----- arguments unchanged -----
      integer nt,np(*),nx(*),ny(*),nz(*)
      integer maxi2,maxp2
      integer istart,iend,jstart,jend
      integer igbegn,igend,jgbegn,jgend
      integer lamax,lbmax
      integer ptcf1(0:ltop,0:ltop,0:ltop)
      real*8 alpp(*),dp(*),exx(*)
      real*8 ca(3),cb(3)
      complex*16 acoef1(len1)
c     ----- scratch arrays -----
      real*8 qrad(maxp2,0:lamax+lbmax,0:lamax+lbmax)
      real*8 angsum(maxp2,0:lamax+lbmax,0:lamax+lbmax)
      real*8 thetak(maxp2),phik(maxp2)
      real*8 q(0:ntop,0:ltop)
      real*8 ang(0:ltop)
      real*8 xab(0:ltop+ltop), yab(0:ltop+ltop), zab(0:ltop+ltop)
c     ----- arguments returned -----
      real*8 sp(maxp2,maxi2)
c     ----- local variables -----
      real*8 ca2,cb2
      real*8 alpa,alpb,eabc,argab
      real*8 alpha,xk,yk,zk,rk
      real*8 pab1,pab2,pab3
      real*8 zero,one,two
      parameter (zero=0.d+00,one=1.0d+00,two=2.0d+00)
c
c
c
c     ----- find distance between sites a,b and ecp center c -----
      ca2=zero
      cb2=zero
      do 10 i=1,3
         ca2=ca2+ca(i)*ca(i)
         cb2=cb2+cb(i)*cb(i)
   10 continue
c
c     ----- compute basic radial integrals qrad over primitives
c           and sum over terms in effecive potential -----
      ltot=lamax+lbmax
      prim=0
      do 100 jgauss=jgbegn,jgend
         alpb=exx(jgauss)
         do 90 igauss=igbegn,igend
            alpa=exx(igauss)
c
            prim=prim+1
            eabc=-alpa*ca2-alpb*cb2
            xk=-two*(alpa*ca(1)+alpb*cb(1))
            yk=-two*(alpa*ca(2)+alpb*cb(2))
            zk=-two*(alpa*ca(3)+alpb*cb(3))
            rk=sqrt(xk*xk+yk*yk+zk*zk)
            if (rk.eq.zero) then
               xk=zero
               yk=zero
               zk=one
               thetak(prim)=zero
               phik(prim)=zero
            else
               xk=xk/rk
               yk=yk/rk
               zk=zk/rk
               thetak(prim)=acos(zk)
               phik(prim)=atan(yk/xk)
            endif
c
c           ----- loop over terms in the projector -----
            do 20 abc=0,ltot
               do 20 lambda=0,abc
                  qrad(prim,abc,lambda)=zero
   20       continue
            do 50 i=1,nt
               alpha=alpa+alpb+alpp(i)
               do 30 abc=0,ltot+np(i)
                  do 30 lambda=0,abc
                     q(abc,lambda)=zero
   30          continue
               call recur1 (np(i),ntop,ltot,q,alpha,rk,argab)
c
c              ----- exponential factor below includes a bias (-argab), 
c              where argab is negative. this is because exp(argab) was
c              included in the  qrad's to prevent overflow.
               do 40 abc=0,ltot
                  lbeg=mod(abc,2)
                  do 40 lambda=lbeg,abc,2
                     qrad(prim,abc,lambda)=qrad(prim,abc,lambda)
     $                        +q(np(i)+abc,lambda)*dp(i)*exp(eabc-argab)
   40          continue
   50       continue
   90    continue
  100 continue
c
c     ----- loop over shell block components.
      intc=0
      do 200 ii=istart,iend
         na=nx(ii)
         la=ny(ii)
         ma=nz(ii)
c
         do 190 iii=jstart,jend
            nb=nx(iii)
            lb=ny(iii)
            mb=nz(iii)
            intc=intc+1
            call facab (na,nb,ca(1),cb(1),xab)
            call facab (la,lb,ca(2),cb(2),yab)
            call facab (ma,mb,ca(3),cb(3),zab)
c
c           ----- compute angular integrals and combine.
c                 zero out the angular region.
            l=na+nb+la+lb+ma+mb
            do 110 abc=0,l
               do 110 lambda=0,abc
                  do 110 prim=1,maxp2
                     angsum(prim,abc,lambda)=zero
  110       continue
            do 150 i=0,na+nb
               pab1=xab(i)
               do 150 j=0,la+lb
                  pab2=pab1*yab(j)
                  do 150 k=0,ma+mb
                     pab3=pab2*zab(k)
                     ijk=i+j+k
                     lbeg=mod(ijk,2)
c
c                    ----- loop over all primitives -----
                     ind=ptcf1(i,j,k)
                     do 140 prim=1,maxp2
                           call omega1(i,j,k,thetak(prim),phik(prim),
     $                                 ang(0),acoef1(ind+1))
                           do 130 lambda=lbeg,ijk,2
                              angsum(prim,ijk,lambda)
     $                        =angsum(prim,ijk,lambda)+ang(lambda)*pab3
  130                      continue
  140                continue
  150       continue
c
c           ----- combine angular and radial integrals -----
            do 160 abc=0,ltot
               lbeg=mod(abc,2)
               do 160 lambda=lbeg,abc,2
                  do 155 prim=1,maxp2
                     sp(prim,intc)=sp(prim,intc)
     $                  +angsum(prim,abc,lambda)*qrad(prim,abc,lambda)
c                 call vmul(ss,angsum(1,abc,lambda),qrad(1,abc,lambda),maxp2)
c                 call vadd(sp(1,intc),sp(1,intc),ss,maxp2)
  155             continue
  160       continue
  190    continue
  200 continue
c     end shell loop.
c
      return
      end
