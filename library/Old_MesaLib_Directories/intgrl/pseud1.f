*deck @(#)pseud1.f	5.1 11/6/94
      subroutine pseud1 (nt,np,alpp,dp,sp,exx,maxi2,maxp2,nx,ny,nz)
      implicit real*8(a-h,o-z)
c
c     computes type 1 pseudopotential integrals.
c
      common/center/xa,ya,za,xb,yb,zb,xc,yc,zc
      common/dist/cax,cay,caz,ca,ca2,cbx,cby,cbz,cb,cb2
      common/qstore/q(13,11),alpha,rk,t
      common/argab/argab,expab
      common/angmax/lamax,lbmax
      common/limit/istart,jstart,iend,jend
      common/const/zero,one,two,three,four,five,six,ten
      common/prims/igbegn,igend,jgbegn,jgend
      common/pifac/pi,twopi,fpi,pi3haf,pi5hf2,piquart,sqpi,sqpi2
      dimension np(1), alpp(1), dp(1), exx(1),nx(1),ny(1),nz(1)
      dimension sp(maxp2,maxi2)
      dimension ang(11), xab(11), yab(11), zab(11)
      dimension qsum(11,11), angsum(11,11)
c
      ltot=lamax+lbmax-2
      ltot1=ltot+1
c
c     compute basic radial integrals q and sum over terms in potential.
c
c     loop over primitive gaussians
      iprim=0
      do 100 jgauss=jgbegn,jgend
         alpb=exx(jgauss)
         do 100 igauss=igbegn,igend
            alpa=exx(igauss)
c
            iprim=iprim+1
c
c           zero out radial and angular integral regions.
            do 10 j=1,ltot1
            do 10 k=1,j
               qsum(j,k)=zero
   10       continue
            eabc=-alpa*ca2-alpb*cb2
            xk=-two*(alpa*cax+alpb*cbx)
            yk=-two*(alpa*cay+alpb*cby)
            zk=-two*(alpa*caz+alpb*cbz)
            rk=sqrt(xk*xk+yk*yk+zk*zk)
            if (rk.eq.zero) then
               xk=zero
               yk=zero
               zk=one
            else
               xk=xk/rk
               yk=yk/rk
               zk=zk/rk
            endif
c
c           loop over terms in the projector.
            do 50 i=1,nt
               alpha=alpa+alpb+alpp(i)
               kstart=np(i)
               call recur1 (np(i),ltot)
c
c              exponential factor below includes a bias (-argab), where argab
c              is negative. this is because exp(argab) was included in the q's
c              to prevent overflow.
               do 40 j=1,ltot1
               do 40 k=1,j,2
                  kb=j-k+1
               qsum(j,kb)=qsum(j,kb)
     $                   +q(kstart+j,kb)*dp(i)*exp(eabc-argab)
   40          continue
   50       continue
c
c           loop over shell block components.
            intc=0
            do 90 ii=istart,iend
               na=nx(ii)
               la=ny(ii)
               ma=nz(ii)
c
               do 90 iii=jstart,jend
                  nb=nx(iii)
                  lb=ny(iii)
                  mb=nz(iii)
                  intc=intc+1
                  call facab (na,nb,cax,cbx,xab)
                  call facab (la,lb,cay,cby,yab)
                  call facab (ma,mb,caz,cbz,zab)
c
c                 compute angular integrals and combine.
                  nanb=na+nb+1
                  lalb=la+lb+1
                  mamb=ma+mb+1
c                 zero out the angular region.
                  ltot1=nanb+lalb+mamb-2
                  do 60 j=1,ltot1
                  do 60 k=1,j
   60                angsum(j,k)=zero
                  do 70 i=1,nanb
                     pab1=xab(i)
                     do 70 j=1,lalb
                        pab2=pab1*yab(j)
                        do 70 k=1,mamb
                           pab3=pab2*zab(k)
                           ijk=i+j+k-2
                           call ang1 (i-1,j-1,k-1,xk,yk,zk,ang)
                           do 70 lambda=1,ijk,2
                              lmb=ijk-lambda+1
                              angsum(ijk,lmb)=angsum(ijk,lmb)
     $                                       +ang(lmb)*pab3
   70             continue
c
c                 now combine angular and radial integrals.
                  ss=zero
                  do 80 j=1,ltot1
                  do 80 k=1,j,2
                     kb=j-k+1
                     ss=ss+angsum(j,kb)*qsum(j,kb)
   80             continue
                  sp(iprim,intc)=sp(iprim,intc)+ss
   90       continue
c           end shell loop.
  100 continue
c     end primitives loop.
c
c
c
      return
      end
