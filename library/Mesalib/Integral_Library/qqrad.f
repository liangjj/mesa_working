*deck  %W% %G%
      subroutine qqrad(l,lmalo,lmahi,lmblo,lmbhi,nhi,np,alpp,dp,ntcnt,
     $                 exx,qq,maxi2,maxp2)
      implicit real*8(a-h,o-z)
c
c     evaluate the type 2 radial integrals.
c
      common/dist/cax,cay,caz,ca,ca2,cbx,cby,cbz,cb,cb2
      common/ptwtdat/ptpow(50,11),f(50,9,9),pt(50)
      common/ptwt/npts
      common/const/zero,one,two,three,four,five,six,ten
      common/prims/igbegn,igend,jgbegn,jgend
      common/pifac/pi,twopi,fpi,pi3haf,pi5hf2,piquart,sqpi,sqpi2
      common/io/in,iout
      dimension np(1), alpp(1), dp(1), exx(1)
      dimension qq(11,9,9,1)
      dimension qtemp(11,9,9)
      parameter (hundrd=1.0d+02,eps1=1.d-15)
c
      ilo=1
      ihi=ntcnt
      iprim=0
c
c     loop over primitive gaussians.
c
      do 150 jgauss=jgbegn,jgend
         alpb=exx(jgauss)
         rkb=two*alpb*cb
c
         do 150 igauss=igbegn,igend
            alpa=exx(igauss)
            rka=two*alpa*ca
c
            iprim=iprim+1
            call rzero (qq(1,1,1,iprim),11*9*9)
c
c           loop over terms in the projector.
            do 140 i=ilo,ihi
               npi=np(i)
               alpha=alpa+alpb+alpp(i)
               argsum=((rka+rkb)*(rka+rkb))/(two*alpha)
c              ----- determine which method is to be used
c                    to compute radial integrals.
               if (argsum.ge.hundrd) then
c
c                 ----- pts and wts method on (-inf,inf) to be used -----
c                 result includes a factor exp(-(xka+xkb)**2/(four*alpha)) 
c                 to prevent overflow.
                  arg=(-alpa*alpb*(ca-cb)*(ca-cb)
     $                  -alpp(i)*(alpa*ca*ca+alpb*cb*cb))/alpha
                  exparg=exp(arg)*dp(i)/sqrt(alpha)
c                 determine limit for radial integrals and skip evaluation 
c                 if less than threshold eps1.
                  rc=(rka+rkb)/(two*alpha)
                  if (rc.le.one) nlim=npi
                  if (rc.gt.one) nlim=nhi+npi-1
                  if (rka.eq.zero) then
                     bessa=one
                  else 
                     bessa=one/(two*rc*rka)
                  endif
                  if(rkb.eq.zero) then
                     bessb=one
                  else
                     bessb=one/(two*rc*rkb)
                  endif
                  qlim=rc**nlim*sqpi*bessa*bessb
                  if ((qlim*fpi*abs(exparg)).lt.eps1) go to 140
               else
c
c                 ----- pts and wts method on (0,inf) to be used -----
                  exparg=exp(-alpa*ca*ca-alpb*cb*cb)*dp(i)/sqrt(alpha)
               endif
c
               lmahm1=max(lmahi-1,1)
               lmbhm1=max(lmbhi-1,1)
               call rzero (qtemp(1,1,1),11*9*9)
c
c              ----- branch for special cases -----
               if (rka.eq.zero.and.rkb.eq.zero) then
c                 ----- one center case -----
                  if ((lmalo.eq.1).and.(lmblo.eq.1)) then
                     call ptprep (npi,nhi,1,1,1,1,alpha,rka,rkb,argsum)
                     call quadr (1,1,1,1,1,nhi,qtemp(1,1,1))
                  endif
               else if(rka.eq.zero) then
c                 ----- two center case -----
                  if (lmalo.eq.1) then
                     call ptprep (npi,nhi,1,1,lmbhm1,lmbhi,alpha,
     $                            rka,rkb,argsum)
                     call quadr (1,1,lmbhm1,lmbhi,1,nhi,qtemp(1,1,1))
                     if (lmbhi.gt.2) then
                        fkb=one/rkb
                        do 60 lb=lmbhi,lmblo+2,-1
                           lbtru=lb-1
                           do 40 j=1,npts
                              f(j,1,lb-2)=f(j,1,lb)
     $                              +(2*lbtru-1)*(fkb/pt(j))*f(j,1,lb-1)
                              qtemp(1,1,lb-2)=qtemp(1,1,lb-2)
     $                                       +ptpow(j,1)*f(j,1,lb-2)
   40                      continue
                           do 50 n=1,nhi-1
                              qtemp(n+1,1,lb-2)=qtemp(n+1,1,lb)
     $                              +(2*lbtru-1)*fkb*qtemp(n,1,lb-1)
   50                      continue
   60                   continue
                     endif
                  endif
               else if(rkb.eq.zero) then
                  call ptprep (npi,nhi,lmahm1,lmahi,1,1,alpha,rka,rkb,
     $                         argsum)
                  if (lmblo.eq.1) then
                     call quadr (lmahm1,lmahi,1,1,1,nhi,qtemp(1,1,1))
                     if (lmahi.gt.2) then
                        fka=one/rka
                        do 100 la=lmahi,lmalo+2,-1
                           latru=la-1
                           do 80 j=1,npts
                              f(j,la-2,1)=f(j,la,1)
     $                            +(2*latru-1)*(fka/pt(j))*f(j,la-1,1)
                              qtemp(1,la-2,1)=qtemp(1,la-2,1)
     $                                       +ptpow(j,1)*f(j,la-2,1)
   80                      continue
                           do 90 n=1,nhi-1
                              qtemp(n+1,la-2,1)=qtemp(n+1,la,1)
     $                            +(2*latru-1)*fka*qtemp(n,la-1,1)
   90                      continue
  100                   continue
                     endif
                  endif
               else
c                 ----- general three center case -----
                  call ptprep (npi,nhi,lmahm1,lmahi,lmbhm1,lmbhi,alpha,
     $                         rka,rkb,argsum)
                  call recurf (lmalo,lmahi,lmblo,lmbhi,rka,rkb)
                  call quadr (lmalo,lmahi,lmblo,lmbhi,1,1,qtemp(1,1,1))
                  if (nhi.gt.1) then
                     call quadr (lmahm1,lmahi,lmbhm1,lmbhi,2,nhi,
     $                           qtemp(1,1,1))
                     if (lmahi.gt.2.or.lmbhi.gt.2) then
                        call recur2 (nhi-1,lmalo,lmahi,lmblo,lmbhi,
     $                               rka,rkb,qtemp(1,1,1))
                     endif
                  endif
               endif
  120          continue
c
               do 130 lama=lmalo,lmahi
               do 130 lamb=lmblo,lmbhi
               do 130 n=1,nhi
                  qq(n,lama,lamb,iprim)=qq(n,lama,lamb,iprim)
     $                                 +qtemp(n,lama,lamb)*exparg
  130          continue
c           end projector loop.
  140       continue
c     end primitives loop.
  150 continue
c
c
      return
      end
