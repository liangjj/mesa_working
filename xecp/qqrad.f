*deck %W%  %G%
      subroutine qqrad(l,lamax,lbmax,nmax,np,alpp,dp,ntcnt,
     $                 igbegn,igend,jgbegn,jgend,exx,
     $                 ca,cb,qq,qtemp,maxi2,maxp2)
      implicit real*8(a-h,o-z)
c
c     evaluate the type 2 radial integrals.
c
c     ----- arguments unchanged -----
      integer l,lamax,lbmax,np(1),ntcnt
      integer igbegn,igend,jgbegn,jgend
      integer maxp2,maxi2
      real*8 alpp(*),dp(*),exx(*)
      real*8 ca,cb
c     ----- arguments returned -----
      real*8 qq(0:nmax,0:lamax,0:lbmax,maxp2)
c     ----- arguments scratch -----
      real*8 qtemp(0:nmax,0:lamax,0:lbmax)
c     ----- local variables -----
      common/ptwtdat/ptpow(50,0:6),f(50,0:6,0:6),pt(50)
      common/ptwt/npts
      common/pifac/pi,twopi,fpi,pi3haf,pi5hf2,piquart,sqpi,sqpi2
      common/io/in,iout
      real*8 zero,one,two,four
      parameter (zero=0.0d+00,one=1.0d+00,two=2.0d+00,four=4.0d+00)
      parameter (hundrd=1.0d+02,eps1=1.d-15)
c
      ilo=1
      ihi=ntcnt
c
c     ----- generate the range of lambda and n -----
      n=lamax+lbmax
      lahi=l+lamax
      lalo=max(l-lamax,0)
      if(mod(l+lamax,2).eq.0) then
         if(mod(lalo,2).ne.0) lalo=lalo+1
      else
         if(mod(lalo,2).eq.0) lalo=lalo+1
      endif
      lbhi=l+lbmax
      lblo=max(l-lbmax,0)
      if(mod(l+lbmax,2).eq.0) then
         if(mod(lblo,2).ne.0) lblo=lblo+1
      else
         if(mod(lblo,2).eq.0) lblo=lblo+1
      endif
c
c     ----- loop over primitive gaussians -----
      prim=0
      do 150 jgauss=jgbegn,jgend
         alpb=exx(jgauss)
         rkb=two*alpb*cb
c
         do 150 igauss=igbegn,igend
            alpa=exx(igauss)
            rka=two*alpa*ca
c
            prim=prim+1
            call rzero (qq(0,0,0,prim),343)
c
c           ----- loop over terms in the projector -----
            do 140 i=ilo,ihi
               npi=np(i)
               alpha=alpa+alpb+alpp(i)
               argsum=((rka+rkb)*(rka+rkb))/(two*alpha)
c              ----- determine which method is to be used
c                    to compute radial integrals.
               if (argsum.ge.hundrd) then
c                 ----- pts and wts method on (-inf,inf) to be used -----
c                 result includes a factor exp(-(xka+xkb)**2/(four*alpha)) 
c                 to prevent overflow.
                  arg=(-alpa*alpb*(ca-cb)*(ca-cb)
     $                  -alpp(i)*(alpa*ca*ca+alpb*cb*cb))/alpha
                  exparg=exp(arg)*dp(i)/sqrt(alpha)
c                 determine limit for radial integrals and skip evaluation 
c                 if less than threshold eps1.
                  rc=(rka+rkb)/(two*alpha)
c                 ??????? check on what to do if rc=1.0
c                 the above should be /sqrt(alpha)?????
                  if (rc.lt.one) nlim=npi
                  if (rc.gt.one) nlim=n+npi
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
c                 ----- pts and wts method on (0,inf) to be used -----
                  exparg=exp(-alpa*ca*ca-alpb*cb*cb)*dp(i)/sqrt(alpha)
               endif
c
               call rzero (qtemp(0,0,0),343)
c
c              ----- branch for special cases -----
               if (rka.eq.zero.and.rkb.eq.zero) then
c                 ----- one center case -----
                  if ((lalo.eq.0).and.(lblo.eq.0)) then
                     call ptprep (npi,n,0,0,0,0,alpha,rka,rkb,argsum)
                     call quadr (0,0,0,0,0,n,qtemp(0,0,0))
                  endif
               else if(rka.eq.zero) then
c                 ----- two center case -----
                  if (lalo.eq.0) then
                     call ptprep (npi,n,0,0,lbhi-1,lbhi,alpha,
     $                            rka,rkb,argsum)
                     call quadr (0,0,lbhi-1,lbhi,0,n,qtemp(0,0,0))
                     if (lbhi.gt.1) then
                        fkb=one/rkb
                        do 60 lb=lbhi,lblo+2,-1
                           do 40 j=1,npts
                              f(j,0,lb-2)=f(j,0,lb)
     $                              +(2*lb-1)*(fkb/pt(j))*f(j,0,lb-1)
                              qtemp(0,0,lb-2)=qtemp(0,0,lb-2)
     $                                       +ptpow(j,0)*f(j,0,lb-2)
   40                      continue
                           do 50 nabc=1,n
                              qtemp(nabc,0,lb-2)=qtemp(nabc,0,lb)
     $                              +(2*lb-1)*fkb*qtemp(nabc-1,0,lb-1)
   50                      continue
   60                   continue
                     endif
                  endif
               else if(rkb.eq.zero) then
                  call ptprep (npi,n,lahi-1,lahi,0,0,alpha,rka,rkb,
     $                         argsum)
                  if (lblo.eq.0) then
                     call quadr (lahi-1,lahi,0,0,0,n,qtemp(0,0,0))
                     if (lahi.gt.1) then
                        fka=one/rka
                        do 100 la=lahi,lalo+2,-1
                           do 80 j=1,npts
                              f(j,la-2,0)=f(j,la,0)
     $                            +(2*la-1)*(fka/pt(j))*f(j,la-1,0)
                              qtemp(0,la-2,0)=qtemp(0,la-2,0)
     $                                       +ptpow(j,0)*f(j,la-2,0)
   80                      continue
                           do 90 nabc=1,n
                              qtemp(nabc,la-2,0)=qtemp(nabc,la,0)
     $                            +(2*la-1)*fka*qtemp(nabc-1,la-1,0)
   90                      continue
  100                   continue
                     endif
                  endif
               else
c                 ----- general three center case -----
                  call ptprep (npi,n,lahi-1,lahi,lbhi-1,lbhi,
     $                         alpha,rka,rkb,argsum)
                  call recurf (lalo,lahi,lblo,lbhi,rka,rkb)
                  call quadr (lalo,lahi,lblo,lbhi,0,0,qtemp(0,0,0))
                  if (n.gt.0) then
                     call quadr (lahi-1,lahi,lbhi-1,lbhi,1,n,
     $                           qtemp(0,0,0))
                     if (lahi.gt.0.or.lbhi.gt.0) then
                        call recur2 (n,lalo,lahi,lblo,lbhi,
     $                               rka,rkb,qtemp(0,0,0))
                     endif
                  endif
               endif
  120          continue
c
               do 130 lama=lalo,lahi
                  do 130 lamb=lblo,lbhi
                     do 130 nabc=0,n
                        qq(nabc,lama,lamb,prim)=qq(nabc,lama,lamb,prim)
     $                                    +qtemp(nabc,lama,lamb)*exparg
  130          continue
c           end projector loop.
  140       continue
c     end primitives loop.
  150 continue
c
c
      return
      end
