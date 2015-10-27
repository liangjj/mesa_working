*deck @(#)lagran.f	5.1  11/6/94
      subroutine lagran(lag,xbcd,sqdm,t1,twopdm,nactiv,nbf,nnpact,
     #                  nactsq,lnxbcd,nabcd,nder,nocc,ncore)
c
c***begin prologue     lagran
c***date written       860819  (yymmdd)
c***revision date      (yymmdd)
c
c***keywords
c***author             lengsfield, byron (brl)
c***source             @(#)lagran.f	5.1   11/6/94
c***purpose            formation of the mcscf lagrangian
c***description
c
c***references
c***routines called
c***end prologue       lagran

      implicit integer (a-z)
c
      real*8 lag(nbf,nocc),xbcd(lnxbcd),sqdm(nactsq,nnpact)
      real*8 t1(nnpact,nnpact),twopdm(nabcd)
      common /io/     inp,iout
c
      call iosys('read real mcscf_mo_2pdm from rwf',nabcd,twopdm,
     #            0,' ')
c
      do 53 ij=1,nabcd
         twopdm(ij)=twopdm(ij)*0.5d+00
   53 continue
c
      ii=0
      do 1 i=1,nnpact
         ii=ii+i
         twopdm(ii)=twopdm(ii)*2.0d+00
    1 continue
c
      call trtosq(t1,twopdm,nnpact,nabcd)
c
      ii=0
      do 3 i=1,nactiv
         ii=ii+i
         do 2 kl=1,nnpact
            t1(ii,kl)=t1(ii,kl)*2.0d+00
    2    continue
    3 continue
c
      do 4 kl=1,nnpact
         call trtosq(sqdm(1,kl),t1(1,kl),nactiv,nnpact)
    4 continue
c
c     ----- form the lagrangian contribution -----
c
      nt=nbf*nactiv*nnpact
      call iosys('read real der_i(xbcd) from dints',
     #            nt,xbcd,(nder-1)*nt,' ')
c
ctemp
cps         write (6,23) der
cps         call matout(lag(1,1,der),nbf,nocc,nbf,nocc,6)
cend
         call apbct(lag(1,ncore+1),xbcd,sqdm,
     #              nbf,nactiv*nnpact,nactiv)
ctemp
cps         if (der.gt.1) go to 10
cps         write (iout,23) der
cps   23    format (//,' ao lagrangian for',i3,/)
cps         call matout(lag(1,1,der),nbf,nocc,nbf,nocc,iout)
cend
   10 continue
c
c
      return
      end
