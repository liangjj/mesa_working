*deck @(#)socig.f	5.1  11/6/94
      subroutine socig(ms,ifiga,ifigb,ifig,nphase,detdet,vnt)
c
c  evaluate the spin-orbit matrix elements for configurations ifig
c  and jfig (which differ by a single excitation and are not in
c  maximum alignment) by evaluating the hamiltonian matrix over the
c  corresponding determinants.
c
      implicit real*8 (a-h,o-z)
c
      common /c4/iex(4),jex(4),jsymi,jsymj,nopeni,nopenj,ndeti,ndetj
      common /c1/ nbf,nel,ntotfg,mxopn,maxdet,ndbgu,mxocc,maxne0,mel
      common /io/ inp,iout
c
      dimension ms(maxdet,2),ifiga(mxocc,maxdet,2),
     1  ifigb(mxocc,maxdet,2),ifig(nbf,2,*),nphase(maxdet,2),
     2  detdet(ndeti,*),vnt(*)
c
      do 400 ii = 1,ndeti
      nbetai=mxocc-ms(ii,1)
      nalfai=nel-nbetai
      do 350 jj = 1,ndetj
      nbetaj=mxocc-ms(jj,2)
      nalfaj=nel-nbetaj
c
c  determine the excitation type
c
      nmarka=0
      do 100 ka=1,nalfaj
      if( ifig( ifiga(ka,jj,2), 1, ii ) .eq. 1 ) go to 100
      nmarka=nmarka+1
      if( nmarka .gt. 1 ) go to 350
  100 continue
      nmarkb=nmarka
      do 110 ka=1,nbetaj
      if( ifig( ifigb(ka,jj,2), 2, ii ) .eq. 1 ) go to 110
      nmarkb=nmarkb+1
      if( nmarkb .gt. 1 ) go to 350
  110 continue
c
c  single excitation
c
c  the excited orbitals were determined in the calling routine
c
c     iex(1)       the excited orbital of configuration i
c     jex(1)       the excited orbital of configuration j
c
c  determine the sign (ncount is the number of sign changes)
c
      ncount = nphase(ii,1) + nphase(jj,2)
      if(nalfai-nalfaj) 120,190,160
c
c  i-beta j-alpha case
c
  120 do 130 ia=1,nbetai
  130 if(ifigb(ia,ii,1) .lt. iex(1)) ncount=ncount+1
      do 140 ia=1,nalfaj
  140 if(ifiga(ia,jj,2) .gt. jex(1)) ncount=ncount+1
      term = vnt(1) + vnt(2)
      go to 260
c
c  i-alpha j-beta case
c
  160 do 170 ia=1,nalfai
  170 if(ifiga(ia,ii,1) .gt. iex(1)) ncount=ncount+1
      do 180 ia=1,nbetaj
  180 if(ifigb(ia,jj,2) .lt. jex(1)) ncount=ncount+1
      term = vnt(1) - vnt(2)
      go to 260
c
c  i-alpha j-alpha and i-beta j-beta cases
c
  190 if(nmarka .eq. 0) ncount=ncount+1
      term = vnt(3)
      if( abs(iex(1)-jex(1)) .le. 1 ) go to 260
      kbs=max(iex(1),jex(1))-1
      ka=iex(1)+jex(1)-kbs
      if( nmarka .ne. 0 ) then
        do 230 ia=1,nalfai
        if( ifiga(ia,ii,1) .gt. kbs ) go to 260
        if( ifiga(ia,ii,1) .ge. ka ) ncount=ncount+1
  230   continue
      else
        do 250 ia=1,nbetai
        if( ifigb(ia,ii,1) .gt. kbs ) go to 260
        if( ifigb(ia,ii,1) .ge. ka ) ncount=ncount+1
  250   continue
      endif
c
  260 if( mod(ncount,2) .eq. 1 ) term = -term
      detdet(ii,jj) = detdet(ii,jj) + term
c
  350 continue
  400 continue
c
      return
      end
