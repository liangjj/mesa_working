*deck @(#)sefcig.f	5.1  11/6/94
      subroutine sefcig(valtwo,iopn,ms,ifiga,ifigb,ifig,isuma,isumb,
     2  nphase,detdet,vnt)
c
c  this routine evaluates the kinetic-energy, coulomb-interaction
c  matrix-element block over configurations ifig and jfig (which
c  differ by a single or double excitation and are not in maximum
c  alignment) by evaluating matrix elements over the corresponding
c  determinants.  the determinants must have the same spin
c  projection to have non-zero matrix elements.
c
      implicit real*8 (a-h,o-z)
c
      common /c4/iex(4),jex(4),jsymi,jsymj,nopeni,nopenj,ndeti,ndetj
      common /c1/ nbf,nel,ntotfg,mxopn,maxdet,ndbgu,mxocc,maxne0,mel
      common /io/ inp,iout
c
      dimension nij(4)
      dimension valtwo(*),iopn(mxopn,2),ms(maxdet,2),
     1  ifiga(mxocc,maxdet,2),ifigb(mxocc,maxdet,2),ifig(nbf,2,*),
     2  isuma(maxdet,2),isumb(maxdet,2),nphase(maxdet,2),
     3  detdet(ndeti,*), vnt(*)
c
      do 400 ii = 1,ndeti
      do 350 jj = 1,ndetj
c
c  determine if dets ii and jj have the same spin projection
c
      if(ms(ii,1).ne.ms(jj,2)) go to 350
c
c  calculate nalpha and nbeta for this pair of determinants.
c
      nbeta = mxocc - ms(ii,1)
      nalpha  = nel - nbeta
c
c  determine the excitation type
c
      nmarka=0
      do 100 ka=1,nalpha
      if( ifig( ifiga(ka,jj,2), 1, ii ) .eq. 1 ) go to 100
      nmarka=nmarka+1
      if( nmarka .gt. 2 ) go to 350
      nij(nmarka)=ka
  100 continue
      nmarkb=nmarka
      do 110 ka=1,nbeta
      if( ifig( ifigb(ka,jj,2), 2, ii ) .eq. 1 ) go to 110
      nmarkb=nmarkb+1
      if( nmarkb .gt. 2 ) go to 350
      nij(nmarkb)=ka
  110 continue
c
      if( nmarkb .eq. 1 ) go to 220
      if( nmarka .eq. 1 ) go to 120
      if( nmarkb .eq. 2 ) go to 170
c
c  no excitation
c
      term=vnt(1)
      if(nopeni.le.1) go to 310
      mm=1
      do 115 iopen=2,nopeni
      iorb=iopn(iopen,1)
      do 115 jopen=1,iopen-1
      mm=mm+1
      jorb=iopn(jopen,1)
      if(ifig(iorb,1,ii).eq.ifig(jorb,1,ii)) term=term-vnt(mm)
  115 continue
      go to 310
c
c  double excitation (product of single excitations from the
c  alpha and beta sets)
c
  120 ncount = nphase(ii,1) + nphase(jj,2)
c
c  determine the excited orbitals
c
c     iex(3) and iex(4)      the excited orbitals of the alpha set
c     jex(3) and jex(4)      the excited orbitals of the beta set
c
      nij1 = nij(1)
      nij2 = nij(2)
      iex(3) = ifiga(nij1,jj,2)
      jex(3) = ifigb(nij2,jj,2)
      iex(4)=isuma(ii,1)-isuma(jj,2)+iex(3)
      jex(4)=isumb(ii,1)-isumb(jj,2)+jex(3)
c
c  determine the signs (ncount is the number of sign changes)
c
      if( abs(iex(3)-iex(4)) .le. 1 ) go to 140
      kbs=max(iex(3),iex(4))-1
      ka=iex(3)+iex(4)-kbs
      do 130 ia=1,nalpha
      if( ifiga(ia,ii,1) .gt. kbs ) go to 140
      if( ifiga(ia,ii,1) .ge. ka ) ncount=ncount+1
  130 continue
  140 if( abs(jex(3)-jex(4)) .le. 1 ) go to 160
      kbs=max(jex(3),jex(4))-1
      ka=jex(3)+jex(4)-kbs
      do 150 ia=1,nbeta
      if( ifigb(ia,ii,1) .gt. kbs ) go to 160
      if( ifigb(ia,ii,1) .ge. ka ) ncount=ncount+1
  150 continue
c
  160 term = valtwo(indexf(jex(3),jex(4),iex(3),iex(4)))
      go to 300
c
c  double excitation
c
c  the excited orbitals were determined in the calling routine
c
c     iex(1) and iex(2)      the excited orbitals of configuration i
c     jex(1) and jex(2)      the excited orbitals of configuration j
c
c  determine the sign (ncount is the number of sign changes)
c
  170 ncount = nphase(ii,1) + nphase(jj,2) + nij(2) - nij(1) - 1
      if( (iex(2)-iex(1)) .le. 1 ) go to 210
      if( nmarka .eq. 0 ) go to 190
      do 180 ia=1,nalpha
      if( ifiga(ia,ii,1) .ge. iex(2) ) go to 210
      if( ifiga(ia,ii,1) .gt. iex(1) ) ncount=ncount+1
  180 continue
      go to 210
  190 do 200 ia=1,nbeta
      if( ifigb(ia,ii,1) .ge. iex(2) ) go to 210
      if( ifigb(ia,ii,1) .gt. iex(1) ) ncount=ncount+1
  200 continue
c
  210 term = vnt(1) - vnt(2)
      go to 300
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
  220 ncount = nphase(ii,1) + nphase(jj,2)
      if( abs(iex(1)-jex(1)) .le. 1 ) go to 260
      kbs=max(iex(1),jex(1))-1
      ka=iex(1)+jex(1)-kbs
      if( nmarka .eq. 0 ) go to 240
      do 230 ia=1,nalpha
      if( ifiga(ia,ii,1) .gt. kbs ) go to 260
      if( ifiga(ia,ii,1) .ge. ka ) ncount=ncount+1
  230 continue
      go to 260
  240 do 250 ia=1,nbeta
      if( ifigb(ia,ii,1) .gt. kbs ) go to 260
      if( ifigb(ia,ii,1) .ge. ka ) ncount=ncount+1
  250 continue
c
  260 term = vnt(1)
c
      kk=1
      if( nmarka .eq. 0 ) kk=2
      do 270 ia=1,nbf
      if( ifig(ia,kk,ii) .eq. 0 ) go to 270
      if( ia .eq. iex(1) .or. ia .eq. jex(1) .or.
     1    ifig(ia,1,ii)+ifig(ia,2,ii) .eq. 2 ) go to 270
      term=term-valtwo(indexf(iex(1),ia,ia,jex(1)))
  270 continue
c
  300 if( mod(ncount,2) .eq. 1 ) term = -term
  310 detdet(ii,jj) = detdet(ii,jj) + term
c
  350 continue
  400 continue
c
      return
      end
