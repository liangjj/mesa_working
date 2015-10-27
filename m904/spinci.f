*deck @(#)spinci.f	5.1  11/6/94
      subroutine spinci(valone,valtwo,valx,valy,valz,joutfg,iopn,
     1  ms,ifiga,ifigb,ifig,isuma,isumb,nphase,idstor,number,
     2  coef,detdet,vnt,row,msymtp)
c
c  for configurations i and j,
c
c       1. determine the excitation number
c       2. determine the excited orbitals
c       3. call socig or sefcig to evaluate blocks of matrix
c            elements over determinants
c       4. call ciout to transform to double-group-adapted
c            functions and write out the matrix elements
c
      implicit real*8 (a-h,o-z)
      character*80 rtitle, ctitle, blabel
c
      common /c4/iex(4),jex(4),jsymi,jsymj,nopeni,nopenj,ndeti,ndetj
      common /c5/ i,j,nsefi,nsefj,isef,jsef,nzero,kmax,kmaxb
      common /detio/ idfrst, idspc, irt, iau, inextu, intmax, ihso
      common/tapes/iw,iunt1a,iunt1b,iunt2a,iunts1,iunts2
      common /c1/ nbf,nel,ntotfg,mxopn,maxdet,ndbgu,mxocc,maxne0,mel
      common /core/ ecore,thresh,ev,mask,nmask,limit1,limit2,intmxo
      common /headng/ rtitle, ctitle, blabel
      common /c3/ ndet,nsef,idspcu,nroots,maxit,nw,maxesc,ianalz
      common /io/ inp,iout
c
      dimension valone(*),valtwo(*),valx(*),valy(*),valz(*),
     1  joutfg(nbf,*),iopn(mxopn,2),ms(maxdet,2),
     2  ifiga(mxocc,maxdet,2),ifigb(mxocc,maxdet,2),ifig(nbf,2,*),
     3  isuma(maxdet,2),isumb(maxdet,2),nphase(maxdet,2),idstor(*),
     4  number(*),coef(*),detdet(*),vnt(*),row(*),ncoef(2)
c
      id(m,n)=((max(m,n))*((max(m,n))-1))/2+min(m,n)
c
      thresh = 1.0d-8
c
      if(mel.eq.1) go to 840
      iijju=(maxdet*ndbgu)
      if(msymtp.ge.2) iijju=iijju+iijju
      do 800 iijj=1,iijju
  800 coef(iijj)=0.0d0
      rthlf=sqrt(0.5d0)
      iijj=0
      if(msymtp.eq.1) then
        do 810 ii=1,ndbgu
        coef(iijj+1)=rthlf
        coef(iijj+2)=-rthlf
        coef(iijj+maxdet+1)=rthlf
        coef(iijj+maxdet+2)=rthlf
        coef(iijj+2*maxdet+1)=rthlf
        coef(iijj+2*maxdet+2)=rthlf
        coef(iijj+3*maxdet+1)=rthlf
        coef(iijj+3*maxdet+2)=-rthlf
  810   iijj=iijj+(4*maxdet+2)
      elseif(msymtp.eq.2) then
        iijju=(maxdet*ndbgu)
        do 820 ii=1,ndbgu/2
        coef(iijj+1)=rthlf
        coef(iijj+2)=rthlf
        coef(iijj+maxdet+1)=rthlf
        coef(iijj+maxdet+2)=-rthlf
        coef(iijju+iijj+1)=rthlf
        coef(iijju+iijj+2)=rthlf
        coef(iijju+iijj+maxdet+1)=rthlf
        coef(iijju+iijj+maxdet+2)=-rthlf
  820   iijj=iijj+(2*maxdet+2)
        ncoef(2)=iijju+1
      else
        iijju=(maxdet*ndbgu)
        do 830 ii=1,ndbgu
        coef(iijj+1)=rthlf
        coef(iijj+2)=rthlf
        coef(iijju+iijj+1)=rthlf
        coef(iijju+iijj+2)=-rthlf
  830   iijj=iijj+(maxdet+2)
        ncoef(2)=iijju+1
      endif
c
  840 ncoef(1)=1
c
      rewind iunts1
      rewind iunt2a
c
      idspc = idfrst
      nzero = 0
      kmax = 1
      kmaxu = 1
      irowu = 2
c
      isef = 1
c
      do 700 i=1,ntotfg
c      write(iw,*) 'doing row ',i
c
      call detin(joutfg,iopn,ms,ifiga,ifigb,isuma,isumb,nphase,
     1  idstor,1,jsymi,nsefi,nopeni,ndeti,ncpi)
      ncpi=ncoef(ncpi)
c
c  zero out the array containing the number of nonzero matrix
c  elements per row
c
      do 100 ii=1,nsefi
  100 number(ii) = 0
c
c  construct the determinants for configuration i.
c
c  ifig(i,1,ii) and ifig(i,2,ii) specify the determinants in
c  configuration i.
c
      do 140 ii=1,ndeti
c
c  determine nalpha and nbeta
c
      nbeta = mxocc - ms(ii,1)
      nalpha  = nel - nbeta
c
c  clear ifig
c
      do 110 kk = 1,2
      do 110 ia = 1,nbf
  110 ifig(ia,kk,ii)=0
c
c  put a 1 in ifig(i,1,ii) if the ith mo is occupied with alpha spin
c
      do 120 ka = 1,nalpha
  120 ifig( ifiga(ka,ii,1), 1, ii ) = 1
c
c  process the list of beta spin orbitals in a similar manner
c
      do 130 ka = 1,nbeta
  130 ifig( ifigb(ka,ii,1), 2, ii ) = 1
c
  140 continue
c
      maxne0 = min((intmax-1)/(nsefi+1),nsef)
      icimat = max(kmax,maxne0) + 2
      maxne0 = min(maxne0,(intmax-icimat+1)/nsefi)
      kmaxb = 0
      kmaxc = kmax
      rewind iunts1
      jsef = 1
      idspc = idfrst
      do 610 j=1,i
c      write(iw,*) 'doing column',j
      call detin(joutfg,iopn,ms,ifiga,ifigb,isuma,isumb,nphase,
     1  idstor,2,jsymj,nsefj,nopenj,ndetj,ncpj)
      ncpj=ncoef(ncpj)
c
c  zero out the detdet hamiltonian matrix element block
c
      iijju=ndeti*ndetj
      do 150 iijj=1,iijju
  150 detdet(iijj) = 0.0d0
c
c  check to see if a diagonal detdet block is being processed
c
      if(i.eq.j) go to 300
c
c  the detdet block is not diagonal
c
c  determine the excitation number and the excited orbitals.  if
c  the excitation number is greater than 2 then the block is zero
c
      now=abs(nopeni-nopenj)
      if( now .gt. 4 ) go to 600
c
      mm=0
      kk=0
      do 190 ia=1,nbf
      jmarr=joutfg(ia,1)-joutfg(ia,2)
      if(jmarr) 160, 190, 170
  160 kk=kk+1
      jex(kk)=ia
      if( jmarr .ne. -2 ) go to 180
      kk=kk+1
      jex(kk)=ia
      go to 180
  170 mm=mm+1
      iex(mm)=ia
      if( jmarr .ne. 2 ) go to 180
      mm=mm+1
      iex(mm)=ia
  180 if( kk .gt. 2 .or. mm .gt. 2 ) go to 600
  190 continue
c
      if( kk .eq. 1 ) go to 250
c
c  double excitation
c
c  the excitations can be classified according to (1) the difference
c  in the number of open shells (now) and (2) the sum of the
c  occupation numbers of the excited orbitals (jmarr) as follows:
c
c                                       excitation type
c
c     now = 0       jmarr = 4           ab = cd
c                         = 8           a(2) = b(2)
c     now = 2       jmarr = 5           a(2) = bc
c                         = 6           a(2)b = acd
c     now = 4       jmarr = 6           a(2)b(2) = abcd
c
c  determine the integrals needed and the alignment
c
      vnt(1)=valtwo(indexf(iex(1),jex(1),iex(2),jex(2)))
      vnt(2)=valtwo(indexf(iex(1),jex(2),iex(2),jex(1)))
      go to 510
c
c  single excitation
c
c  determine vnt(1) and the alignment
c
c     vnt(1) = (iex1,h,jex1) + sum(core)(iex1,2jc-kc,jex1)
c                            + sum(open)(iex1,ji,jex1)
c                            - (iex1,jiex1,jex1)
c
  250 iex1 = iex(1)
      jex1 = jex(1)
      ijex1 = id(iex1,jex1)
c
c  spin-orbit matrix element
c
      vnt(1) = valx(ijex1)
      vnt(2) = valy(ijex1)
      vnt(3) = valz(ijex1)
      if( vnt(1) .eq. 0.0d0 .and. vnt(2) .eq. 0.0d0
     1                      .and. vnt(3) .eq. 0.0d0 ) go to 260
      if( jex1 .gt. iex1 ) then
        vnt(1) = -vnt(1)
        vnt(2) = -vnt(2)
        vnt(3) = -vnt(3)
      endif
      vnt(1) = 0.5d0 * vnt(1)
      vnt(2) = 0.5d0 * vnt(2)
      vnt(3) = 0.5d0 * vnt(3)
      call socig(ms,ifiga,ifigb,ifig,nphase,detdet,vnt)
c
  260 if( jsymi .ne. jsymj ) go to 600
      ijex = indexf(iex1,jex1,iex1,iex1)
      vnt(1) = valone(ijex1) - valtwo(ijex)
      do 270 ia=1,nbf
      if( joutfg(ia,1) .eq. 0 ) go to 270
      ijex = indexf(iex1,jex1,ia,ia)
      vnt(1) = vnt(1) + valtwo(ijex)
      if( joutfg(ia,1) .eq. 2 )  vnt(1) = vnt(1) + valtwo(ijex)
      iaijex = indexf(iex1,ia,ia,jex1)
      if( joutfg(ia,1) .eq. 2 .and. joutfg(ia,2) .eq. 2)
     2        vnt(1) = vnt(1) - valtwo(iaijex)
  270 continue
      go to 510
c
c  diagonal block
c
c  check to see if any kints are needed
c
  300 if( nopeni .le. 1 ) go to 320
c
c  kints obtained from transpositions on the open shells of
c  configuration i.
c
      mm=1
      do 310 ii=2,nopeni
      do 310 jj=1,ii-1
      mm=mm+1
      iopnii = iopn(ii,1)
      iopnij = iopn(jj,1)
  310 vnt(mm) = valtwo(indexf(iopnii, iopnij, iopnij, iopnii))
c
c  determine vnt(1)
c
c     vnt(1) = sum(core+opens)*(ii,h,ii)*occnum(ii)
c                       + sum(core)*sum(open)*(2jci-kci)
c                       + sum(c1.le.c2)*(2jc1c2-kc1c2)
c                       + sum(c1.lt.c2)*(2jc1c2-kc1c2)
c                       + sum(i1.lt.i2)ji1i2
c
  320 vnt(1) = 0.0d0
      do 380 ibf=1,nbf
      if( joutfg(ibf,1) .eq. 0 ) go to 380
      ih = id(ibf,ibf)
      if( joutfg(ibf,1) .eq. 1 ) then
        vnt(1)=vnt(1)+valone(ih)
      else
        idir = indexf(ibf,ibf,ibf,ibf)
        vnt(1)=vnt(1)+valone(ih)+valone(ih)+valtwo(idir)
      endif
      if( ibf .eq. 1 ) go to 380
      jbfu = ibf-1
      do 370 jbf=1,jbfu
      if( joutfg(jbf,1) .eq. 0 ) go to 370
      idir = indexf(ibf,ibf,jbf,jbf)
      ninj = joutfg(ibf,1)*joutfg(jbf,1)
      if( ninj .eq. 1 ) then
        vnt(1)=vnt(1)+valtwo(idir)
      else
        iexch = indexf(ibf,jbf,ibf,jbf)
        ps = valtwo(idir)+valtwo(idir)-valtwo(iexch)
        if( ninj .eq. 2 ) then
          vnt(1)=vnt(1)+ps
        else
          vnt(1)=vnt(1)+ps+ps
        endif
      endif
  370 continue
  380 continue
c
      if(i.eq.1) then
        vnt1 = vnt(1)
        ecore = ecore + vnt1
        write (iunt2a) rtitle, blabel, ecore, thresh
      endif
      vnt(1) = vnt(1) - vnt1
c
  510 call sefcig(valtwo,iopn,ms,ifiga,ifigb,ifig,isuma,isumb,
     2  nphase,detdet,vnt)
c
  600 ids=1
      iss=1
      if(mel.eq.0) then
        if(nopenj.ne.0) ids=ndeti*ndetj+1
        iss=ids
        if(nopeni.ne.0) iss=ids+ndeti*nsefj
      endif
c      write(iw,*) 'before ciout'
c      write(iw,*) 'icimat',icimat
  610 call ciout(coef(ncpi),coef(ncpj),number,detdet,detdet(ids),
     1  detdet(iss),row,row(icimat),msymtp)
      kmaxu = max(kmaxu,kmaxb)
      irowu = max(irowu,max(kmaxc,kmaxb)+nsefi*kmaxb+1)
c
  700 continue
c
      rewind iunt2a
c
      maxne0 = kmaxu + 1
      intmax = irowu
      write (iw,1000) ecore, thresh, nzero, maxne0
 1000 format(//
     1  '0core energy, including diagonal displacement      =',f14.8/
     2  '0threshold for the hamiltonian matrix elements     =',1p,d10.1/
     3  '0number of matrix elements above threshold         =',i10/
     4  '0maximum number of nonzero matrix elements per row =',i10)
c
      return
      end
