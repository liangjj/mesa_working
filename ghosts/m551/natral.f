*deck %W%  %G%
      subroutine natral(lfab,f,den,temp,tv,c,ct,evec,eval,
     $                  nsym,nfo,nco,nao,nvo,nbf,nato,ifock,
     $                  opnscf,iobcas)
c
c***begin prologue
c***date written       871022   (yymmdd)
c***revision date      871987   (yymmdd)
c   19 november 19871  bhl      lanl
c   code inserted to compute naot and to skip the read
c   of the 1-particle dm if naot=0
c
c***keywords
c***author             lengsfield, byron (brl)
c***source             %W%   %G%
c
c***purpose
c
c***description
c
c***references
c
c***routines called    (none)
c
c***end prologue
c
c
      implicit real*8(a-h,o-z)
c
      dimension temp(2),c(2),ct(2),evec(2),eval(2),tv(2),
     $     nfo(2),nco(2),nao(2),nvo(2),nbf(2),nato(2),nobt(11),
     $     isymc(11),isyma(11),isymv(11),isymf(11),lden(11),lfab(2),
     $     f(2),den(2),iobcas(*)
c
      common /io/ inp,iout
c
      integer opnscf
c
c
      nvc=1
      nden=1
      nobtot=0
      nbftot=0
      ncot=0
      do 10 i=1,nsym
         nfoi=nfo(i)
         ncoi=nco(i)
         naoi=nao(i)
         nvoi=nvo(i)
         nbfi=nbf(i)
         nobt(i)=nfoi+ncoi+naoi+nvoi
         nobtot=nobtot+nobt(i)
         isymf(i)=nvc
         nvc=nvc+nfoi*nbfi
         isymc(i)=nvc
         nvc=nvc+ncoi*nbfi
         isyma(i)=nvc
         nvc=nvc+naoi*nbfi
         isymv(i)=nvc
         nvc=nvc+nvoi*nbfi
         lden(i)=nden
         nden=nden+(naoi*(naoi+1))/2
         ncot=ncot+ncoi
         nbftot=nbftot+nbfi
 10   continue
c
      do 15 i=1,nobtot
         eval(i)=0.d0
 15   continue
ccc
c     read fock matrix
ccc
      lnfock=nbftot*nbftot
      if(ifock.eq.0) then
         write(iout,16)
 16      format(' ao total fock operator used to rotate orbitals')
         call iosys('read real mcscf_ao_total_fock from rwf',
     $        infock,f,0,' ')
      else
         write(iout,17)
 17      format(' ao core  fock operator used to rotate orbitals')
         call iosys('read real mcscf_ao_core_fock from rwf',
     $        lnfock,f,0,' ')
      endif
ccc
c     read density matrix
ccc
      ia=1
      naot=0
      do 52 i=1,nsym
         naot=naot+nao(i)
         if(nao(i).eq.0)go to 52
         ia=ia+(nao(i)*(nao(i)+1))/2
 52   continue
c
      norb1=ia-1
      if(naot.ne.0) then
      call iosys('read real mcscf_mo_1pdm from rwf',
     $     norb1,den,0,' ')
      endif
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cc                                                           cc
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      inao=1
      ix=1
      do 100 i=1,nsym
         nbfi=nbf(i)
         if(nbfi.eq.0)go to 100
         nfoi=nfo(i)
         if(nfoi.ne.0) then
c            move frozen orbitals
            ntotm=nfoi*nbfi
            call vmove(ct(isymf(i)),c(isymf(i)),ntotm)
            ix=ix+nfoi
         endif
ccc
ccc
         ncoi=nco(i)
         naoi=nao(i)
         if(ncoi+naoi.eq.0)go to 100
         if(ncoi.ne.0) then
ccc
c            rotate core orbitals
ccc
            ncoit=(ncoi*(ncoi+1))/2
            call ebc(temp,f(lfab(i)),c(isymc(i)),nbfi,nbfi,ncoi)
            call ebtc(tv,c(isymc(i)),temp,ncoi,nbfi,ncoi)
            call sqtotr(temp,tv,ncoi,ncoit)
            call givens(ncoi,ncoi,ncoi,temp,tv,eval(ix),evec)
            call ebc(ct(isymc(i)),c(isymc(i)),evec,nbfi,ncoi,ncoi)
ccc
c           reorder the eigenvectors and eigenvalues
ccc
            call reordr(ct(isymc(i)),eval(ix),nbfi,ncoi,c(isymc(i)))
            ix=ix+ncoi
         endif
c
         if(naoi.eq.0)go to 70
         if(nato(i).eq.0)go to 60
         if(opnscf.eq.0) then
c
            write(iout,29)
  29        format(' 1pdm matrix used to rotate active orbitals')
c
ccc
c           zero density matrix elements between invariant orbital sets
ccc
            if(nato(i).ne.-1) then
               call zerodm(den(lden(i)),naoi,iobcas(inao))
               inao=inao+naoi
            else
               write(iout,31)
  31           format(' mcscf density matrix test ...............',/,
     $                  ' invariant sets may be mixed..............',/,
     $                  ' resulting orbitals may not be variational')
            end if
ccc
c     transform active orbitals to natural orbitals
ccc
            call mcfden(den(lden(i)),temp,naoi)
c
            call givens(naoi,naoi,naoi,temp,tv,eval(ix),evec)
            call ebc(ct(isyma(i)),c(isyma(i)),evec,nbfi,naoi,naoi)
ccc
c     reorder the eigenvectors and eigenvalues
ccc
            call reordr(ct(isyma(i)),eval(ix),nbfi,naoi,c(isyma(i)))
c
            do 55 li=1,naoi
               eval(ix+li-1)=-eval(ix+li-1)
 55         continue
c
         else
c
            write(iout,56)
 56         format(' fock matrix  used to rotate active orbitals')
ccc
c     rotate active orbitals
ccc
            naoit=(naoi*(naoi+1))/2
            call ebc(temp,f(lfab(i)),c(isyma(i)),nbfi,nbfi,naoi)
            call ebtc(tv,c(isyma(i)),temp,naoi,nbfi,naoi)
            call sqtotr(temp,tv,naoi,naoit)
            call givens(naoi,naoi,naoi,temp,tv,eval(ix),evec)
            call ebc(ct(isyma(i)),c(isyma(i)),evec,nbfi,naoi,naoi)
ccc
c     reorder the eigenvectors and eigenvalues
ccc
            call reordr(ct(isyma(i)),eval(ix),nbfi,naoi,c(isyma(i)))
c
         end if
c
         ix=ix+naoi
         go to 70
 60      continue
ccc
c     move active orbitals
ccc
         ntotm=naoi*nbfi
         call vmove(ct(isyma(i)),c(isyma(i)),ntotm)
c
c
         ix=ix+naoi
c
 70      continue
         nvoi=nvo(i)
         if(nvoi.ne.0.and.ncot.ne.0) then
ccc
c           rotate virtual orbitals
ccc
            nvoit=(nvoi*(nvoi+1))/2
            call ebc(temp,f(lfab(i)),c(isymv(i)),nbfi,nbfi,nvoi)
            call ebtc(tv,c(isymv(i)),temp,nvoi,nbfi,nvoi)
            call sqtotr(temp,tv,nvoi,nvoit)
            call givens(nvoi,nvoi,nvoi,temp,tv,eval(ix),evec)
            call ebc(ct(isymv(i)),c(isymv(i)),evec,nbfi,nvoi,nvoi)
ccc
c           reorder the eigenvectors and eigenvalues
ccc
            call reordr(ct(isymv(i)),eval(ix),nbfi,nvoi,c(isymv(i)))
cc
cc
            ix=ix+nvoi
         else
cc
c           move virtual orbitals ..fock operator is not computed if ncot=0
ccc

            ntotm=nvoi*nbfi
            call vmove(ct(isymv(i)),c(isymv(i)),ntotm)
c
c
            ix=ix+nvoi
         endif
c
 100  continue
c
c
      return
c
      end
