*deck @(#)prjovr.f	5.1  11/6/94
      subroutine prjovr(gestyp,nbasi,nnp,c,eigval,smhalf,a,z,maxcor,
     $                  altges,ops)
c***begin prologue     prjges
c***date written       850601  yymmdd
c***revision date      880408  yymmdd
c
c    8 april    1988   rlm at lanl
c        changing test on top and maxcor prior to getscm call.
c    1 december 1986   pws at lanl
c        changing 'namchk' and iosys open to character
c
c***keywords           orbital overlap
c***author             martin, richard (lanl); braunstein, matt (lanl)
c***source             @(#)prjovr.f	5.1   11/6/94
c***purpose           
c***description
c
c***references         (none)
c***routines called    iosys(io)
c                      rdbas(m401), rdmo(m401) 
c                      minbas(m401), vmove(math), ovrlap(m401), sqtotr(math),
c                      sinv(util), hukmat(m401), alter(m401), ebc(math),
c                      trtosq(math), ebtc(math), ebct(math), sizmb(m401),
c                      smul(math), givens(math)
c***end prologue       prjges
      implicit integer(a-z)
      real*8 z(maxcor),c(nbasi,nbasi),eigval(nbasi),smhalf(nnp)
      character*(*) gestyp,ops
      character*128 namchk
      logical altges
      logical logkey
      dimension a(*)
c
      common /io/ inp,iout
c
c     note that a,z are equivalenced in the calling routine.
c
c
c     retrieve the current basis parameters from the rwf.
      call iosys('read integer "number of atoms" from rwf',1,nati,0,' ')
      call iosys('length of exponents on rwf',nprimi,0,0,' ')
      call iosys('length of "contraction coefficients" on rwf',
     $     nconti,0,0,' ')
      call iosys('length of "number of pure functions" on rwf',
     $     ntypes,0,0,' ')
      call iosys('read integer "number basis types" from rwf',
     $     1,nbtype,0,' ')
      call iosys('length of "power of x" on rwf',lenxyz,0,0,' ')
      call iosys('read integer "number of alpha electrons" from rwf',
     $     1,nae,0,' ')
c
c
c
      if(gestyp.eq.'order') then
         base=1
c
c        prepare the checkpoint file for action.
         call iosys('read character "checkpoint filename" from rwf',
     $        0,0,0,namchk)
         call iosys('open chk as old',0,0,0,namchk)
c
c        retrieve the old basis parameters from the chk.
c
         call iosys('read integer "number of atoms" from chk',
     $        1,natj,0,' ')
         call iosys('read integer "number of basis functions" from chk',
     $        1,nbasj,0,' ')
         call iosys('length of exponents on chk',nprimj,0,0,' ')
         call iosys('length of "contraction coefficients" on chk',
     $        ncontj,0,0,' ')
c
c        ----- the lengths are returned in integer words -----
c
      endif
c
c     ----- divide core for basis set information -----
c
c
      ci=base
      cj=ci+3*nati
      conti=cj+3*natj
      contj=conti+nconti
      exi=contj+ncontj
      exj=exi+nprimi
      ediag=exj+nprimj
      ptprmi=wpadti(ediag+nbasj)
      ptprmj=ptprmi+ntypes*nati
      noprmi=ptprmj+ntypes*natj
      noprmj=noprmi+ntypes*nati
      nocnti=noprmj+ntypes*natj
      nocntj=nocnti+ntypes*nati
      ptcnti=nocntj+ntypes*natj
      ptcntj=ptcnti+ntypes*nati
      strti=ptcntj+ntypes*natj
      strtj=strti+ntypes*nati
      pstart=strtj+ntypes*natj
      corflg=pstart+ntypes*natj
      nctype=corflg+nbasj
      nocart=nctype+ntypes
      nobf=nocart+ntypes
      minmom=nobf+ntypes
      maxmom=minmom+ntypes
      mintyp=maxmom+ntypes
      nx=mintyp+ntypes
      ny=nx+lenxyz
      nz=ny+lenxyz
c
      nnj=nbasj*(nbasj+1)/2
      nbsq=max(nbasi*nbasi,nbasj*nbasj)
      s=iadtwp(nz+lenxyz)
      smo=s+nbsq
      newmo=smo+nbsq
      oldmo=newmo+nbsq
      oldeig=oldmo+nbsq
      triang=oldeig+nbasj
      t1=triang+max(nnp,nnj)
      t2=t1+nbsq
      xorbs=wpadti(t2+nbsq)
      top=xorbs+nbasj
      top=iadtwp(top)
c
      if((top+30000).gt.maxcor) then
         call getscm(top+30000,z(1),ngot,'chkges',0)
         maxcor=maxcor+iadtwp(ngot)
      endif
      nleft=maxcor-top+1
c
c     read in old basis set information.
      call rdbas('chk',z(exi),z(conti),z(ci),a(ptprmi),a(noprmi),
     $            a(ptcnti),a(nocnti),a(strti),a(nctype),a(nocart),
     $            a(nobf),a(maxmom),a(minmom),a(mintyp),a(nx),a(ny),
     $            a(nz),nati,nprimi,nconti,ntypes,lenxyz)
c
c read in current basis set information
         call rdbas('rwf',z(exj),z(contj),z(cj),a(ptprmj),a(noprmj),
     $            a(ptcntj),a(nocntj),a(strtj),a(nctype),a(nocart),
     $            a(nobf),a(maxmom),a(minmom),a(mintyp),a(nx),a(ny),
     $            a(nz),natj,nprimj,ncontj,ntypes,lenxyz)
c
c        ----- compute the overlap matrix with old basis at old geometry
c
         call ovrlap(z(cj),z(ci),z(exj),z(exi),z(top),z(contj),z(conti),
     $            z(s),a(ptprmj),a(ptprmi),a(noprmj),a(noprmi),
     $            a(nocntj),a(nocnti),a(ptcntj),a(ptcnti),natj,nati,
     $            nprimj,nprimi,nleft,ntypes,nbtype,ncontj,nconti,
     $            a(strtj),a(strti),nbasj,nbasi,a(nocart),a(nobf),
     $            a(maxmom),a(minmom),a(mintyp),a(nx),a(ny),a(nz),
     $            lenxyz)
c
       if(logkey(ops,'print=contracted-overlaps',.false.,' ')) then
          write(iout,*)'contracted basis overlaps'
          call matout(z(s),nbasj,nbasi,nbasj,nbasi,iout)
       end if
c
c
c        ----- get the old mo's. -----
c
       call iosys('read real "mcscf vector" from chk',-1,z(oldmo),0,' ')
c
c        ----- get the new mo's. -----
c
       call iosys('read real "mcscf vector" from rwf',-1,z(newmo),0,' ')
c        ----- close the checkpoint file. -----
c
         call iosys('close chk',namchk,0,0,' ')
c
c matrix multiplication
c
c C^t (newmo) O(basis) C (oldmo)
c
      call ebc(z(t1),z(s),z(oldmo),nbasj,nbasi,nbasi)
      call ebtc(z(smo),z(newmo),z(t1),nbasj,nbasj,nbasi)
c
c    print overlaps
c    
       if(logkey(ops,'print=molecular-overlaps',.false.,' ')) then
          write(iout,*)'molecular orbital overlaps'
          call matout(z(smo),nbasj,nbasi,nbasj,nbasi,iout)
       end if
c
c store diagonal elements of overlap matrix
c
       call intarr(ops,'overlap=order',a(xorbs),nbasj,' ')
       write(iout,*)'diagonal elements of overlap matrix'
       do 10 i=a(xorbs),a(xorbs+1),1
          diag=smo+(i-1)*(nbasj+1)
          write(iout,*)i,z(diag)
10     continue
c
c check for pattern | 0 1 |
c                   | 1 0 |
c
c read relevant orbital list
c
       do 20 i=a(xorbs),a(xorbs+1),1
         diag=smo+(i-1)*(nbasj+1)
         diagp1=diag+nbasj+1
         diagl=diag+1
         diags=diagp1-1
         if(abs(z(diag)).le.0.9d0.and.abs(z(diagl)).gt.0.9d0
     >     .and.abs(z(diags)).gt.0.9d0
     >     .and.abs(z(diagp1)).le.0.9d0) then
c
c switch orbitals i and i+1 from rwf
c
        iswitch=(i-1)*nbasj+newmo
        jswitch=iswitch+nbasj
        call vmove(z(t1),z(iswitch),nbasj)
        call vmove(z(iswitch),z(jswitch),nbasj)
        call vmove(z(jswitch),z(t1),nbasj)
c
c write switched mo's to rwf file
c
         call iosys('write real "mcscf vector" to rwf',-1,
     > z(newmo),0,' ')
c
c switch these diagonal elements of the overlap matrix
c
         z(diag)=z(diagl)
         z(diagp1)=z(diags)
c
c message of which orbitals were switched
c
       write(iout,*)'switching orbitals ',i,' and ',i+1
       go to 25
       end if
20     continue
25     continue
c
c warning if we do not have this pattern
c
       do 30 i=a(xorbs),a(xorbs+1),1
         diag=smo+(i-1)*(nbasj+1)
         if(abs(z(diag)).le.0.9d0)then
           write(iout,*)'warning - unrecognized pattern for overlaps'
           call lnkerr
         else
           write(iout,*)'overlaps o.k.'
         end if
30     continue
c
c change phase of orbital
c
       do 40 i=1,nbasj
         diag=smo+(i-1)*(nbasj+1)
         if (z(diag).lt.0.9d0.and.abs(z(diag)).gt.0.9d0) then
            icount=(i-1)*nbasj+newmo
            call smul(z(t1),z(icount),-1.d0,nbasj)
            call vmove(z(icount),z(t1),nbasj)
            write(iout,*)'changing phase of orbital ',i
         end if
40      continue
        call iosys('write real "mcscf vector" to rwf',-1,
     >      z(newmo),0,' ')
c
c
      return
      end
