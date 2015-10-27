*deck @(#)prjges.f	5.1  11/6/94
      subroutine prjges(gestyp,nbasi,nnp,c,eigval,smhalf,altges,ops)
c***begin prologue     prjges.f
c***date written       850601  yymmdd
c***revision date      11/6/94
c
c   30 june     1993   pjh at lanl
c        removed expges call for rdinp case
c        cannot simultaneously use guess=rdinp and expand now
c   12 march    1991   rlm at lanl
c        removing some benign errors associated with use of
c        intowp.
c    8 april    1988   rlm at lanl
c        changing test on top and maxcor prior to getscm call.
c    1 december 1986   pws at lanl
c        changing 'namchk' and iosys open to character
c
c***keywords           initial guess, projection
c***author             martin, richard (lanl)
c***source             @(#)prjges.f	5.1   11/6/94
c***purpose            driver for guesses which require projection
c                      from one basis set to another.
c***description
c
c***references         (none)
c***routines called    iosys(io)
c                      rdbas(m412), rdmo(m412), schmdt(math),
c                      minbas(m412), vmove(math), ovrlap(m412), sqtotr(math),
c                      sinv(util), hukmat(m412), alter(m412), ebc(math),
c                      trtosq(math), ebtc(math), ebct(math), sizmb(m412),
c                      smul(math), givens(math)
c***end prologue       prjges.f
      implicit integer(a-z)
      real*8 c(nbasi,nbasi),eigval(nbasi),smhalf(nnp)
      real*8 maxerr
      real*8 z, y
      character*(*) gestyp,ops
      character*128 namchk
      logical altges,rveceq,iveceq,logkey
      logical debug
      dimension need(2), ngot(2)
c
      parameter (debug=.false.)
c
      common /io/ inp,iout
      pointer (py,y(1)), (py,ay(1))
      pointer (pz,z(1)), (pz,az(1))
c
 1000 format(5x,'original huckel guess orbitals')
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
      if(gestyp.eq.'rdchk') then
c         base=1
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
      else if(gestyp.eq.'rdinp') then
c         base=1
         natj=nati
         nbasj=nbasi
         nprimj=nprimi
         ncontj=nconti
      else if(gestyp.eq.'huckel') then
         natj=nati
         atomno=1
         atomz=iadtwp(atomno+natj)
         need(1)=wpadti(atomz+natj)
         call getmem(need(1),py,ngot(1),'huckel',0)
c         base=atomz+natj
         call iosys('read integer "atomic numbers" from rwf',
     $        -1,ay(atomno),0,' ')
         call iosys('read real "nuclear charges" from rwf',
     $               -1,y(atomz),0,' ')
         call sizmb(natj,nprimj,ncontj,nbasj,ay(atomno),y(atomz))
      endif
c
c     ----- divide core for basis set information -----
c
c
c     ci=base
      ci=1
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
      oldmo=smo+nbsq
      oldeig=oldmo+nbsq
      triang=oldeig+nbasj
      scr=triang+max(nnp,nnj)
      eigvec=scr+max(5*nbasi,nbsq)
      smbhlf=eigvec+max(5*nbasi,nbsq)
      hmat=smbhlf+nnj
      u=hmat+nnj
      u0=u+nbsq
      s0ph=u0+nbsq
      smh=s0ph+nbsq
      t1=smh+nbsq
      t2=t1+nbsq
      t3=t2+nbsq
      t4=t3+nbsq
      need(2)=wpadti(t4+nbsq+30000)
      top=need(2)
      call getmem(need(2),pz,ngot(2),'prjges',0)
c
c     read in current basis set information.
      call rdbas('rwf',z(exi),z(conti),z(ci),az(ptprmi),az(noprmi),
     $            az(ptcnti),az(nocnti),az(strti),az(nctype),az(nocart),
     $            az(nobf),az(maxmom),az(minmom),az(mintyp),
     $            az(nx),az(ny),az(nz),nati,nprimi,nconti,ntypes,lenxyz)
c
c     generate guess vectors.
      if(gestyp.eq.'rdchk') then
         call rdbas('chk',z(exj),z(contj),z(cj),az(ptprmj),az(noprmj),
     $               az(ptcntj),az(nocntj),az(strtj),az(nctype),
     $               az(nocart),az(nobf),az(maxmom),az(minmom),
     $               az(mintyp),az(nx),az(ny),az(nz),natj,nprimj,
     $               ncontj,ntypes,lenxyz)
c
c        ----- transform old vector to new coordinates, a la
c              louis carlacci and james w. mciver, jr.
c              j. chem. phys. 85, 634 (1986)
c              'a rotationally invariant orbital transformation'
c
c
c        ----- compute the overlap matrix with old basis at old geometry
c
         call ovrlap(z(cj),z(cj),z(exj),z(exj),z(contj),z(contj),
     $               z(s),az(ptprmj),az(ptprmj),az(noprmj),az(noprmj),
     $               az(nocntj),az(nocntj),az(ptcntj),az(ptcntj),
     $               natj,natj,nprimj,nprimj,ntypes,nbtype,
     $               ncontj,ncontj,az(strtj),az(strtj),nbasj,nbasj,
     $               az(nocart),az(nobf),az(maxmom),az(minmom),
     $               az(mintyp),az(nx),az(ny),az(nz),lenxyz)
c
c        ----- get s-1/2 for basis. eigvec holds eigenvectors of s -----
c
         call sqtotr(z(triang),z(s),nbasj,nnj)
         iprint=0
         call sinv(z(triang),z(smbhlf),z(u0),z(oldeig),z(smo),
     $             z(oldmo),nbasj,nnj,z(scr),iprint)
c         write (iout,3366)
c 3366    format ('old s and u')
c         call matout(z(s),nbasj,nbasj,nbasj,nbasj,iout)
c         call matout(z(u0),nbasj,nbasj,nbasj,nbasj,iout)
c
c        ----- now multiply by s to get s+1/2. -----
c
         call trtosq(z(scr),z(smbhlf),nbasj,nnj)
         call ebc(z(s0ph),z(s),z(scr),nbasj,nbasj,nbasj)
c
c        ----- get the old mo's. -----
c
         call rdmo('chk',nbasj,nummo,z(oldmo),z(oldeig),ops)
c
c        ----- close the checkpoint file. -----
c
         call iosys('close chk',namchk,0,0,' ')
c
c        ----- calculate the overlap matrix and s**-1/2 at new geometry
c
         call ovrlap(z(ci),z(ci),z(exj),z(exj),z(contj),z(contj),
     $               z(s),az(ptprmj),az(ptprmj),az(noprmj),az(noprmj),
     $               az(nocntj),az(nocntj),az(ptcntj),az(ptcntj),
     $               natj,natj,nprimj,nprimj,ntypes,nbtype,
     $               ncontj,ncontj,az(strtj),az(strtj),nbasj,nbasj,
     $               az(nocart),az(nobf),az(maxmom),az(minmom),
     $               az(mintyp),az(nx),az(ny),az(nz),lenxyz)
c
c        ----- get s-1/2 for basis. eigvec holds eigenvectors of s -----
c
         call sqtotr(z(triang),z(s),nbasj,nnj)
         iprint=0
         call sinv(z(triang),z(smbhlf),z(u),z(oldeig),z(smo),
     $             z(t1),nbasj,nnj,z(scr),iprint)
         call trtosq(z(smh),z(smbhlf),nbasj,nnj)
c         write (iout,4466)
c 4466    format ('new s and u')
c         call matout(z(s),nbasj,nbasj,nbasj,nbasj,iout)
c         call matout(z(u),nbasj,nbasj,nbasj,nbasj,iout)
c
c        ----- make sure that u --> u0 -----
c
         call fixup(nbasj,z(u),z(u0),z(t1),z(t2))
c
c        ----- form u . u0(t) -----
c
         call ebct(z(t1),z(u),z(u0),nbasj,nbasj,nbasj)
c         write (iout,9922)
c 9922    format (' u u dagger')
c         call matout(z(t1),nbasj,nbasj,nbasj,nbasj,iout)
c
c        ----- form u . u(t)s**1/2 -----
c
         call ebc(z(t2),z(t1),z(s0ph),nbasj,nbasj,nbasj)
c
c        ----- and s**-1/2 . uu(t)s**1/2 -----
c
         call ebc(z(t3),z(smh),z(t2),nbasj,nbasj,nbasj)
c 9933    format (' suus')
c         write (iout,9933)
c         call matout(z(t3),nbasj,nbasj,nbasj,nbasj,iout)
c
c        ----- and c = s**-1/2 . u . u0(t) . s0**1/2 . c0 -----
c
         call ebc(z(t4),z(t3),z(oldmo),nbasj,nbasj,nbasj)
c         write (iout,9944)
c 9944    format (' new vector:')
c         call matout(z(t4),nbasj,nbasj,nbasj,nbasj,iout)
c
         call vmove(z(oldmo),z(t4),nbasj**2)
c
c        ----- check if the two basis are identical, in which case
c              transfer old to new and be done. this avoids the
c              arbitrary rotations invoked by the diagonalization
c
         if (nconti.eq.ncontj.and.nprimi.eq.nprimj.and.nati.eq.natj)
     #                                                   then
            if (rveceq(z(exi),z(exj),nprimi).and.
     #          rveceq(z(conti),z(contj),nconti).and.
     #          iveceq(az(ptprmi),az(ptprmj),ntypes*nati).and.
     #          iveceq(az(noprmi),az(noprmj),ntypes*nati).and.
     #          iveceq(az(nocnti),az(nocntj),ntypes*nati).and.
     #          iveceq(az(ptcnti),az(ptcntj),ntypes*nati)) then
c
cps               call rdmo('chk',nbasj,nummo,c,eigval,ops)
               call vmove(eigval,z(oldeig),nbasj)
               call vmove(c,z(oldmo),nbasj**2)
               if (altges) call alter(c,eigval,z(triang),nbasj)
c
               return
            end if
         end if
c
      else if(gestyp.eq.'rdinp') then
         call rdbas('rwf',z(exj),z(contj),z(cj),az(ptprmj),az(noprmj),
     $               az(ptcntj),az(nocntj),az(strtj),az(nctype),
     $               az(nocart),az(nobf),az(maxmom),az(minmom),
     $               az(mintyp),az(nx),az(ny),az(nz),natj,nprimj,
     $               ncontj,ntypes,lenxyz)
c
c        get s-1/2 for basis.
         call iosys('read real "overlap integrals" from rwf',
     $               nnp,z(triang),0,' ')
         iprint=0
         call sinv(z(triang),z(smbhlf),z(oldmo),z(oldeig),z(smo),
     $             z(eigvec),nbasj,nnj,z(scr),iprint)
c        now multiply by s to get s+1/2.
         call trtosq(z(scr),z(triang),nbasj,nnj)
         call trtosq(z(oldmo),z(smbhlf),nbasj,nbasj)
         call ebc(z(eigvec),z(scr),z(oldmo),nbasj,nbasj,nbasj)
         call sqtotr(z(smbhlf),z(eigvec),nbasj,nnj)
c        get the input mo's.
         call rdmo('inp',nbasj,nae,z(oldmo),z(oldeig),ops)
         maxerr=1.0d-06
         call schmdt(z(oldmo),z(triang),z(smo),z(scr),z(eigvec),
     $               nae,nbasj,nnj,maxerr)
c        call wvec(z(oldmo),z(oldeig),nbasj,nae,' ',' ')
         call vmove(c,z(oldmo),nbasj**2)
c

c
      else if(gestyp.eq.'huckel') then
c
c        generate the minimum basis set, huckel diagonal elements,
c        and a flag for the core orbitals.
         call minbas(natj,z(exj),z(contj),az(ptprmj),az(noprmj),
     $               az(ptcntj),az(nocntj),az(strtj),az(pstart),
     $               az(nctype),ay(atomno),y(atomz),az(nocart),az(nobf),
     $               az(maxmom),az(minmom),az(mintyp),az(nx),az(ny),
     $               az(nz),natj,nprimj,ncontj,nbasj,ntypes,nbtype,
     $               lenxyz,z(ediag),az(corflg))
         call getmem(-ngot(1),py,idum,'huckel',idum)
c
c        compute the minimum basis overlap matrix.
c
         call vmove(z(cj),z(ci),3*nati)
         call ovrlap(z(cj),z(cj),z(exj),z(exj),z(contj),z(contj),
     $               z(s),az(ptprmj),az(ptprmj),az(noprmj),az(noprmj),
     $               az(nocntj),az(nocntj),az(ptcntj),az(ptcntj),
     $               natj,natj,nprimj,nprimj,ntypes,nbtype,ncontj,
     $               ncontj,az(strtj),az(strtj),nbasj,nbasj,az(nocart),
     $               az(nobf),az(maxmom),az(minmom),az(mintyp),
     $               az(nx),az(ny),az(nz),lenxyz)
c
c        get s-1/2 for basis.
         iprint=0
         call sqtotr(z(scr),z(s),nbasj,nnj)
         call sinv(z(scr),z(smbhlf),z(oldmo),z(oldeig),z(smo),
     $             z(eigvec),nbasj,nnj,z(triang),iprint)
c
c        get the huckel vectors.
         call hukmat(z(s),z(smbhlf),z(ediag),az(corflg),z(hmat),
     $               z(oldmo),z(oldeig),z(triang),z(smo),z(scr),
     $               nbasj,nnj)
c
c        --- print the original huckel guess
         if(logkey(ops,'print=guess',.false.,' ')) then
            if(logkey(ops,'print=guess=all',.false.,' ')) then
               write(iout,1000)
               call matprt(z(oldmo),nbasj,nbasj,nbasj,nbasj,0,0,
     $                     ' ',' ',0,z(oldeig),.true.)
            else
               numprt=min(nbasj,nae+5)
               write(iout,1000)
               call matprt(z(oldmo),nbasj,nbasj,nbasj,numprt,0,0,
     $                     ' ',' ',0,z(oldeig),.true.)
            endif
         endif
c        now multiply s-1/2 by s to get s+1/2.
         call trtosq(z(scr),z(smbhlf),nbasj,nnj)
         call ebc(z(eigvec),z(s),z(scr),nbasj,nbasj,nbasj)
         call sqtotr(z(smbhlf),z(eigvec),nbasj,nnj)
      endif
c
c     possibly alter the guess.
      if(altges) call alter(z(oldmo),z(oldeig),z(triang),nbasj)
c
c     compute the overlap between the two bases.
c     note that for the purposes of the projection we effectively
c     translate the current basis to the old coordinates.
c
      if(gestyp.ne.'rdinp') then
         call ovrlap(z(cj),z(cj),z(exi),z(exj),z(conti),z(contj),
     $               z(s),az(ptprmi),az(ptprmj),az(noprmi),az(noprmj),
     $               az(nocnti),az(nocntj),az(ptcnti),az(ptcntj),nati,
     $               natj,nprimi,nprimj,ntypes,nbtype,nconti,ncontj,
     $               az(strti),az(strtj),nbasi,nbasj,az(nocart),
     $               az(nobf),az(maxmom),az(minmom),az(mintyp),
     $               az(nx),az(ny),az(nz),lenxyz)
c
c        now form the overlap between the current basis and the
c        (possibly) altered guess.
c        call ebc(z(scr),z(s),z(oldmo),nbasi,nbasj,nbasj)
c        now transform to the orthonormal basis.
c        call trtosq(z(eigvec),z(smbhlf),nbasj,nnj)
c        call ebc(z(smo),z(scr),z(eigvec),nbasi,nbasj,nbasj)
c
c        form and diagonalize ss(dagger) corresponding to the
c        occupied portion of the guess.
c        call ebct(z(scr),z(smo),z(smo),nbasi,nae,nbasi)
c        call sqtotr(z(triang),z(scr),nbasi,nnp)
c        call trtosq(z(eigvec),z(smbhlf),nbasj,nnj)
c        call ebc(z(scr),z(s),z(eigvec),nbasi,nbasj,nbasj)

         call trtosq(z(eigvec),smhalf,nbasi,nnp)
         call ebc(z(scr),z(eigvec),z(s),nbasi,nbasi,nbasj)
c
         call ebc(z(smo),z(scr),z(oldmo),nbasi,nbasj,nbasj)
         call ebct(z(scr),z(smo),z(smo),nbasi,nae,nbasi)
         call sqtotr(z(triang),z(scr),nbasi,nnp)
         if(debug) then
            write(iout,*) 'scr'
            call matout(z(scr),nbasi,nbasi,nbasi,nbasi,iout)
         endif
c
c
c        diagonalize -s*s(dagger) so eigenvalues bounded by (-1,0).
         call smul(z(triang),z(triang),-1.0d0,nnp)
         call rsp(nbasi,nbasi,nnp,z(triang),eigval,+1,z(eigvec),
     $            z(scr),z(scr+nbasi),ierr)
         if(debug) then
            write(iout,*) 'eigenvalues',eigval
            write(iout,*) 'eigvec',ierr
            call matout(z(eigvec),nbasi,nbasi,nbasi,nbasi,iout)
         endif
c
c        back transform the eigenvectors.
         call trtosq(z(smo),smhalf,nbasi,nnp)
         call ebtc(c,z(smo),z(eigvec),nbasi,nbasi,nbasi)
         if(debug) then
            write(iout,*) 'prjges c'
            call matout(c,nbasi,nbasi,nbasi,nbasi,iout)
         endif
      endif
c
c
      call getmem(-ngot(2),pz,idum,'prjges',idum)
      return
      end
