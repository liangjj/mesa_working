*deck  @(#)m921.f	1.6 7/30/91
      program m921
c
c***begin prologue     m921
c***date written       yymmdd   (yymmdd)
c***revision date      891127   (yymmdd)
c
c  27 november 1989   rlm at lanl
c     removing hubbard options; fixing up  problem with transition
c     density matrices and the kohn option. adding flags to discern
c     whether m903 was used, in which case the standard order of
c     the orbitals is preserved. mcscf/ci is determined by the
c     flag "mcorci".  this is set in m551 or various ci links
c     (m901,m902,m903,m914). the use of m903 is determined
c     by either "mcscf: ci used", or "ci used".
c
c  19 january  1989   bhl at llnl
c     only call the routine enforce if ci=no=enforce
c
c  02 february 1988   bhl at brl
c     read mcscf unconverged vectors if mcscf failed
c
c  11 august 1987     pws at lanl
c     1. fixed bug causing frozen core density to be 0 for diagonal
c     densities when also running transition density matrices.
c     2. made cutoff an input option, and set so if small all
c     configurations can be searched.
c
c  10 august 1987     pws at lanl
c     fixing up so don't try more than the number-of-configurations
c     roots.
c
c  25 february 1987   pws at lanl
c     adding transition density capabilities.
c
c   2 february 1987   rlm at lanl
c     modifying to handle average natural orbitals with associated
c     weights.
c
c   28 august 1986  pws at lanl
c     don't do natural orbitals for mcscf case, since they have
c     already been done in mcscf
c
c
c***keywords           one-particle density matrix, most important
c                      configurations, ci analysis
c***author             saxe, paul (lanl)
c***source             @(#)pm921.f	1.6 7/30/91
c
c***purpose            to print the most important configurations in a
c                      ci wavefunction and form the one-particle density
c                      matrix.
c
c***description
c       m921 will analyze the ci vectors from one or more
c       roots of a ci, as well as form and print the approximate
c       natural orbitals for each root. currently the following options
c       are recognized:
c
c       nroots=n         the number of roots to analyze.
c       nimportant=n     the number of important configurations to
c                        print, with a default of 15.
c       starting_row=n   the row in the drt to use as the top of the
c                        graph. this is used to pick out a smaller ci
c                        from a larger list.
c       print_nos        print the natural orbitals.
c       print_occ        print only the occupancies of the natural
c                        orbitals.
c       average=(n,m..)  do average natural orbitals from the densities
c                        of roots n,m,...
c       transition=(n,m) calculate one-particle transition density matrix
c                        between roots n and m.
c
c          at the minute, this link requires the following constants
c       from the iosys routine:
c
c       symmetries in ci,  norbs,  nbasis,  nrows,  nwks,  nlevs,
c       orbfrm
c
c       and requires the following files on the read-write file:
c
c       arc,  weight,  orbtbf,  iout,  bfsym,  nlwks,  b,  orbsym,
c       "xform vector",  "ci root n"
c
c       where n is the numbers of the ci roots to analyze. this link
c       places the following files on the read-write file:
c
c       "mo 1pdm n", "no vector n", "no occ n"
c
c       if average natural orbitals are evaluated, they are denoted
c       as vector 0.
c
c
c***references
c
c***routines called    (none)
c
c***end prologue       pm921
c
      implicit integer (a-z)
c
      parameter (maxnbf=2000)
      real*8 one
      real*8 z,cutoff
      real*8 fpkey,wgt
      character*4096 ops,itoc*4
      character*128 nmcnfg, nmfile, namchk
      character*16 bflabl(maxnbf)
      character*3 symlab(8), ans
      character*4 vn
      character*8 citype, scatyp, filtyp, dsk
      character refops*80, chrkey*16
      character xform*32
      character*16 key
      logical doavg,avgall,logkey,trans,scat
      logical ci,mcscf,m903
      logical iskey
      integer a
      pointer(p,z(1)),(p,a(1))
c
      common /io/ inp,iout
c
      parameter (one=1.0d+00)
c
      data symlab /'a','b','c','d','e','f','g','h'/
      data refops/'no average print'/
      data mcscf/.false./, ci/.false./, m903/.false./
      data scat/.false./
c
 1001 format (//,t5,'transition density between roots',i4,' and',i4)
 1002 format (/,t5,'average natural orbitals:')
c
c     ----- recover the options string -----
c
      call drum
      write(iout,*)' m921:wfn and density codes'
c
      call iosys('read character options from rwf',-1,0,0,ops)
c
c     ----- determine which file is to be used to read
c           the configuration and CI data.
c
      key=chrkey(ops,'int=drt=key','drt',' ')
      call pakstr(key,lenkey)
      call iosys('read character "checkpoint filename" from rwf',
     $            0,0,0,namchk)
      call iosys('open chk as old',0,0,0,namchk)
      call iosys('read character "drt file name '//key(1:lenkey)
     #            //'" from chk',0,0,0,dsk)
      call iosys('read character "drt unit file name '//dsk
     #               //'" from chk',0,0,0,nmcnfg)
      call iosys('open '//dsk//' as unknown',0,0,0,nmcnfg)
      write(iout,*) '          reading information from '//dsk
      call iosys('close chk',0,0,0,namchk)
c
c ---- where do we put the scattering information.
c
      scatyp=chrkey(ops,'scattering','none',' ')
      if(scatyp(1:4).eq.'kohn') then
         call iosys('read character "kohn filename" from rwf',
     #               -1,0,0,nmfile)
         filtyp='kohn'
         filtyp=filtyp(1:4)
         scat=.true.
      else if(scatyp(1:8).eq.'r-matrix') then
         call iosys('read character "r-matrix filename" from rwf',
     #               -1,0,0,nmfile)
         filtyp='rmtrx'
         filtyp=filtyp(1:5)
         scat=.true.
      endif
c
c     what kind of wavefunction is being processed?
      call iosys('read character mcorci from rwf',-1,0,0,vn)
      call pakstr(vn,ilen)
      if(vn(1:ilen).eq.'mc') then
         mcscf=.true.
         call iosys('read character "mcscf: ci used" from rwf',
     $               0,0,0,citype)
c        check for mcscf running an scf, in which case exit.
         if(citype.eq.'m903') m903=.true.
         if(citype.eq.'hf') go to 5000
      else
         ci=.true.
         call iosys('read character "ci used" from rwf',0,0,0,citype)
         if(citype.eq.'m903') m903=.true.
      endif
c
c     ----- get the constants for dividing up core -----
c           based on the information on dsk.
c
      call iosys('read integer "number of drt functions" from '
     #            //dsk,1,nbfdrt,0,' ')
      call iosys('read integer nrows from '//dsk,1,nrows,0,' ')
      call iosys('read integer nwks from '//dsk,1,nwks,0,' ')
      call iosys('read integer nlevs from '//dsk,1,nlevs,0,' ')
      call iosys('read integer orbfrm from '//dsk,1,orbfrm,0,' ')
      call iosys('read integer "symmetries in ci" from '//dsk,
     $           1,nsym,0,' ')
      call iosys('read integer norbs from '//dsk,1,norbs,0,' ')
      call iosys('read integer "number of basis functions"'//
     $           'from rwf',1,nbf,0,' ')
      nnp=norbs*(norbs+1)/2
      nnpbf=nbf*(nbf+1)/2
c
c     ----- check for other options -----
c
      nroots=intkey(ops,'ci=nroots',1,' ')
      call iosys('read integer roots from rwf',1,nroots,0,' ')
      nimprt=intkey(ops,'ci=nimportant',15,' ')
      srow=intkey(ops,'ci=starting_row',1,' ')
      cutoff=fpkey(ops,'ci=coefficient=cutoff',0.01d+00,' ')
      if(.not.iskey(ops,'ci=no=average',refops)) then
c        the average keyword does not exist.
         doavg=.false.
         avgall=.false.
      else if(chrkey(ops,'ci=no=average',' ',refops).eq.' ') then
c        the average keyword exists as a standalone.
         doavg=.true.
         avgall=.true.
      else
c        the average keyword exists as a replacement string.
         doavg=.true.
         avgall=.false.
      endif
      trans=logkey(ops,'properties=transition',.false.,' ')
c
c     ----- prepare for scattering option -----
c
      if( scat ) then
         trans=.true.
c
c        read the number of internal orbitals from input
c        or from the rwf file.
c 
         nsmall=intkey(ops,'scattering=nsmall',0,' ') 
         call iosys('does "internal orbitals" exist on rwf',0,0,0,ans)
         if(ans.eq.'yes') then
            call iosys('read integer "internal orbitals" from rwf',
     1                   1,nsmall,0,' ')
         endif
         if(nsmall.eq.0) then
            write(iout,*)' nsmall is needed for scattering density'
            call lnkerr(' m921: input error')
         end if
         npvec=intkey(ops,'scattering=npvec',1,' ')
         call iosys('open '//filtyp//' as new',0,0,0,nmfile)
      endif
c
c     retrieve the basis function labels.
      call iosys('read character "basis function labels" from rwf',
     $           -1,0,0,bflabl)
c
c     ----- the number of important configurations must be no larger
c              than the number of configurations!
c
      if(logkey(ops,'print=walks=full',.false.,' ')) then
        nimprt=nwks
      end if
c
      nroots=min(nroots,nwks)
      nimprt=min(nimprt,nwks)
      nlarge=min(10000,nwks)
      if (cutoff.lt.0.0001d+00) nlarge=nwks
c
c     ----- allocate core for drt arrays, vectors, etc. --
c
      walk=1
      brkdn=walk+nimprt
      arc=brkdn+norbs
      weight=arc+4*nrows
      orbtbf=weight+4*nrows
      ref=orbtbf+norbs
      alpha=ref+norbs
      beta=alpha+norbs
      bfsym=beta+norbs
      symnum=bfsym+nbfdrt
      nlwks=symnum+nbfdrt
      orbsym=nlwks+nrows
      b=orbsym+norbs
      irowsv=b+nrows
      segsv=irowsv+nlevs
      pagesv=segsv+nlevs
      iwtsv=pagesv+nlevs
      traksv=iwtsv+nlevs
      jrowsv=traksv+nlevs
      jwtsv=jrowsv+nlevs
      hrowsv=jwtsv+nlevs
      hwt=hrowsv+nlevs
      harcsv=hwt+nlevs
      bftorb=harcsv+nlevs
      lgwk=bftorb+nbf
c
      coef=iadtwp(lgwk+nlarge)
      c=coef+nimprt
      if (trans) then
         s=c+nwks
      else
         s=c
      end if
      lgcoef=s+nwks
      acoef=lgcoef+nlarge
      bcoef=acoef+nlevs
      dave=bcoef+nlevs
      if(doavg) then
         dscr=dave+nnp
         d1=dscr+nnp
      else
         d1=dave
      endif
      natocc=d1+nnp
      eigvec=natocc+norbs
      t1=eigvec+norbs**2
      t2=t1+nbf**2
      aotono=t2+nbf**2
      occ=aotono+nbf**2
      aotomo=occ+nbf
      motono=aotomo+nbf**2
      smat=motono+nbf**2
      smhalf=smat+nnpbf
      tocc=smhalf+nnpbf
      eigval=tocc+nbf
      u=eigval+nbf
      triang=u+nbf**2
      top=wpadti(triang+nnpbf)
c
c     ----- get core, and then read in arrays -----
c
      call getmem(top,p,ngot,'m921',0)  
c
      call iosys('read integer arc from '//dsk,4*nrows,a(arc),0,' ')
      call iosys('read integer weight from '//dsk,
     $            4*nrows,a(weight),0,' ')
      call iosys('read integer orbtbf from '//dsk,norbs,a(orbtbf),0,' ')
      call iosys('read integer iout from '//dsk,nbfdrt,a(bftorb),0,' ')
      call iosys('read integer bfsym from '//dsk,nbfdrt,a(bfsym),0,' ')
      call iosys('read integer nlwks from '//dsk,nrows,a(nlwks),0,' ')
      call iosys('read integer b from '//dsk,nrows,a(b),0,' ')
      call iosys('read integer orbsym from '//dsk,norbs,a(orbsym),0,' ')
      do 1 i=orbsym,orbsym+norbs-1
         a(i)=a(i)-1
    1 continue
c
c     read in the appropriate vector set.
      call iosys('read character "transformation vector" from rwf',
     $           -1,0,0,xform)
      if(ci) then
         call iosys('read real '//xform//' from rwf',nbf**2,
     $           z(aotomo),0,' ')
c        handle the case of an mcscf-ci calculation.
         if(xform(1:14).eq.'"mcscf vector"') then
            call iosys('read real "mcscf orbital energies" from rwf',
     $                  nbf,z(tocc),0,' ')
         else
            call iosys('read real "orbital energies" from rwf',
     $                  nbf,z(tocc),0,' ')
         end if
      else if (mcscf) then
         if(xform.eq.'"mcscf failed"') then
c
c            ---- read unconverged mcscf vectors
c
             call iosys('read real "mcscf unconverged vector" from rwf',
     $                   nbf**2,z(aotomo),0,' ')
             call iosys('read real "mcscf unconverged orbital energies"'
     $                  //'  from rwf',nbf,z(tocc),0,' ')
         else
c
c           ---- read   converged mcscf vectors
c
            call iosys('read real "mcscf vector" from rwf',nbf**2,
     $                  z(aotomo),0,' ')
            call iosys('read real "mcscf orbital energies" from rwf',
     $                  nbf,z(tocc),0,' ')
         end if
      else
         call iosys('read real "scf vector" from rwf',nbf**2,
     $               z(aotomo),0,' ')
         call iosys('read real "orbital energies" from rwf',
     $               nbf,z(tocc),0,' ')
      end if
      call rzero(z(dave),nnp)
c
c     ----- fix up bftorb in the case of an mcscf -----
c
      if (mcscf) then
         if(nbfdrt.lt.nbf) then
            call iosys('read integer mc_nfrozen from rwf',1,nfob,0,' ')
            call iosys('read integer mc_ncore from rwf',1,ncob,0,' ')
            call iosys('read integer mc_nactive from rwf',1,naob,0,' ')
c
            if (naob.ne.nbfdrt) call lnkerr('in ci analysis code '//
     #           'number of mcscf active orbitals does not match '//
     #           'the number of orbitals in the drt')
c
            call izero(a(bftorb),nbf)
            do 23 i=bftorb,bftorb+nfob+ncob-1
               a(i)=-1
   23       continue
            do 24 i=bftorb+nfob+ncob,bftorb+nfob+ncob+naob-1
               a(i)=i-bftorb-nfob-ncob+1
   24       continue
         end if
      end if
c
c     section to fix up orbtbf/bftorb in case of m903.
c     this is independent of whether it was run by mcscf or ci.
c     m903 assumes the original vector list, even though m811
c     may have written a transposed one.
c      if(m903) then
c         ncore=0
c         do 21 i=1,nbf
c            if(a(bftorb+i-1).eq.-1) then
c               ncore=ncore+1
c            else if (a(bftorb+i-1).ne.0) then
c               a(bftorb+i-1)=i
c            endif
c   21    continue
c         do 22 i=1,norbs
c            a(orbtbf+i-1)=ncore+i
c   22    continue
c      endif
c
c
c     ----- find the most important configurations in each root -----
c
      if(.not.scat) then
         do 1000 root=1,nroots
            write (iout,2) root
    2       format (/,'    most important configurations for root',i3,/)
c
c           reset the number of important configurations to be printed
            imprt=nimprt
            call iosys('read real "'//vn(1:ilen)//' root '
     #                             //itoc(root)//'"'
     #                            //' from rwf',nwks,z(c),0,' ')
            if(logkey(ops,'print=ci-vector',.false.,' ')) then
               write(iout,99012) root
               write(iout,99013)(z(c-1+i),i=1,nwks)
99012          format(//,' ci vector for root ',i6,/)
99013          format(5(2x,f8.5))
            end if
c
            if(logkey(ops,'print=walks=full',.false.,' ')) then
               do 99014 i=1,nwks
                  z(coef-1+i)=z(c-1+i)
                  a(walk-1+i)=i
99014          continue
            else
               call import(z(c),nwks,z(coef),a(walk),imprt,cutoff,
     #                     nlarge,z(lgcoef),a(lgwk))
            end if
c
            call prwalk(ops,z(coef),a(walk),a(brkdn),a(arc),a(weight),
     $           imprt,nrows,a(orbtbf),norbs,a(ref),a(alpha),a(beta),
     #           a(bfsym),a(symnum),nsym,nbfdrt,symlab,iout)
c
            call onepdm(a(arc),a(weight),a(nlwks),a(orbsym),a(b),srow,
     #                  z(c),z(c),z(t1),
     #                  nrows,norbs,nlevs,orbfrm,nsym,nwks,nnp,
     #                  a(irowsv),a(segsv),a(pagesv),a(iwtsv),a(traksv),
     #                  a(jrowsv),a(jwtsv),a(hrowsv),a(hwt),a(harcsv),
     #                  z(acoef),z(bcoef),trans)
c
            call fold(z(d1),z(t1),norbs,nnp)
c
c           add the density to the average, if requested.
            if(doavg) then
               if(avgall) then
c                 we are including all roots in the average.
                  call vadd(z(dave),z(dave),z(d1),nnp)
               else
c                 is this in the list of roots to include in the average?
                  if(logkey(ops,'ci=no=average='//itoc(root),.false.,
     $                      refops)) then
c                 it is included. extract associated weight (default=1.0).
                     wgt=fpkey(ops,'ci=no=average='//itoc(root),one,
     $                         refops)
                     call smul(z(dscr),z(d1),wgt,nnp)
                     call vadd(z(dave),z(dave),z(dscr),nnp)
                  else
c                    this root is not included in the average.
                  endif
               endif
            endif
c
c
c           reform the density matrix into the original transformation
c           vector sequence.
            call trtosq(z(t1),z(d1),norbs,nnp)
            call dordr(z(t2),z(t1),nbf,norbs,a(bftorb),.false.)
            call iosys('write real "mo 1pdm '//itoc(root)//'" to rwf',
     #                  nbf*nbf,z(t2),0,' ')
c..bhl
            if(logkey(ops,'print=ci=density',.false.,' ')) then
                write(iout,99011)
99011           format(/,' mo 1-particle density matrix',/)
                call matout(z(t2),nbf,nbf,nbf,nbf,iout)
            end if
c..bhl
c           transform it to the ao basis.
            call ebct(z(t1),z(t2),z(aotomo),nbf,nbf,nbf)
            call ebc(z(t2),z(aotomo),z(t1),nbf,nbf,nbf)
            call iosys('write real "ao 1pdm '//itoc(root)//'" to rwf',
     #                  nbf*nbf,z(t2),0,' ')
c
            call natorb(z(d1),nnp,z(natocc),z(eigvec),z(t1),z(t2),norbs,
     #                  iout,a(bftorb),z(aotono),z(occ),z(aotomo),
     #                  z(motono),nbf,ops,bflabl,mcscf,
     #                  z(smat),z(smhalf),z(tocc),z(u),z(eigval),nnpbf,
     #                  z(triang))
c
            call iosys('write real "no vector '//itoc(root)//'" to rwf',
     #                  nbf**2,z(aotono),0,' ')
            call iosys('write real "no occ '//itoc(root)//'" to rwf',
     #                  nbf,z(occ),0,' ')
c
 1000    continue
      end if
c
c     generate average natural orbitals, if requested.
      if(doavg) then
         write(iout,1002)
         call natorb(z(dave),nnp,z(natocc),z(eigvec),z(t1),z(t2),norbs,
     #               iout,a(bftorb),z(aotono),z(occ),z(aotomo),
     #               z(motono),nbf,ops,bflabl,mcscf,
     #               z(smat),z(smhalf),z(tocc),z(u),z(eigval),nnpbf,
     #               z(triang))
c
         call iosys('write real "no vector '//itoc(0)//'" to rwf',
     $               nbf*nbf,z(aotono),0,' ')
c
         call iosys('write real "no occ '//itoc(0)//'" to rwf',
     $               nbf,z(occ),0,' ')
      endif
c
c     ----- set up for transition density matrix if requested -----
c
c...bhl modified for optical potential 2/28/89
c
      if ( scat ) then
         call iosys ('write integer "p-space nwks" to '//filtyp,
     $                1,nwks,0,' ')
c        do full lower triangle of transition density matrices.
         do 200 root1=1,nroots
            call iosys('read real "ci root '//itoc(root1)//
     #                 '" from rwf',nwks,z(c),0,' ')
            call iosys('write real "p-space ci root '//itoc(root1)//
     #                 '" to '//filtyp,nwks,z(c),0,' ')
            if(logkey(ops,'print=ci-vector',.false.,' ')) then
               write(iout,99012) root1
               write(iout,99013)(z(c-1+i),i=1,nwks)
            end if
            do 190 root2=1,root1
c
               call iosys('read real "ci root '//itoc(root2)//
     #                    '" from rwf',nwks,z(s),0,' ')
c
               write(iout,99015) root1, root2
               call onepdm(a(arc),a(weight),a(nlwks),a(orbsym),
     #                     a(b),srow,z(c),z(s),z(t1),
     #                     nrows,norbs,nlevs,orbfrm,nsym,nwks,nnp,
     #                     a(irowsv),a(segsv),a(pagesv),a(iwtsv),
     #                     a(traksv),a(jrowsv),a(jwtsv),a(hrowsv),
     #                     a(hwt),a(harcsv),z(acoef),z(bcoef),trans)
c
c
c              reform the density matrix into the original transformation
c              vector sequence.
c
               call dordr(z(t2),z(t1),nbf,norbs,a(bftorb),trans)
c
               call iosys('write real "mo t1pdm:'//itoc(root1)//
     #                     itoc(root2)//'" to rwf',nbf**2,z(t2),
     $                     0,' ')
99015          format(/,'mo1pdm for root 1 = ',i3,
     1                /,'           root 2 = ',i3)     
c
               if(logkey(ops,'print=ci=tdensity',.false.,' ')) then
                  write(iout,1001) root1,root2
                  call wmat(z(t2),nbf,nbf,bflabl,bflabl)
               endif
c
               call kzero(z(t2),nsmall,nbf)
               if(logkey(ops,'print=scattering=tdensity',
     #                                         .false.,' '))then
                  write(iout,1001) root1,root2
                  call wmat(z(t2),nbf,nbf,bflabl,bflabl)
               endif
               call iosys('write real "mo t1pdm:'//itoc(root1)//
     #                     itoc(root2)//'" to '//filtyp,nbf**2,z(t2),
     $                     0,' ')
c
c              transform it to the ao basis.
               call ebct(z(t1),z(t2),z(aotomo),nbf,nbf,nbf)
               call ebc(z(t2),z(aotomo),z(t1),nbf,nbf,nbf)
c
               call iosys('write real "ao t1pdm:'//itoc(root1)//
     #                     itoc(root2)//'" to '//filtyp,nbf**2,z(t2),
     #                     0,' ')
c
  190       continue
  200    continue
      else if (trans) then
c        do only requested transition densities.
        do 300 root1=1,nroots
            if(logkey(ops,'properties=transition='//itoc(root1),
     $         .false.,' ')) then
               call iosys('read real "ci root '//itoc(root1)//
     #                     '" from rwf',nwks,z(c),0,' ')
c              note that we do the upper triangular entries requested.
               do 290 root2=root1+1,nroots
c 
                  if(logkey(ops,'properties=transition='//itoc(root2),
     $                      .false.,' ')
     $           .or.logkey(ops,'properties=transition=*',.false.,' '))
     $            then
                     call iosys('read real "ci root '//itoc(root2)//
     #                          '" from rwf',nwks,z(s),0,' ')
c
                     call onepdm(a(arc),a(weight),a(nlwks),a(orbsym),
     #                     a(b),srow,z(c),z(s),z(t1),
     #                     nrows,norbs,nlevs,orbfrm,nsym,nwks,nnp,
     #                     a(irowsv),a(segsv),a(pagesv),a(iwtsv),
     #                     a(traksv),a(jrowsv),a(jwtsv),a(hrowsv),
     #                     a(hwt),a(harcsv),z(acoef),z(bcoef),trans)
c 
c                     reform the density matrix into the original transformation
c                     vector sequence.
c
                     call dordr(z(t2),z(t1),nbf,norbs,a(bftorb),trans)
c 
                     call iosys('write real "mo t1pdm:'//itoc(root1)//
     #                        itoc(root2)//'" to rwf',nbf**2,z(t2),
     $                        0,' ')
c 
                     if(logkey(ops,'print=ci=tdensity',.false.,' '))then
                        write(iout,1001) root1,root2
                        call wmat(z(t2),nbf,nbf,bflabl,bflabl)
                     endif
c
c                    transform it to the ao basis.
                     call ebct(z(t1),z(t2),z(aotomo),nbf,nbf,nbf)
                     call ebc(z(t2),z(aotomo),z(t1),nbf,nbf,nbf)
c
                     call iosys('write real "ao t1pdm:'//itoc(root1)//
     #                       itoc(root2)//'" to rwf',nbf**2,z(t2),
     $                       0,' ')
c
                  endif
  290          continue
            endif
  300    continue
      endif
c
c     ----- and exit with grace -----
c
 5000 continue
c
c     ----- close the files
c   
      if(scat) then
         call iosys('close '//filtyp,0,0,0,' ')
      endif

      call getmem(-ngot,p,idum,'m921',idum)
      call chainx(0)
c
c
      stop
      end
