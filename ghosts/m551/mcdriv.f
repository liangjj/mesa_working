*deck %W%  %G%
      subroutine mcdriv(cr,icr,nwint,prtflg)
c
c***begin prologue     mcdriv
c***date written       871022   (yymmdd)
c
c***revision date      910605   (yymmdd)
c
c    5 june     1991  rlm at lanl
c      adding option to converge on excited mcscf root.
c   19 april    1990  rlm at lanl
c      working on 32 bit integer version.
c   28 november 1989  rlm at lanl
c      modifying to print basis function labels as row indices
c      for vector printing.
c   16 december 1987
c      ipfdm set to 1 on last iteration
c   16 december 1987
c      mo mcscf 1pdm and 2pdm written to moden
c      if opnscf=1
c   26 october 1987    pws at lanl
c      changing the input section to avoid the namelist and
c      to supply as many of the variables as possible from
c      previous input.
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
c***end prologue       mcdriv
c
      implicit real*8 (a-h,o-z)
      parameter (maxnbf=2000)
      parameter (maxsym=20)
c
c
      character*3 answer
      character*8 chrkey
      character*8 citype
      character*8 dcaptl
      character*4096 mcinp
      character*4096 ops
      character*80 cjunk
      character*8 prtflg,mcscf
      character*16 bflabl(maxnbf)
      logical dollar
      logical logkey
      logical bypass
      logical abort
      logical otout
c
      integer opnscf,wpadti,intkey
      integer icas(maxsym), nob(maxsym), nbf(maxsym)
      integer nfob(maxsym), ncob(maxsym),naob(maxsym)
      integer nato(maxsym), nvob(maxsym),nocc(maxsym)
      integer nobt(maxsym), isymm(maxsym+1),locsym(maxsym+1)
      integer nmixs(maxsym)
      integer iobsym(500),iobcas(100)
      integer mcmthd(50)
c
      real*8 fpkey,temp,energy,sqcdf
      real*8 smlld(25),shftd(25)
c
      common /number/ zero, pt5, one, two, four, eight
      common /io/ inp,iout
      real*8 cr(1)
      integer icr(2)
c
      dimension kksym(256)
c
c
c
c
c
      data nbf/maxsym*0/, nfob/maxsym*0/, ncob/maxsym*0/
      data naob/maxsym*0/, nob/maxsym*0/
      data nato/maxsym*1/
      data nstate/1/, istate/1/, itcorh/1/
      data iocore/1/, itflag/1/, itran/1/, thrsh/1.d-12/
      data iobsym/500*0/, iobcas/100*0/, icas/maxsym*0/
      data throv/1.d-12/,thmcdv/10.d0/, nporb/0/
      data mcmthd/50*1/, nmcit/10/, thrmc/1.d-9/, nvcset/0/
      data nleqps/1/, nleqit/15/, thrleq/1.d-4/, nheig/1/
      data thre/1.d-15/, nviter/3/, innorm/0/,iprtci/0/
      data lbuf84/20000/,nhd/750/,ldar/3840/
      data ipfdm/0/,igeom/0/,opnscf/0/,iprtg/0/
      data smlld/25*-25.d0/,shftd/25*0.0d0/,icutah/40/,nmaxah/20/,
     $     thrcah/1.d-6/,iauto/1/,osqcdf/1.d0/,autoso/1.d-3/,mcroot/0/
      data icphf/1/,lbcphf/3840/,iorder/-1/,jktrn/0/,lsect/512/
      data lbufso /20000/,iprtmo/0/,ksym/0/,ifock/1/,iprtlag/0/
      data lendia/20000/, nphes/0/, lbufh/20000/
c
      external dcaptl
c
c      namelist /mcinp/ name, nsym, nbf, nob, nfob, ncob, naob,
c     $     iobsym, iobcas, icas, nstate, istate, nvcset,
c     $     ind, throv, iocore, itflag, thrsh, itran, mcmthd, nphes,
c     $     nleqps, nleqit, thrleq, jnd, knd, lnd, mnd, itcorh,
c     $     nheig, thre, nviter, innorm, iprtg,smlld,shftd,
c     $     lbufh, lbuf84, nhd, ldar, ldsk60, lbf11, iorder,
c     $     nmcit, thrmc, nvcset, nporb, ipfdm, igeom, opnscf, nato,
c     $     icutah, nmaxah, thrcah, iauto,autoso, icphf,lbcphf,jktrn,
c     $     lsect,lbufso,iprtmo,ksym,kksym,ifock,iprtlag,iprtci,mcroot
c
c
c     load common/number/.
      zero=0.0d+00
      pt5=0.5d+00
      one=1.0d+00
      two=2.0d+00
      four=4.0d+00
      eight=8.0d+00
c
c
c-------------------------------
c     read input
c-------------------------------
c
      call iosys('read character options from rwf',-1,0,0,ops)
      call iosys('read character "basis function labels" from rwf',
     $            -1,0,0,bflabl)
c
c     ----- first set options from options string -----
c
      if (logkey(ops,'mcscf=cas',.false.,' ')) icas(1)=1
      call intarr(ops,'mcscf=cas',icas,maxsym,' ')
      nmcit=intkey(ops,'mcscf=cycles',10,' ')
      mcroot=intkey(ops,'mcscf=root',1,' ')
      call intarr(ops,'mcscf=method',mcmthd,50,' ')
      call fparr(ops,'mcscf=small-diagonals',smlld,25,' ')
      nleqit=intkey(ops,'mcscf=linear-equations=cycles',30,' ')
      thrleq=fpkey(ops,'mcscf=linear-equations=convergence',
     $             1.0d-04,' ')
      if (logkey(ops,'print=mcscf=density-matrix-energy',.false.,' '))
     $           ipfdm=1
      if (logkey(ops,'mcscf=out-of-core-transformation',.false.,' '))
     $           jktrn=1
      if (logkey(ops,'mcscf-out-of-core-hessian',.false.,' ')) itcorh=0
      if (logkey(ops,'mcscf=use-total-fock-operator',.false.,' '))
     $           ifock=0
      call fparr(ops,'mcscf=level-shift',shftd,25,' ')
      thrmc=fpkey(ops,'mcscf=convergence',1.0d-9,' ')
      if (logkey(ops,'hf=quadratic',.false.,' ')) iauto=0
      if (logkey(ops,'mcscf=nonautomatic',.false.,' ')) iauto=0
      if (logkey(ops,'print=mcscf=hessian',.false.,' ')) nphes=1
      if (logkey(ops,'print=mcscf=gradient',.false.,' ')) iprtg=1
      if (logkey(ops,'print=mcscf=vector',.false.,' ')) iprtmo=1
      if (logkey(ops,'print=mcscf=lagrangian',.false.,' ')) iprtlag=1
      if (logkey(ops,'print=mcscf=ci',.false.,' ')) iprtci=1
      icutah=intkey(ops,'mcscf=augmented-hessian=cut',40,' ')
      abort=.not.logkey(ops,'mcscf=noabort',.false.,' ')
      call iosys('does "symmetries in ci" exist on rwf',0,0,0,answer)
      if (answer.eq.'no') then
         junk=1
      else
         call iosys('read integer "symmetries in ci" from rwf',
     $               1,junk,0,' ')
      end if
      nsym=intkey(ops,'mcscf=number-of-symmetries',junk,' ')
      call iosys('read integer "number of basis functions" from rwf',
     $            1,nob,0,' ')
      call intarr(ops,'mcscf=number-of-orbitals',nob,nsym,' ')
      call icopy(nob,nbf,nsym)
      call intarr(ops,'mcscf=number-of-basis-functions',nbf,nsym,' ')
      call intarr(ops,'mcscf=number-of-frozen-orbitals',nfob,nsym,' ')
      call intarr(ops,'mcscf=number-of-core-orbitals',ncob,nsym,' ')
      call intarr(ops,'mcscf=natural-orbitals',nato,nsym,' ')
c
      call iosys('does "number of drt functions" exist on rwf',0,0,0,
     $            answer)
      if (answer.eq.'no') then
         naob(1)=1
      else
         call iosys('read integer "number of drt functions" from rwf',
     $               1,naob,0,' ')
      end if
      call intarr(ops,'mcscf=number-of-active-orbitals',naob,nsym,' ')
      call intarr(ops,'mcscf=invariant-orbital-sets',iobcas,100,' ')
      call intarr(ops,'mcscf=orbital-symmetry-label',iobsym,500,' ')
      bypass=logkey(ops,
     $              'mcscf=bypass-final-orbital-rotation',.false.,' ')
      if (logkey(ops,'hf=quadratic',.false.,' ')) opnscf=1
      if (logkey(ops,'mcscf=hf',.false.,' ')) opnscf=1
c
c     ----- read from a $mcscf section if it is there -----
      if (dollar('$mcscf',mcinp,cjunk,inp)) then
         call locase(mcinp,mcinp)
c
         if (logkey(mcinp,'mcscf=cas',.false.,' ')) icas(1)=1
         call intarr(mcinp,'mcscf=cas',icas,maxsym,' ')
         nmcit=intkey(mcinp,'mcscf=cycles',10,' ')
         mcroot=intkey(mcinp,'mcscf=root',1,' ')
         call intarr(mcinp,'mcscf=method',mcmthd,50,' ')
         call fparr(mcinp,'mcscf=diagonals=small',smlld,25,' ')
         nleqit=intkey(mcinp,'mcscf=linear-equations=cycles',30,' ')
         thrleq=fpkey(mcinp,'mcscf=linear-equations=convergence',
     $                1.0d-04,' ')
         if (logkey(mcinp,'print=mcscf=density-matrix-energy',
     $             .false.,' ')) ipfdm=1
         if (logkey(mcinp,'mcscf=out-of-core-transformation',
     $             .false.,' ')) jktrn=1
         if (logkey(mcinp,'mcscf-out-of-core-hessian',.false.,' '))
     $              itcorh=0
         call fparr(mcinp,'mcscf=level-shift',shftd,25,' ')
         thrmc=fpkey(mcinp,'mcscf=convergence',1.0d-9,' ')
         if (logkey(mcinp,'hf=quadratic',.false.,' ')) iauto=0
         if (logkey(mcinp,'mcscf=nonautomatic',.false.,' ')) iauto=0
         if (logkey(mcinp,'print=mcscf=hessian',.false.,' ')) nphes=1
         if (logkey(mcinp,'print=mcscf=gradient',.false.,' ')) iprtg=1
         if (logkey(mcinp,'print=m551=vector',.false.,' ')) iprtmo=1
         if (logkey(mcinp,'print=mcscf=lagrangian',.false.,' '))
     $              iprtlag=1
         if (logkey(mcinp,'print=mcscf=ci',.false.,' ')) iprtci=1
         icutah=intkey(mcinp,'mcscf=augmented-hessian=cut',40,' ')
         nsym=intkey(mcinp,'mcscf=number-of-symmetries',junk,' ')
         call intarr(mcinp,'mcscf=number-of-orbitals',nob,nsym,' ')
         call icopy(nob,nbf,nsym)
         call intarr(mcinp,'mcscf=number-of-basis-functions',
     $               nbf,nsym,' ')
         call intarr(mcinp,'mcscf=number-of-frozen-orbitals',
     $               nfob,nsym,' ')
         call intarr(mcinp,'mcscf=number-of-core-orbitals',
     $               ncob,nsym,' ')
         call intarr(mcinp,'mcscf=number-of-active-orbitals',
     $               naob,nsym,' ')
         call intarr(mcinp,'mcscf=invariant-orbital-sets',
     $               iobcas,100,' ')
         if (logkey(mcinp,'hf=quadratic',.false.,' ')) opnscf=1
      end if
c
c     ----- hf case fixup -----
c
      if (opnscf.eq.1) then
         thrmc=fpkey(ops,'scf=convergence',1.0d-5,' ')
         thrmc=thrmc**2
         nmcit=intkey(ops,'scf=cycles',10,' ')
         call iosys('read integer "number of alpha electrons" from rwf',
     $               1,nae,0,' ')
         call iosys('read integer "number of beta electrons" from rwf',
     $               1,nbe,0,' ')
c
         naob(1)=nae-nbe
         ncob(1)=nbe
         icas(1)=1
      end if
c
      if(ncob(1).eq.0) ifock=0
c
c-------------------------------
c     set up index arrays
c-------------------------------
c
      im = 0
      ls = 0
      do 800 i = 1, nsym
         nocc(i) = ncob(i) + naob(i)
         nvob(i) = nob(i) - nfob(i) - nocc(i)
         nobt(i) = nob(i) - nfob(i)
         isymm(i) = im
         im = im + nbf(i) * nocc(i)
         locsym(i) = ls
         ls = ls + nocc(i)
 800  continue
      isymm(nsym+1) = im
      locsym(nsym+1) = ls
c
c     ----- store numbers of orbitals in each type on rwf -----
c
      call iosys('write integer mc_nfrozen to rwf',nsym,nfob,0,' ')
      call iosys('write integer mc_ncore to rwf',nsym,ncob,0,' ')
      call iosys('write integer mc_nactive to rwf',nsym,naob,0,' ')
c
c     ----- put these critical arrays on the scratch unit -----
c
      junk=(nocc(1)**2*nbf(1)**2+nocc(1)*nbf(1)**3)/10
      junk=max(min(junk,2000000),256000)
      call iosys('open mcscr as scratch on ssd',junk,0,0,' ')
c
      call iosys('write character mcscf_type to mcscr',0,0,0,mcscf)
      call iosys('write integer nsym to mcscr',1,nsym,0,' ')
      call iosys('write integer nbf to mcscr',nsym,nbf,0,' ')
      call iosys('write integer nob to mcscr',nsym,nob,0,' ')
      call iosys('write integer nfob to mcscr',nsym,nfob,0,' ')
      call iosys('write integer ncob to mcscr',nsym,ncob,0,' ')
      call iosys('write integer naob to mcscr',nsym,naob,0,' ')
      call iosys('write integer nocc to mcscr',nsym,nocc,0,' ')
      call iosys('write integer nvob to mcscr',nsym,nvob,0,' ')
      call iosys('write integer nobt to mcscr',nsym,nobt,0,' ')
      call iosys('write integer icas to mcscr',nsym,icas,0,' ')
      call iosys('write integer method to mcscr',nmcit,mcmthd,0,' ')
      call iosys('write integer maximum_iterations to mcscr',1,
     $            nmcit,0,' ')
      call iosys('write real convergence_threshold to mcscr',
     $           1,thrmc,0,' ')
      call iosys('write real divergence_threshold to mcscr',1,
     $           thmcdv,0,' ')
      call iosys('write integer max_linear_eq_iterations to mcscr',
     $           1,nleqit,0,' ')
      call iosys('write real linear_eq_convergence to mcscr',1,
     $           thrleq,0,' ')
      call iosys('write integer number_aug_hess_roots to mcscr',
     $           1,nheig,0,' ')
      call iosys('write real aug_hess_convergence to mcscr',
     $           1,thre,0,' ')
      call iosys('write integer number_inverse_iterations to mcscr',
     $           1,nviter,0,' ')
c
c
      call iosys('read character mcscf_type from rwf',-1,0,0,mcscf)
      if(opnscf.ne.0) then
         mcscf='opnscf'
      endif
      call iosys('write character mcscf_type to rwf',0,0,0,mcscf)
      call iosys('write character mcscf_type to mcscr',0,0,0,mcscf)
c
c
       if(opnscf.eq.0) then
          call iosys('does "h diagonals" exist on rwf',0,0,0,answer)
          if(answer.eq.'no') then
             call iosys('read integer nwks from rwf',1,nwks,0,' ')
             maxdia=intkey(ops,'mcscf=maxdia',lendia,' ')
             maxdia=max(maxdia,nwks)
             call iosys('create real "h diagonals" on rwf',maxdia,0,0,
     $                  ' ')
             call iosys('create real "mc root 1" on rwf',maxdia,0,0,' ')
             call iosys('create real "mc root 2" on rwf',maxdia,0,0,' ')
          end if
       end if
c
c     ----- decide which ci to use -----
c
      if (mcscf.eq.'casscf') then
c
c        ----- for small casscf calculations, still use m902
c
         call iosys('read integer nwks from rwf',1,nwks,0,' ')
         minwks=intkey(ops,'mcscf=size902',-2,' ')
         if (nwks.le.minwks) then
            citype='m902'
         else
            citype='m903'
         end if
      else if (mcscf.eq.'mcscf') then
         citype='m902'
      else
         citype='hf'
      end if
c
c     ----- and check for an override -----
c
      cjunk=citype
      citype=chrkey(ops,'mcscf=ci-to-use',cjunk,' ')
      cjunk=citype
      citype=dcaptl(cjunk)
      citype=cjunk
c
      call iosys('write character "mcscf: ci used" to rwf',
     $            0,0,0,citype)
c
      if (citype.ne.'hf') then
         write (iout,3678) citype,mcscf,nwks
 3678    format (5x,'using ',a4,' for ',a6,' of ',i5,' configurations')
      else
         write (iout,3679)
 3679    format(5x,'mcscf used for quadratic scf')
      end if
c
c----------------------------
c     print input data
c----------------------------
c
      if(iprtci.ne.0) then
         write(iout,17011)
17011    format(/,' ****** the print flag has been turned on ******',/)
         prtflg=' '
         call iosys('write character "print flag" to rwf',0,0,0,prtflg)
      endif
      if (prtflg.ne.'minimum') then
         write (iout,7011)
 7011    format('0******** basis set and orbital definitions ')
c
         i2 = 0
         do 300 isym = 1, nsym, 15
            i1 = i2 + 1
            i2 = min(i2+15,nsym)
            write (iout,7012) (i,i=i1,i2)
            write (iout,7013) (nbf(i),i=i1,i2)
            write (iout,7014) (nob(i),i=i1,i2)
            write (iout,7015) (nfob(i),i=i1,i2)
            write (iout,7016) (ncob(i),i=i1,i2)
            write (iout,7017) (naob(i),i=i1,i2)
            write (iout,7018) (icas(i),i=i1,i2)
 300     continue
 7012    format(' sym   ',15(1x,i3))
 7013    format(' nbf   ',15(1x,i3))
 7014    format(' nob   ',15(1x,i3))
 7015    format(' nfob  ',15(1x,i3))
 7016    format(' ncob  ',15(1x,i3))
 7017    format(' naob  ',15(1x,i3))
 7018    format(' icas  ',15(1x,i3))
c
         write (iout,7020)
 7020    format(' *** orbital sub-symmetry labels(iobsym)')
         l0 = 0
         do 320 isym = 1, nsym
            write (iout,7021) isym
 7021       format(' symmetry       ',1x,i3)
            n = nob(isym)
            write (iout,7022) (i,i=1,n)
 7022       format(' orbital        ',15(1x,i3)/(16x,15(1x,i3)))
            write (iout,7023) (iobsym(l0+i),i=1,n)
 7023       format(' iobsym         ',15(1x,i3)/(16x,15(1x,i3)))
            l0 = l0 + n
 320     continue
c
         write (iout,7025)
 7025    format(' *** orbital free-rotation group labels(iobcas)')
         l0 = 0
         do 350 isym = 1, nsym
            write (iout,7026) isym
 7026       format(' symmetry       ',1x,i3)
            n = naob(isym)
            write (iout,7027) (i,i=1,n)
 7027       format(' active orbital ',15(1x,i3)/(16x,15(1x,i3)))
            write (iout,7028) (iobcas(l0+i),i=1,n)
 7028       format(' iobcas         ',15(1x,i3)/(16x,15(1x,i3)))
            l0 = l0 + n
 350     continue
c
         write (iout,7030)
 7030    format(' ******** mcscf run control and thresholds')
         write (iout,7031)
 7031    format(' *** mcscf method(mcmthd): 1, augmented hessian',
     $          ' ; 2, second order with ci coupling')
         write (iout,7032) (i,i=1,nmcit)
         write (iout,7033) (mcmthd(i),i=1,nmcit)
 7032    format(' iteration ',30(1x,i1))
 7033    format(' mcmthd    ',30(1x,i1))
c
         write (iout,7035)
 7035    format(' *** run controls and thresholds')
         write (iout,7036) nmcit
 7036    format(' ncmit  ',4x,i5,3x,
     $          ' maximum number of mcscf iterations ')
         write (iout,7037) thrmc
 7037    format(' thrmc  ',1d9.2,3x,
     $          ' mcscf convergence threshold        ')
         write (iout,7038) thmcdv
 7038    format(' thmcdv ',1d9.2,3x,
     $          ' mcscf divergence threshold         ')
         write (iout,7040) nleqit
 7040    format(' nleqit ',4x,i5,3x,
     $          ' no of linear equation iterations   ')
         write (iout,7041) thrleq
 7041    format(' thrleq ',1d9.2,3x,
     $          ' lin. eq. convergence threshold     ')
         write (iout,7045) nheig
 7045    format(' nheig  ',4x,i5,3x,
     $          ' no of augmented hessian roots      ')
         write (iout,7046) thre
 7046    format(' thre   ',1d9.2,3x,
     $          ' augmented hessian conv. threshold  ')
         write (iout,7047) nviter
 7047    format(' nviter ',4x,i5,3x,
     $          ' no of aug. hes. inverse iterations ')
c
         write (iout,7050)
 7050    format(' *** hessian construction controls')
         write (iout,7051) iocore
 7051    format(' iocore ',4x,i5,3x)
         write (iout,7052) itflag
 7052    format(' itflag ',4x,i5,3x)
         write (iout,7053) itran
 7053    format(' itran  ',4x,i5,3x)
         write (iout,7054) thrsh
 7054    format(' thrsh  ',1d9.2,3x)
c
         write (iout,7060)
 7060    format(' ******** dataset information')
         write (iout,7061) lbufh
 7061    format(' lbufh  ',4x,i5,3x,
     $          'orbital hessian record length, nf   ')
         write (iout,7062) lbuf84
 7062    format(' lbuf84 ',4x,i5,3x,
     $          'orbital hessian record length, nf   ')
         write (iout,7063) nhd
 7063    format(' nhd    ',4x,i5,3x,
     $          'diag. element record length on nf81 ')
         write (iout,7064) ldar
 7064    format(' ldar   ',4x,i5,3x,
     $          'record length of da files 16 and 60 ')
      endif
c
c     allocate some core.
c
      imix=1
c
ccccc      if(icphf.ne.0) imix=imix+lbcphf
c
      imixh = imix + im
      imixi = imixh + im
      imixo = imixi + im
      imixl = imixo + im
      imixu = imixl + im
      ilen  = imixu + im
      iloc  = ilen + ls
      ilocp = iloc + ls
      need  = ilocp+ls
crlm     lstor = iadtwp(ilocp + ls)
crlm      need=lstor+nob(1)**2
crlm      need=wpadti(need)
      call getscm(need,icr,ngot,'mcmxdr',0)
crlm     nwwp = iadtwp(ngot) - lstor
c
c-------------------------------------------
c     generate non-redundant orbital mixings
c-------------------------------------------
c
      call mcmxdr(icr(imix),icr(imixo),
     $            icr(imixh),icr(imixi),
     $            icr(imixl), icr(imixu), icr(ilen),
     $            icr(iloc), icr(ilocp),
     $            nocc, nobt, ncob, naob, nvob,iobsym,
     $            iobcas, icas, nmix, nsym, isymm, locsym, nmixs,nbf,
     $            prtflg)
c
c--------------------------
c     order integrals
c--------------------------
c
c     get some more core
      lstor=iadtwp(need)
      need=wpadti(lstor+nob(1)**2)
      call getscm(need,icr,ngot,'mcdriv',0)
      nwwp=iadtwp(ngot)-lstor
      if(iorder.ge.0) then
         write(iout,*) '   ordering integrals '
cos      call mn330(cr(lstor),cr(lstor),nwwp)
         if(iorder.gt.0) call lnkerr('iorder gt 0')
      endif
c
c-------------------------------
c     read input orbitals
c-------------------------------
c
      ntot=nob(1)*nob(1)
      call iosys('read real "guess vector" from rwf',
     $            ntot,cr(lstor),0,' ')
      call iosys('create real mo_new on mcscr',ntot,0,0,' ')
      call iosys('create real mo_old on mcscr',ntot,0,0,' ')
      call iosys('write real mo_new on mcscr',ntot,cr(lstor),
     $            0,' ')
c
      otout = .true.
      iter = 0
      ilast=0
      write (iout,4100)
 4100 format (t8,'iteration',t20,'mcscf energy',t40,'convergence')
c
c----------------------------------
c----------------------------------
c     mcscf iteration loop
c----------------------------------
c----------------------------------
c
 1000 continue
      iter = iter + 1
      iaugh = 1
      icicup = 0
c
      if(iauto.ne.0.and.osqcdf.lt.autoso) mcmthd(iter)=2
      if(ilast.eq.1.and.icphf.ne.0) mcmthd(iter)=2
      if (mcmthd(iter) .ne. 1) then 
         iaugh = 0
         icicup = 1
      end if
c
      if(iter.le.25) then
         smlldd=smlld(iter)
         shftdd=shftd(iter)
      end if
c
c
      call mciter(nsym, nbf, nfob, ncob, naob, nvob, nmix,
     $     icr(imix), icr(imixh), icr(imixi), icr(imixo),
     $     icr(ilen), icr(iloc),  icr(ilocp), icas,
     $     nstate, istate,
     $     iocore, itflag, thrsh, itran,
     $     icicup, nleqps, nleqit, thrleq,
     $     iaugh, nheig, thre, nviter, innorm,
     $     nf14, nf16, nf30, nf36, nf39,
     $     nf46, nf49, nf81, nhd, nf82, nf83, nf84, lbuf84,
     $     nf91, nf92, nf93, lbufh, nf94,
     $     sqcdf, otout, ipfdm, igeom, opnscf, iprtg,smlldd,shftdd,
     $     cr(lstor),cr(lstor), wpadti(nwwp), nphes, ldar,
     $     icutah,nmaxah,thrcah, osqcdf,icphf,iter,jktrn,lsect,
     $     lbufso,ilast,itcorh,mcroot,ops)
c
      call iosys('read real "mcscf energy" from rwf',1,energy,0,' ')
      if(ilast.ne.1) write (iout,9102) iter,energy,sqcdf
 9102 format (t10,i5,t20,f15.9,t40,1pe10.3)
      endfile iout
      backspace iout
c
      otout = .false.
c
c     convergence test.
      if (sqcdf .le. thrmc)  then
         if(ilast.eq.0) then
            ilast=1
            ipfdm=1
            write(iout,1903)
 1903       format(/,' ** mcscf converged -- computing final energy ',
     $           'and hessian **')
            num=nob(1)
            nnp=num*(num+1)/2
c
c           allocate core, again.
            is=lstor
            ism=is+nnp
            iu=ism+nnp
            it1=iu+num**2
            it2=it1+num**2
            ieig=it2+num**2
            itrng=ieig+num
            icmo=itrng+nnp
            llstor=icmo+num**2
            need=wpadti(llstor)
            call getscm(need,icr,ngot,' ',0)
            nwwp=iadtwp(ngot)-llstor
c
c     ----- read in overlap integrals -----
c
            call iosys('read real "overlap integrals" from rwf',
     $                  nnp,cr(is),0,' ')
c
c     ----- form s**(1/2) and s**(-1/2) -----
c
            call sinv(cr(is),cr(ism),cr(iu),cr(ieig),cr(it1),cr(it2),
     $                num,nnp,cr(itrng),iprint)
c
            call iosys('read real mo_old from mcscr',
     $                  ntot,cr(icmo),0,' ')
c
            call pschmd(cr(icmo),cr(is),cr(iu),cr(it1),cr(it2),
     $                  num,num,nnp,1.0d-04)
c
            if(.not.bypass) then
               call natodr(cr(llstor),cr(icmo),nwwp,nfob,ncob,
     $                     naob,nvob,nob,nbf,nsym,nato,cr(ieig),
     $                     ifock,opnscf,iobcas)
            end if
c
            call pschmd(cr(icmo),cr(is),cr(iu),cr(it1),cr(it2),
     $                  num,num,nnp,1.0d-04)
c
c     ----- enforce symmetry by rotating degenrate orbitals -----
c
cbl
c      call enforc(num,nnp,cr(is),cr(ism),cr(ieig),cr(icmo),
c     #            cr(it1),cr(it2))
cbl
            if (opnscf.eq.1) then
               call iosys('write real "orbital energies" to rwf',
     $              nbf(1),cr(ieig),0,' ')
            else
               call iosys('write real "mcscf orbital energies" to rwf',
     $              nbf(1),cr(ieig),0,' ')
            end if
c
            if(iprtmo.ne.0) then
               write(iout,88188)
88188          format(//,' *** final mcscf vectors ***',/)
               call wvec(cr(icmo),cr(ieig),num,num,bflabl,' ')
            endif
            call iosys('write real mo_new on mcscr',
     $                  ntot,cr(icmo),0,' ')
c
         else
c           we have converged.
            go to 2000
         end if
c
      end if
c
      if (nporb.ne.0) then
         write (iout,9103)
         icmo = lstor
         call iosys('read real mo_new from mcscr',ntot,cr(icmo),0,' ')
c
         write(iout,9103)
         call wmat(cr(icmo),nbf(1),nob(1),'bflabl',' ')
 9103    format(1x,'output orbitals')
      end if
c
c------------------------------
c     divergence 
c     test and return if diverging. 
c------------------------------
      if (sqcdf.ge.thmcdv) then
         write (iout,9400)
 9400    format(1x,'****** mcscf diverged '/
     $          ' improved trial vectors or damping are recommended')
         return
      end if
      if (iter .lt. nmcit) go to 1000
c
c-----------------------------------------------
c-----------------------------------------------
c     end of iteration loop
c-----------------------------------------------
c-----------------------------------------------
c
c
c------------------------------
c     convergence failure
c     falls through to here if not converged; cleans up and returns
c------------------------------
      icmo=lstor
c
      call iosys('read real mo_old from mcscr',ntot,cr(icmo),0,' ')
      cjunk='"mcscf failed"'
      call iosys('write character "transformation vector" to rwf',
     $            0,0,0,cjunk(1:16))
      call iosys('write real "mcscf unconverged vector" to rwf',
     $            nbf(1)**2,cr(icmo),0,' ')
      call iosys('write real "mcscf unconverged orbital energies" '//
     $           'to rwf',nbf(1),cr(icmo),0,' ')
c
      write (iout,9300) nmcit
 9300 format(/,5x,'mcscf failed to converge in ',i5,' iterations')
c
      if(abort) then
         call lnkerr('mcscf failed to converge')
      else
         write (iout,9301)
 9301    format(/,5x,'mcscf vectors will be saved ')
      end if
      return 
c
c----- comes here if converged -----
c------------------------------------------------
c     convergence
c------------------------------------------------
c
 2000 continue
      call iosys('read real mo_old from mcscr',
     $            ntot,cr(icmo),0,' ')
c
      if(ksym.ne.0) then
         write(iout,20001)
20001    format(/,'  reorder the orbitals for symmetry in the ci')
         call symord(cr(icmo),cr(llstor),nbf(1),nob(1),ksym,kksym)
         if(iprtmo.ne.0) then
            write(iout,99088)
99088       format(//,' *** final symmetrized mcscf orbitals *** ',/)
            call wmat(cr(icmo),nbf(1),nob(1),'bflabl',' ')
         endif
      endif
c
      call iosys('write real mcscf_orbitals to rwf',nbf(1)**2,
     $            cr(icmo),0,' ')
c
      if (opnscf.eq.0) then
         cjunk='"mcscf vector"'
         call iosys('write character "transformation vector" to rwf',
     $               0,0,0,cjunk(1:16))
         call iosys('write real "mcscf vector" to rwf',nbf(1)**2,
     $               cr(icmo),0,' ')
      else
         cjunk='"scf vector"'
         call iosys('write character "transformation vector" to rwf',
     $               0,0,0,cjunk(1:16))
         call iosys('write real "scf vector" to rwf',nbf(1)**2,
     $               cr(icmo),0,' ')
         call iosys('read real "mcscf one-electron energy" from rwf',
     $               1,temp,0,' ')
         call iosys('write real "hf 1e energy" to rwf',1,temp,0,' ')
         call iosys('read real "mcscf two-electron energy" from rwf',
     $               1,temp,0,' ')
         call iosys('write real "hf 2e energy" to rwf',1,temp,0,' ')
         call iosys('read real "mcscf energy" from rwf',1,temp,0,' ')
         call iosys('write real "hf energy" to rwf',1,temp,0,' ')
      end if
c
c     ----- fix up a hf like density matrix -----
c
      num=nob(1)
      nnp=num*(num+1)/2
      ndoc=ncob(1)
      nalp=naob(1)
      noc=ndoc+nalp
c
      icmo=1
      idhf=icmo+num**2
      need=wpadti(idhf+nnp*2)
c
      call getscm(need,icr,ngot,'hfden',0)
c
      call hfden(cr(icmo),cr(idhf),num,nnp,ndoc,nalp,opnscf,ops)
c
c     ----- transform the appropriate lagrangian and density
c              to the ao basis
c
      nact=naob(1)
      ncor=ncob(1)
c
      icmo=1
      iao=icmo+num**2
      imo=iao+num**2
      it1=imo+num**2
      need=wpadti(it1+num**2)

c
      call getscm(need,icr,ngot,'lagrangian',0)
c
      call trlag(cr(icmo),cr(iao),cr(imo),cr(it1),
     $           noc,num,opnscf,ops,nact,ncor)
c
      if(iprtlag.ne.0) then
         call iosys('read real mcscf_mo_lagrangian from rwf',
     $               nbf(1)*nocc(1),cr(icmo),0,' ')
         write(iout,99099)
99099    format(//,' *** final mcscf lagrangian *** ',/)
         call wmat(cr(icmo),nbf(1),nocc(1),'bflabl',' ')
      endif
c
c     ----- copy the density matrices over to ints -----
c
      call getscm(0,icr,maxcor,'?',0)
c
      call iosys('copy "mo 1pdm" from mcscr to '//
     $           '"mcscf mo 1pdm" on moden',maxcor,cr,0,' ')
      call iosys('copy "mo 2pdm" from mcscr to '//
     $           '"mcscf mo 2pdm" on moden',maxcor,cr,0,' ')
c
c
      return
      end
