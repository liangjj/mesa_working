*deck @(#)kohnop.f	1.1 9/8/91
c***begin prologue     m6009
c***date written       880423   (yymmdd)
c***revision date      890620   (yymmdd)
c***                   extensive revision to make compatible with polyatomic
c***                   calculations. all single center matrix elements
c***                   removed. they are now calculated in m6005 and final
c***                   matrices assembled in m6008.             
c***keywords           m6009, link 6009, kohn variational
c***author             schneider, barry (lanl), rescigno, tom(llnl)  
c***source             m6009
c***purpose            perform kohn variational calculations
c***description        muti-channel kohn variational calculations using
c***                   separable exchange and optical potentials for
c***                   molecules.
c
c***references         schneider and rescigno, physical review
c
c***routines called    iosys, util and mdutil
c***end prologue       m6009
      program kohnop
      implicit integer(a-z)
      parameter (dime=100)
      logical logkey, logky
      real *8 z, dumr
      complex *16 dumc
      real *8 energy, eau, eev, egrnd
      character *4096 ops
      character *8 cpass, typem
      character *128 filint, filtmt
      character *1600 card
      character *24 filnm
      character *15 type
      character *3 ans, copt
      character *16 fptoc
      common a(1)
      dimension z(1)
      common /io/ inp, iout
      common / memory / ioff
      dimension energy(dime), logky(10)
      equivalence (z,a)
      call drum
      write (iout,1)
c----------------------------------------------------------------------c
c                  recover options string                              c
c----------------------------------------------------------------------c
      call iosys ('read character options from rwf',-1,0,0,ops)
      logky(1)=logkey(ops,'m6009=no-partitioning',.false.,' ')
      logky(2)=logkey(ops,'m6009=real',.false.,' ')
      logky(3)=logkey(ops,'print=m6009=all',.false.,' ')
      logky(4)=logkey(ops,'print=m6009=complex-rhs',.false.,' ')
      logky(5)=logkey(ops,'print=m6009=new-rhs',.false.,' ')
      logky(6)=logkey(ops,'print=m6009=hfull',.false.,' ')
      logky(7)=logkey(ops,'print=m6009=full-solution',.false.,' ')
      logky(8)=logkey(ops,'print=m6009=complex-numerator',.false.,' ')
      logky(9)=logkey(ops,'print=m6009=rhs',.false.,' ')
      if (logky(3)) then
          do 100 i=4,9
             logky(i)=.true.
  100     continue
      endif
c----------------------------------------------------------------------c
c                  position input file                                 c 
c            read in title ,  filenames and                            c
c                   various parameters                                 c
c----------------------------------------------------------------------c
      call posinp ('$kohnopt',cpass)
      call cardin (card)
      if (logky(1)) then
          write (iout,2)
      else
          write (iout,3)
      endif
      if (logky(2)) then
          write (iout,4)
      else
          write (iout,5)
      endif
      call iosys ('read character "kohn integral filename" from '//
     1            'rwf',-1,0,0,filint)
      call iosys ('open kohnint as old',0,0,0,filint)
      call iosys ('read character "kohn tmatrix filename" from rwf',
     1             -1,0,0,filtmt)
      call iosys ('open tmat as unknown',262144,0,0,filtmt) 
c----------------------------------------------------------------------c
c          read in information to determine memory for calculation     c
c----------------------------------------------------------------------c
      call iosys ('read integer "no. energies" from kohnint',1,nen,
     1            0,' ')
      call iosys ('read real "scatt energies" from kohnint',nen,energy,
     1            0,' ')
      call iosys ('read integer "total channels" from kohnint',1,ntchn,
     1            0,' ')
      call iosys ('read integer "total bound" from kohnint',1,matbv,0,
     1            ' ')
      call iosys ('read real "chan energies" from kohnint',1,egrnd,0,
     1            ' ')
      call iosys ('does "complex optical potential" exist on kohnint',
     1             -1,0,0,ans)
      typem='real'
      fac=1
      copt='no'
      if (ans.eq.'yes') then
          copt='yes'
          fac=2
          typem='complex'
          write(iout,50)
      else
          write(iout,60)
      endif
      write (iout,6) ntchn, matbv
c----------------------------------------------------------------------c
c                    lets get the needed memory                        c
c----------------------------------------------------------------------c
      call iosys ('read integer maxsiz from rwf',1,maxcor,0,' ')
      ioff=1
      do 10 trips=1,2
c----------------------------------------------------------------------c 
c                    parcel out the arrays                             c
c----------------------------------------------------------------------c
         hpp=ioff
         hpm=hpp+2*ntchn*ntchn
         hpb=hpm+2*ntchn*ntchn
c----------------------------------------------------------------------c
c                leave enough room beginning at vrr to store           c
c                two real or two complex variables                     c
c                depending on the presence of a real or complex        c
c                           optical potential                          c
c----------------------------------------------------------------------c
         vrr=hpb+2*ntchn*matbv
         vrrc=vrr
         if (fac.eq.2) then 
             vrrc=vrr+fac*ntchn*ntchn
         endif
         vii=vrrc+fac*ntchn*ntchn
         vir=vii+2*ntchn*ntchn
         vri=vir+2*ntchn*ntchn
         tvir=vri
         vfrb=vri+2*ntchn*ntchn
         vfib=vfrb+ntchn*matbv
         vbfr=vfib+2*ntchn*matbv
c----------------------------------------------------------------------c
c                leave enough room beginning at vbfr to store          c
c                two real or two complex variables                     c
c                depending on the presence of a real or complex        c
c                           optical potential                          c
c                                                                      c
c                vbfr is used as a real or a complex array             c
c                depending on the case after the call to               c
c                             slvlin                                   c
c----------------------------------------------------------------------c
         vbfrt=vbfr+fac*ntchn*matbv
         vbfi=vbfrt+fac*ntchn*matbv
c----------------------------------------------------------------------c
c             hambb is used as a real or complex array depending       c
c             depending on the presence of a real or complex           c
c                          optical potential                           c
c----------------------------------------------------------------------c

         hambb=vbfi+2*ntchn*matbv
c----------------------------------------------------------------------c
c    branch here depending on whether you partition off the            c
c                        bound-bound matrix                            c
c----------------------------------------------------------------------c
         if (logky(1)) then
             ncomp=matbv+ntchn
             matcmp=hambb+fac*matbv*matbv
             crhs=matcmp+2*ncomp*ncomp
             clhs=crhs+2*ncomp*ntchn
             hld=clhs+2*ncomp*ntchn
             ipvt=wpadti(hld)
         else
             hld=hambb+fac*matbv*matbv      
             ipvt=wpadti(hld)
         endif
c     for extraction of t-matrix
         tmat=hld+add
         seig=tmat+2*ntchn*ntchn
         svec=seig+2*ntchn
         dum=svec+2*ntchn*ntchn
         words=dum+3*ntchn
         if (trips.eq.1) then
             if (maxcor.lt.words) then
                 call lnkerr ('cannot get enough memory:will quit')
             endif
            call iosys ('write integer maxsiz to rwf',1,words,0,' ')
            call getscm(words,z(1),ngot,'kohn',0)
            write (iout,7) words
         endif
   10 continue
c----------------------------------------------------------------------c
c              memory obtained: lets do calculation                    c
c----------------------------------------------------------------------c
c----------------------------------------------------------------------c
c                  begin energy dependent step                         c
c----------------------------------------------------------------------c
      do 20 ene=1,nen
         eau=.5d+00*energy(ene)
         eev=eau*27.21d0
         write (iout,8) eau, eev
c----------------------------------------------------------------------c
c               read in the bound-bound hamiltonian matrix             c 
c               containing exchange and correlation                    c
c----------------------------------------------------------------------c
         call rddeno(z(hambb),z(hambb),energy(ene),matbv,copt)
c----------------------------------------------------------------------c
c               read in all free-free and bound-free matrix            c
c                           elements                                   c
c----------------------------------------------------------------------c
         filnm='hpp-'//fptoc(energy(ene))
         call iosys ('read real '//filnm//' from kohnint',
     1               2*ntchn*ntchn,z(hpp),0,' ')
         filnm='hpm-'//fptoc(energy(ene))
         call iosys ('read real '//filnm//' from kohnint',
     1               2*ntchn*ntchn,z(hpm),0,' ')
         filnm='hpb-'//fptoc(energy(ene))
         call iosys ('read real '//filnm//' from kohnint',
     1               2*ntchn*matbv,z(hpb),0,' ')
c----------------------------------------------------------------------c
c                  get the derived matrices                            c
c----------------------------------------------------------------------c
         call mkfree(z(hpp),z(hpm),z(vrr),z(vii),z(vri),z(vir),ntchn) 
         call mkfrbn(z(hpb),z(vfrb),z(vfib),z(vbfr),z(vbfi),matbv,ntchn)
c----------------------------------------------------------------------c
c          store vbfr in a real or complex array depending on          c
c                            typem                                     c
c----------------------------------------------------------------------c
         call cvbfr(z(vbfr),z(vbfrt),z(vbfrt),matbv,ntchn,typem)
c----------------------------------------------------------------------c
c           construct the final matrix and right hand side             c 
c----------------------------------------------------------------------c
         if (.not.logky(1)) then
c----------------------------------------------------------------------c
c           solve the linear equations using partitioning              c
c----------------------------------------------------------------------c
             call slvlin(z(hambb),z(hambb),a(ipvt),dumr,dumc,matbv,
     1                   ntchn,'factor',typem,'complex',logky(9))
             call slvlin(z(hambb),z(hambb),a(ipvt),z(vbfi),z(vbfi),
     1                   matbv,ntchn,'no factor',typem,'complex',
     2                   logky(9))
c----------------------------------------------------------------------c
             if (typem.eq.'real') then
                 call slvlin(z(hambb),dumc,a(ipvt),z(vbfrt),dumc,matbv,
     1                       ntchn,'no factor',typem,'real',logky(9))
             elseif (typem.eq.'complex') then
                 call slvlin(dumr,z(hambb),a(ipvt),dumr,z(vbfrt),matbv,
     1                       ntchn,'no factor',typem,'real',logky(9))
             endif
c----------------------------------------------------------------------c
c                     copy it back to vbfr                             c
c----------------------------------------------------------------------c
             call copy(z(vbfrt),z(vbfr),fac*matbv*ntchn)
c----------------------------------------------------------------------c
c                  make effective free-free matrix                     c 
c----------------------------------------------------------------------c
             call ffmat(z(vii),z(vfib),z(vbfi),dumr,ntchn,matbv,
     1                  'complex',logky(5))
             call ffmat(z(vir),z(vfib),z(vbfr),z(vbfr),ntchn,matbv,
     1                  typem,logky(5))
             call matm(z(vir),ntchn*ntchn)
c----------------------------------------------------------------------c
c                    solve complex linear system                       c
c                    to get free coefficients                          c
c----------------------------------------------------------------------c
             filnm='mpp-'//fptoc(energy(ene))
             call clvcmp(z(vii),a(ipvt),z(vir),ntchn,ntchn,'factor',
     1                   logky(4),filnm)
             call cc2opy(z(vir),z(tvir),ntchn*ntchn)
             filnm='rpp-'//fptoc(energy(ene))
             call clvcmp(z(vii),a(ipvt),z(vir),ntchn,ntchn,'no factor',
     1                   logky(4),filnm)
c----------------------------------------------------------------------c
c                extract non variational t matrices                    c
c----------------------------------------------------------------------c
             call tnonvr (z(vir),z(tmat),z(seig),z(svec),z(dum),
     1                    logky(2),ntchn)
c----------------------------------------------------------------------c
c              extract variational t matrices                          c
c----------------------------------------------------------------------c
             type='partitioned'
             filnm='tmat-'//fptoc(energy(ene))
             call tvar(z(vrr),z(vrrc),z(vir),z(tvir),z(vfrb),z(vbfr),
     1                 z(vbfr),z(tmat),z(seig),z(svec),z(dum),logky(2),
     2                 type,ntchn,matbv,filnm,typem)
c----------------------------------------------------------------------c
c           solve the linear equations without partitioning            c
c           and extract required t matrices.                           c
c----------------------------------------------------------------------c
         else
             call frmcmp(z(vrr),z(vrrc),z(vii),z(vir),z(vri),z(vbfr),
     1                   z(vbfr),z(vbfi),z(vfrb),z(vfib),z(hambb),
     2                   z(hambb),z(matcmp),z(crhs),z(clhs),a(ipvt),
     3                   z(tmat),z(seig),z(svec),z(dum),logky(2),
     4                   matbv,ntchn,ncomp,logky(6),logky(7),logky(8),
     5                   filnm,typem)
         endif 
   20 continue
      call iosys ('rewind all on kohnint read-and-write',0,0,0,' ')
      call iosys ('close kohnint',0,0,0,' ')
      call iosys ('rewind all on tmat read-and-write',0,0,0,' ')
      call iosys ('close tmat',0,0,0,' ')
      call iosys ('write integer maxsiz to rwf',1,maxcor,0,' ')
      call chainx (0)
      stop
    1 format (//,20x,'***** m6009:kohn variational scattering program **
     1***')
    2 format (/,15x,'calculation performed without partitioning')
    3 format (/,15x,'calculation performed with partitioning')
    4 format (/,5x,'use the kohn with real boundary conditions')
    5 format (/,5x,'use the kohn with complex boundary conditions')
    6 format (/,5x,'no. channels',1x,i3,2x,'size of bound-bound matrix',
     1        1x,i4)
    7 format (/,5x,i8,1x,'words obtained by getscm')
    8 format (/,5x,'incident electron energy',1x,f15.8,'(hartrees)',
     1        1x,f15.8,'(ev)')
   50 format(//,15x,'using a complex optical potential')  
   60 format(//,15x,'using a real optical potential')  
      end
