*deck @(#)makrte.f	5.2  11/28/95
      subroutine makrte(route,ops,iout,namdat)
c***begin prologue     makrte.f
c***date written       850601  yymmdd
c***revision date      11/28/95
c
c         2 november 1995 rlm at lanl
c             adding solvent routes
c        21 january 1994  rlm at lanl
c             adding direct dft routes.
c        17 june    1993  rlm at lanl
c             adding dft routes.
c        20 january 1993  rlm at lanl
c             opening the data file as unit 1.
c             units 0 and 7 are now reserved for standard error.
c         2 march 1991    rlm at lanl
c             retrieving available options from mesadat.
c        24 march 1988    bhl at llnl
c             removed restrictions on ci optimizations 
c        24 july 1986     pws at lanl
c             added options for m1951 (mulliken populations)
c        10 october 1986     rlm at lanl
c             added options for m1902 (one-electron properties)
c
c***keywords           route, options
c***author             martin, richard (lanl)
c***source             @(#)makrte.f	5.2   11/28/95
c***purpose            sets route defaults, edits options, and prepares
c                      the skeleton route.
c***description
c     call makrte(route,ops,iout,namdat)
c       route   the nonstandard route.
c       ops     the options string
c       iout    the output file.
c       namdat  the data set file.
c
c     this module sets some defaults for the run, and looks for errors
c     in the options specified in the $route section.
c***references
c***routines called    dollar(util),lnkerr(mdutil),logkey(chr)
c                      lodrte(m1)
c***end prologue       makrte.f
      implicit none
c     --- input variables -----
      character*(*) route,ops,namdat
      integer iout
c     --- input arrays (unmodified) ---
c     --- input arrays (scratch) ---
c     --- output arrays ---
c     --- output variables ---
c     --- scratch arrays ---
c     --- local variables ---
      integer numdat,iostat,cskipb
      character*4096 opsavl
      character*80 card
      character*16 junk
      character opttyp*16,wf1typ*16,wf2typ*16,prptyp*16
      logical logkey,dollar
c
 1000 format(' unrecognized option in $route: ',a16)
c
c     --- open the data file, giving it a unit number of 1.
      numdat=1
      open (unit=numdat,file=namdat,access='sequential',
     $      form='formatted',iostat=iostat,status='old')
c
c     --- retrieve the acceptable options from the data file.
      if(dollar('options',opsavl,card,numdat)) then
      else
         junk=namdat
         call lnkerr('reference option list not found on '
     $               //junk)
      end if
c
c     --- default the skeleton route ---
      opttyp='sp'
      wf1typ='hf'
      wf2typ='-'
      prptyp='-'
c
c     --- set up the basic fields requested by the route card ---
      if (logkey(ops,'opt',.false.,' ')) then
         opttyp='murtagh-sargent'
         if (logkey(ops,'opt=fletcher-powell',.false.,' '))
     $        opttyp='fletcher-powell'
         if (logkey(ops,'opt=berny',.false.,' '))
     $         opttyp='berny'
      end if
c
      if(logkey(ops,'force-constants',.false.,' ')) then
         if(logkey(ops,'force-constants=numerical',.false.,' ')) then
            opttyp='d2e-numerical'
         else
            opttyp='d2e-analytic'
         endif
      endif
c
      if (logkey(ops,'hf',.false.,' ')) then
         if (logkey(ops,'hf=quadratic',.false.,' ')) then
            wf1typ='hf-quadratic'
         else
            wf1typ='hf'
            if(logkey(ops,'solvent',.false.,' ')) then
               wf1typ='hf-solvent'
            endif
         end if
      else if (logkey(ops,'mcscf',.false.,'mcscf')) then
         wf1typ='mcscf'
      else if (logkey(ops,'dft',.false.,'dft')) then
         wf1typ='dft'
         if(logkey(ops,'scf=directj',.false.,' ').or.
     $        logkey(ops,'scf=poissonj',.false.,' ')) then
            wf1typ='dft-direct'
         endif
         if(logkey(ops,'parallel',.false.,' ')) then
            wf1typ=wf1typ(1:cskipb(wf1typ,' '))//'-p'
         endif
         if(logkey(ops,'solvent',.false.,' ')) then
            wf1typ=wf1typ(1:cskipb(wf1typ,' '))//'-solvent'
         endif
      end if
c
      if (logkey(ops,'ci',.false.,' ')) then
         wf2typ='ci'
         if (logkey(ops,'ci=full',.false.,' ')) wf2typ='full-ci'
         if (logkey(ops,'ci=formula',.false.,' ')) wf2typ='formula-ci'
      end if
c
      if (logkey(ops,'properties',.false.,'properties')) then
         prptyp='full'
      end if
c
c     --- check for inconsistencies ---
      if(wf1typ.eq.'dft') then
         if(wf2typ.ne.'-') then
            call lnkerr('cannot do dft with ci')
         endif
      endif
c
c     --- retrieve the skeleton route from disk ---
      call lodrte(opttyp,wf1typ,wf2typ,prptyp,route,numdat,namdat)
c
c     --- edit the options list ---
c         check for errors in route.
      if((index(ops,'coord').ne.0).and.(opttyp.ne.'sp')) then
         if (opttyp.eq.'d2e-analytic'.or.opttyp.eq.'d2e-numerical') then
         else
            call lnkerr('coord cannot be used for optimizations.')
         end if
      endif
c
c
      return
      end
