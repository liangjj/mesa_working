*deck @(#)lodrte.f	5.1  11/6/94
      subroutine lodrte(opttyp,wf1typ,wf2typ,prptyp,route,numdat,
     $                  namdat)
c***begin prologue     lodrte.f
c***date written       850601  yymmdd
c***revision date      11/6/94
c
c         2 march 1991     rlm at lanl
c                  opening dat file in makrte.
c        24 july 1986      pws at lanl
c                  added mulliken population analysis (m1951)
c        10 october 1986   rlm at lanl
c                  added one-electron property integrals (m1902)
c
c***keywords           route
c***author             martin, richard (lanl)
c***source             @(#)lodrte.f	5.1   11/6/94
c***purpose            loads the skeleton route for a particular run type.
c***description
c     call lodrte(opttyp,wf1typ,wf2typ,prptyp,route,numdat,namdat)
c       opttyp  the type(name) of optimization run.
c       wf1typ  the type(name) of one-electron wavefunction run.
c       wf2typ  the type(name) of many-electron wavefunction run.
c       prptyp  the type(name) of property run.
c       route   the nonstandard route.
c       numdat  the data set unit number.
c       namdat  the data set file name.
c
c     lodrte takes as input various route types, generates a standard
c     route name, and returns a nnstandard route.
c***references
c***routines called    (none)
c***end prologue       lodrte.f
      implicit none
c     --- input variables -----
      integer numdat
      character*(*) opttyp,wf1typ,wf2typ,prptyp,namdat
c     --- input arrays (unmodified) ---
c     --- input arrays (scratch) ---
c     --- output arrays ---
      character*(*) route
c     --- output variables ---
c     --- scratch arrays ---
c     --- local variables ---
      integer pos,cskipb
      character*8 temp
      character*80 rtname,card
      logical dollar
c
c     --- generate and pack the route name.
      rtname='$'//opttyp//'/'//wf1typ//'/'//wf2typ//'/'//prptyp
      call rmvnb(rtname,rtname)
c
c     --- look for the route.
      pos=cskipb(rtname,' ')
      if(dollar(rtname(:pos),route,card,numdat)) then
      else
         temp=namdat
         call lnkerr('could not find the standard route '//rtname(:pos)
     $             //' on the data file: '//temp)
      endif
c
c
      return
      end
