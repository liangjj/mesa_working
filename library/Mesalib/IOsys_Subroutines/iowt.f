*deck @(#)iowt.f	5.1  11/6/94
      subroutine iowt(un)
c
c***begin prologue     iowt
c***date written       850125   (yymmdd)
c***revision date      860112   (yymmdd)
c***keywords           iosys dependent routines
c
c***author             saxe, paul,    (lanl)
c***source             @(#)iowt.f	5.1   11/6/94
c***purpose            to wait for asynchronous i/o completion and
c                         call lnkerr if an error occurs.
c***description        #
c
c
c***references
c
c***routines called    unit   (cftlib)
c
c   common blocks:     (none)
c
c***end prologue       iowt
c
      implicit integer (a-z)
c
c      if (unit(un)) 1,2,3
    1 continue
      return
    2 continue
    3 continue
      call lnkerr('io system: system i/o error occurred')
      end
