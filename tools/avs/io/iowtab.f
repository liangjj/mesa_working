*deck @(#)iowtab.f	4.1  7/7/93
      subroutine iowtab(un,error)
c
c***begin prologue     iowtab
c***date written       850125   (yymmdd)
c***revision date      860112   (yymmdd)
c***keywords           iosys dependent routines
c
c***author             saxe, paul,    (lanl)
c***source             @(#)iowtab.f	4.1   7/7/93
c***purpose            to wait for asynchronous i/o to finish, and
c                         return error=1 if an error occurred.
c***description        #
c
c
c***references
c
c***routines called    unit   (cftlib)
c
c   common blocks:     (none)
c
c***end prologue       iowtab
c
      implicit integer (a-z)
c
cpws      if (unit(un)) 1,2,3
    1 continue
      error=0
      return
    2 continue
    3 continue
      error=1
      return
      end
