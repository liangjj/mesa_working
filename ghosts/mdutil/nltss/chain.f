*deck %W%  %G%
 
      subroutine chain(next)
c***begin prologue     chain
c***date written       870207  (yymmdd)
c***revision date      870708  (yymmdd)
c
c     8 july 1987      pws at lanl
c         changing to a cray-ctss version. the name now has a 'c'
c         appended.
c
c***keywords           chain, links, exit
c***author             saxe, paul (lanl)
c***source             %W%   %G%
c***purpose            chain invkes the next link for mesa
c                      on the ctss cray's, this is done by writing the name
c                      of the next link to the controller (terminal)
c                      the ccl shell takes care of the rest.
c***description
c***references
c***routines called
c***end prologue       chain
c
      character*(*) next,temp*8
c
      call link("unit59=terminal//")
c
      temp='x'//next//'c'
      write(59,*)temp
c
      call exit
c
      return
      end
