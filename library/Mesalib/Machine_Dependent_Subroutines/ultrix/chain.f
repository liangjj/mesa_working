*deck @(#)chain.f	5.1   11/6/94
      subroutine chain(next)
c***begin prologue     chain
c***date written       870207  (yymmdd)
c***revision date      yymmdd  (yymmdd)
c
c***keywords           chain, links, exit
c***author             saxe, paul (lanl)
c***source             @(#)chain.f	5.1   11/6/94
c***purpose            chain invokes the next link for mesa.
c                      on the sun's, this is done by writing the name
c                      of the next link to the standard output where
c                      the shell takes care of the rest.
c***description
c***references
c***routines called
c***end prologue       chain
c
      character*(*) next
c
      write (6,*) next
c
c
      return
      end
