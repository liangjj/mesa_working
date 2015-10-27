*deck @(#)stderr.f	5.1 11/6/94
      function stderr()
c***begin prologue     stderr
c***date written       930201  (yymmdd)
c***revision date              (yymmdd)
c
c***keywords           standard error, i/o
c***author             martin, richard (lanl)
c***source             @(#)stderr.f	5.1 11/6/94
c***purpose            return the unit number corresponding to the
c                      pre-assigned standard error
c***description
c***references
c***routines called    (none)
c***end prologue       stderr
      implicit integer(a-z)
      integer stderr
c
      stderr=0
c
      return
      end
