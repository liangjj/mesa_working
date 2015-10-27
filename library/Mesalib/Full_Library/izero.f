*deck @(#)izero.f	5.1  11/6/94
      subroutine izero(a,n)
c***begin prologue     izero
c***date written       850601  (yymmdd)
c***revision date      yymmdd  (yymmdd)
c***keywords           clear, zero
c***author             martin, richard (lanl)
c***source             @(#)izero.f	5.1   11/6/94
c***purpose            vectorized initialize:  a=0 .
c***description
c                      call izero(a,n)
c                        a        input vector, declared integer.
c                        n        length of a.
c
c***references
c***routines called    (none)
c***end prologue       izero
      integer a(n)
c
      do 1 i=1,n
           a(i)=0
    1 continue
c
      return
      end
