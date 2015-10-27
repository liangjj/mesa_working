*deck @(#)vfill.f	5.1  11/6/94
      subroutine vfill(v,s,n)
c***begin prologue     vfill
c***date written       850601  (yymmdd)
c***revision date      yymmdd  (yymmdd)  .
c***keywords           vector, fill, load
c***author             saxe, paul (lanl).
c***source             @(#)vfill.f	5.1   11/6/94
c***purpose            vectorized scalar load:  v=s.
c***description
c                      call vfill(v,s,n)
c                        v        output vector of length n.
c                        s        scalar to load.
c                        n        length of vector.
c
c***references
c***routines called    (none)
c***end prologue       vfill
      real*8 v(n),s
      do 1 i=1,n
         v(i)=s
    1 continue
      return
      end
