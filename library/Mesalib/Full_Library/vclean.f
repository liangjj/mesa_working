*deck @(#)vclean.f	5.1  11/6/94
      subroutine vclean(v,toler,n)
c***begin prologue     vfill
c***date written       850601  (yymmdd)
c***revision date      yymmdd  (yymmdd)  .
c***keywords           vector, fill, load
c***author             martin,richard (lanl).
c***source             @(#)vclean.f	5.1   11/6/94
c***purpose            vector clean-up.
c***description
c                      call vclean(v,toler,n)
c                        v        output vector of length n.
c                        toler    maxium magnitude distinct from zero.
c                        n        length of vector.
c
c      this routine replaces all elements of the vector v whose absolute
c      magnitude is less than toler with zero.
c***references
c***routines called    (none)
c***end prologue       vclean
      real*8 v(n),toler,zero
      parameter (zero=0.0d+00)
c
c
      do 10 i=1,n
         if(abs(v(i)).lt.toler) v(i)=zero
   10 continue
c
c
      return
      end
