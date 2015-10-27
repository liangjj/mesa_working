*deck @(#)iadtwp.f	5.1  11/6/94
      function iadtwp(n)
c***begin prologue     iadtwp
c***date written       850601  (yymmdd)
c***revision date      yymmdd  (yymmdd)
c***keywords           working precision, address, integer
c***author             saxe, paul (lanl)
c***source             @(#)iadtwp.f	5.1 11/6/94 
c***purpose            core allocation transfer from integer to real address.
c***description
c                      iadtwp is an integer function used as:
c                        radd=iadtwp(n)
c                          n      integer address.
c                        64 bit integer/64 bit real version
c***references
c***routines called    (none)
c***end prologue       iadtwp
      implicit integer(a-z)
      integer iadtwp
c
      iadtwp=n
c
c     ----- 32 bit integer/64 bit real modification -----
c     iadtwp=(n+2)/2
c     -----
c
c
      return
      end
