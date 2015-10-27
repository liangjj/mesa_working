*deck @(#)wpadti.f	5.1  11/6/94
      function wpadti(n)
c***begin prologue     wpadti
c***date written       850601  (yymmdd)
c***revision date      yymmdd  (yymmdd)
c***keywords           working precision, address, integer
c***author             saxe, paul (lanl)
c***source             @(#)wpadti.f	5.1   11/6/94
c***purpose            core allocation transfer from real to integer address.
c***description
c                      wpadti is an integer function used as:
c                        iadd=wpadti(n)
c                          n      real*8 address.
c
c                        64 bit integer/ 64 bit real
c***references
c***routines called    (none)
c***end prologue       wpadti
      implicit integer(a-z)
      integer wpadti
c
c
      wpadti=n
c
c     ----- 32 bit integer/64 bit real modification -----
c     wpadti=2*n-1
c     -----
c
c
      return
      end
