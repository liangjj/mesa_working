*deck %W%  %G%
      integer function wpadti(n)
c***begin prologue     wpadti
c***date written       850601  (yymmdd)
c***revision date      870207  (yymmdd)
c
c   7 february 1987    pws at lanl
c        64 bit integer / 64 bit real version
c
c***keywords           working precision, address, integer
c***author             saxe, paul (lanl)
c***source             %W%   %G%
c***purpose            core allocation transfer from real to integer address.
c***description
c                      wpadti is an integer function used as:
c                        iadd=wptoin(n)
c                          n      real*8 address.
c
c***references
c***routines called    (none)
c***end prologue       wpadti
      implicit integer(a-z)
c
c
      wpadti=n
c
c32       wpadti=2*n-1
c
      return
      end
