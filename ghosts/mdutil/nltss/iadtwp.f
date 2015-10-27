*deck %W%  %G%
      integer function iadtwp(n)
c***begin prologue     iadtwp
c***date written       850601  (yymmdd)
c***revision date      870207  (yymmdd)
c
c   7 february 1987   pws at lanl
c      version appropriate for a 64 bit integer / 64 bit real machine
c
c***keywords           working precision, address, integer
c***author             saxe, paul (lanl)
c***source             %W%   %G%
c***purpose            core allocation transfer from integer to real address.
c***description
c                      iadtwp is an integer function used as:
c                        radd=iadtwp(n)
c                          n      integer address.
c
c***references
c***routines called    (none)
c***end prologue       iadtwp
      implicit integer(a-z)
c
      iadtwp=n
c
c32      iadtwp=(n+2)/2
c
      return
      end
