*deck @(#)iadtwp.f	5.1  11/6/94
      function iadtwp(n)
c***begin prologue     iadtwp
c***date written       850601  (yymmdd)
c***revision date      870207  (yymmdd)
c
c   7 february 1987   pws at lanl
c      version appropriate for a 32 bit integer / 64 bit real machine
c
c***keywords           working precision, address, integer
c***author             saxe, paul (lanl)
c***source             @(#)iadtwp.f	5.1   11/6/94
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
      integer iadtwp
c
c***************************************
c     for 64 bit integer machines (cray)
c     iadtwp=n
c***************************************
c     for 32 bit integer machines (sun,stardent,ibm risc)
c
      iadtwp=(n+2)/2
c
      return
      end
