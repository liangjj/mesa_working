*deck %W%  %G%
      integer function nchrpw(dummy)
c
c***begin prologue     nchrpw
c***date written       850601  (yymmdd)
c***revision date      870201  (yymmdd)
c
c   7 february 1987   pws at lanl
c       64 bit integer version
c
c***keywords           characters, word length
c***author             martin, richard (lanl)
c***source             %W%   %G%
c***purpose            returns the number of characters which can be stored
c                      in an integer word.
c***description
c                      nchrpw is an integer function used as:
c                        len=nchrpw(dummy)
c                          dummy    a dummy argument.
c
c***references
c***routines called    (none)
c***end prologue       nchrpw
c
      implicit integer(a-z)
c
c
      nchrpw=8
c
c
      return
      end
