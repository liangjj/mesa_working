*deck %W%   %G%
      function bico(n,k)
c***begin prologue     bico
c***date written       910610  (yymmdd)
c***revision date      yymmdd
c
c***keywords           
c***author             press, flannery, teulosky, and vetterling
c***source             %W%   %G%
c***purpose            compute binomial coefficient
c***description
c                                
c                      returns the binomial coefficient (n
c                                                        k)
c
c***references
c                      w.h.press, b.p.flannery, s.a.teulosky,
c                      and w.t.vetterling, "numerical recipes", 
c                      cambridge university press, p.180,1987.
c
c***routines called
c                      factln
c                      anint  ...converts floating point to nearest integer
c                                and stores as floating point
c
c***end prologue       bico
c
c
      implicit integer(a-z)
c     ----- function returned -----
      real*8 bico
c     ----- arguments unchanged -----
      integer n,k
c     ----- external functions -----
      real*8 factln,anint,exp
c
c     test for bad arguments
      if(k.gt.n) then
         call lnkerr('bad arguments to bico:')
      endif
c
      bico=anint(exp(factln(n)-factln(k)-factln(n-k)))
c
c
      return
      end
