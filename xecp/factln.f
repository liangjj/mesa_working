*deck %W%   %G%
      function factln(n)
c***begin prologue     factln
c***date written       910610  (yymmdd)
c***revision date      yymmdd
c
c***keywords           
c***author             press, flannery, teulosky, and vetterling
c***source             %W%   %G%
c***purpose            compute factorial
c***description
c                                
c                      returns ln(n!) as a floating point number
c
c***references
c                      w.h.press, b.p.flannery, s.a.teulosky,
c                      and w.t.vetterling, "numerical recipes", 
c                      cambridge university press, p.180,1987.
c
c***routines called
c                      gammln
c
c***end prologue       factln
c
c
      implicit integer(a-z)
c     ----- function returned ----
      real*8 factln
c     ----- arguments unchanged -----
      integer n
c     ----- local variables -----
      real*8 a(0:99)
      real*8 zero,one
c     ----- external functions -----
      real*8 gammln
      external gammln
c
      parameter (zero=0.0d+00,one=1.0d+00)
c
      data a/100*-1.0d+00/
      save a
c
c
      if(n.lt.0) then
c        bad argument
         call lnkerr('bad argument to factln')
      else if(n.le.99) then
         if(a(n).lt.zero) then
c           fill in the table
            a(n)=gammln(n+one)
            factln=a(n)
         else
c           value already in the table
            factln=a(n)
         endif
      else
c        value out of range of table.
         factln=gammln(n+one)
      endif
c
c
      return
      end
