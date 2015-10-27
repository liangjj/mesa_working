*deck %W%   %G%
      function factrl(n)
c***begin prologue     factrl
c***date written       910610  (yymmdd)
c***revision date      yymmdd
c
c***keywords           
c***author             press, flannery, teulosky, and vetterling
c***source             %W%   %G%
c***purpose            compute factorial
c***description
c                                
c                      returns n! as a floating point number
c
c***references
c                      w.h.press, b.p.flannery, s.a.teulosky,
c                      and w.t.vetterling, "numerical recipes", 
c                      cambridge university press, p.180,1987.
c
c***routines called
c                      gammln
c
c***end prologue       factrl
c
c
      implicit integer(a-z)
c     ----- function returned -----
      real*8 factrl
c     ----- argument unchanged -----
      integer n
c     ----- local variables -----
      real*8 a(0:32),np1
      integer ntop
c     ----- external functions -----
      real*8 gammln
      external gammln
c
      data ntop/0/
      data a(0)/1.0d+00/
      save ntop,a
c
c
      if(n.lt.0) then
c        check for bad arguments
         call lnkerr('bad argument to factrl')
      else if(n.le.ntop) then
c        value already in the table
         factrl=a(n)
      else if(n.le.32) then
c        fill in table
         do 10 i=ntop+1,n
            a(i)=i*a(i-1)
   10    continue
         ntop=n
         factrl=a(n)
      else
c        this is a big number; take exp(ln(gamma function)).
         np1=float(n)+1.0d0
         factrl=exp(gammln(np1))
      endif
c
c
      return
      end
