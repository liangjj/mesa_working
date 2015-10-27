*deck %W%   %G%
      function gammln(xx)
c***begin prologue     gammln
c***date written       910610  (yymmdd)
c***revision date      yymmdd
c
c***keywords           
c***author             press, flannery, teulosky, and vetterling
c***source             %W%   %G%
c***purpose            compute gamma function
c***description
c                                
c                      returns ln(gamma(xx)) for xx>0.
c                      full accuracy is obtained for xx>1.
c                      for 0<xx<1, use a reflection formula
c                        gamma(1-z)=pi/(gamma(z)*sin(pi*z))
c                                  =pi*z/(gamma(1+z)*sin(pi*z))
c
c***references
c                      w.h.press, b.p.flannery, s.a.teulosky,
c                      and w.t.vetterling, "numerical recipes", 
c                      cambridge university press, p.180,1987.
c
c***routines called
c                      none
c
c***end prologue       gammln
c
c
      implicit integer(a-z)
c     ----- function returned -----
      real*8 gammln
c     ----- arguments unchanged -----
      real*8 xx
c     ----- local variables -----
      real*8 cof(6)
      real*8 stp,zero,half,one,fpf,x,tmp,ser
c     ----- external functions -----
      real*8 log
c
      data cof/76.18009173d0, -86.50532033d0,    24.01409822d0,
     $         -1.231739516d0,  0.120858003d-02, -0.536382d-05/
      data stp/2.50662827465d0/
      save cof,stp
c
      parameter (zero=0.0d+00,half=0.5d+00,one=1.0d+00,fpf=5.5d+00)
c
c
      if(xx.le.zero) then
c        check for bad arguments
         call lnkerr('bad argument to gammln')
      endif
c
      x=xx-one
      tmp=x+fpf
      tmp=(x+half)*log(tmp)-tmp
      ser=one
      do 10 i=1,6
         x=x+one
         ser=ser+cof(i)/x
   10 continue
      gammln=tmp+log(stp*ser)
c
c
      return
      end
