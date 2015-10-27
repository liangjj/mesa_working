*deck dfehl
      subroutine dfehl (df, neq, t, y, h, yp, f1, f2, f3, f4, f5, ys,
     +   rpar, ipar)
c***begin prologue  dfehl
c***subsidiary
c***purpose  subsidiary to dderkf
c***library   slatec
c***type      double precision (defehl-s, dfehl-d)
c***author  watts, h. a., (snla)
c***description
c
c     fehlberg fourth-fifth order runge-kutta method
c **********************************************************************
c
c    dfehl integrates a system of neq first order
c    ordinary differential equations of the form
c               du/dx = df(x,u)
c    over one step when the vector y(*) of initial values for u(*) and
c    the vector yp(*) of initial derivatives, satisfying  yp = df(t,y),
c    are given at the starting point x=t.
c
c    dfehl advances the solution over the fixed step h and returns
c    the fifth order (sixth order accurate locally) solution
c    approximation at t+h in the array ys(*).
c    f1,---,f5 are arrays of dimension neq which are needed
c    for internal storage.
c    the formulas have been grouped to control loss of significance.
c    dfehl should be called with an h not smaller than 13 units of
c    roundoff in t so that the various independent arguments can be
c    distinguished.
c
c    this subroutine has been written with all variables and statement
c    numbers entirely compatible with drkfs. for greater efficiency,
c    the call to dfehl can be replaced by the module beginning with
c    line 222 and extending to the last line just before the return
c    statement.
c
c **********************************************************************
c
c***see also  dderkf
c***routines called  (none)
c***revision history  (yymmdd)
c   820301  date written
c   890831  modified array declarations.  (wrb)
c   891214  prologue converted to version 4.0 format.  (bab)
c   900328  added type section.  (wrb)
c   910722  updated author section.  (als)
c***end prologue  dfehl
c
      integer ipar, k, neq
      double precision ch, f1, f2, f3, f4, f5, h, rpar, t, y, yp, ys
      dimension y(*),yp(*),f1(*),f2(*),f3(*),f4(*),f5(*),
     1          ys(*),rpar(*),ipar(*)
c
c***first executable statement  dfehl
      ch = h/4.0d0
      do 10 k = 1, neq
         ys(k) = y(k) + ch*yp(k)
   10 continue
      call df(t+ch,ys,f1,rpar,ipar)
c
      ch = 3.0d0*h/32.0d0
      do 20 k = 1, neq
         ys(k) = y(k) + ch*(yp(k) + 3.0d0*f1(k))
   20 continue
      call df(t+3.0d0*h/8.0d0,ys,f2,rpar,ipar)
c
      ch = h/2197.0d0
      do 30 k = 1, neq
         ys(k) = y(k)
     1           + ch
     2             *(1932.0d0*yp(k) + (7296.0d0*f2(k) - 7200.0d0*f1(k)))
   30 continue
      call df(t+12.0d0*h/13.0d0,ys,f3,rpar,ipar)
c
      ch = h/4104.0d0
      do 40 k = 1, neq
         ys(k) = y(k)
     1           + ch
     2             *((8341.0d0*yp(k) - 845.0d0*f3(k))
     3               + (29440.0d0*f2(k) - 32832.0d0*f1(k)))
   40 continue
      call df(t+h,ys,f4,rpar,ipar)
c
      ch = h/20520.0d0
      do 50 k = 1, neq
         ys(k) = y(k)
     1           + ch
     2             *((-6080.0d0*yp(k)
     3                + (9295.0d0*f3(k) - 5643.0d0*f4(k)))
     4               + (41040.0d0*f1(k) - 28352.0d0*f2(k)))
   50 continue
      call df(t+h/2.0d0,ys,f5,rpar,ipar)
c
c     compute approximate solution at t+h
c
      ch = h/7618050.0d0
      do 60 k = 1, neq
         ys(k) = y(k)
     1           + ch
     2             *((902880.0d0*yp(k)
     3                + (3855735.0d0*f3(k) - 1371249.0d0*f4(k)))
     4               + (3953664.0d0*f2(k) + 277020.0d0*f5(k)))
   60 continue
c
      return
      end
