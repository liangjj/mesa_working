*deck defehl
      subroutine defehl (f, neq, t, y, h, yp, f1, f2, f3, f4, f5, ys,
     +   rpar, ipar)
c***begin prologue  defehl
c***subsidiary
c***purpose  subsidiary to derkf
c***library   slatec
c***type      single precision (defehl-s, dfehl-d)
c***author  watts, h. a., (snla)
c***description
c
c     fehlberg fourth-fifth order runge-kutta method
c **********************************************************************
c
c    defehl integrates a system of neq first order
c    ordinary differential equations of the form
c               du/dx = f(x,u)
c    over one step when the vector y(*) of initial values for u(*) and
c    the vector yp(*) of initial derivatives, satisfying  yp = f(t,y),
c    are given at the starting point x=t.
c
c    defehl advances the solution over the fixed step h and returns
c    the fifth order (sixth order accurate locally) solution
c    approximation at t+h in the array ys(*).
c    f1,---,f5 are arrays of dimension neq which are needed
c    for internal storage.
c    the formulas have been grouped to control loss of significance.
c    defehl should be called with an h not smaller than 13 units of
c    roundoff in t so that the various independent arguments can be
c    distinguished.
c
c    this subroutine has been written with all variables and statement
c    numbers entirely compatible with derkfs. for greater efficiency,
c    the call to defehl can be replaced by the module beginning with
c    line 222 and extending to the last line just before the return
c    statement.
c
c **********************************************************************
c
c***see also  derkf
c***routines called  (none)
c***revision history  (yymmdd)
c   800501  date written
c   890831  modified array declarations.  (wrb)
c   891009  removed unreferenced statement label.  (wrb)
c   891214  prologue converted to version 4.0 format.  (bab)
c   900328  added type section.  (wrb)
c   910722  updated author section.  (als)
c***end prologue  defehl
c
c
      dimension y(*),yp(*),f1(*),f2(*),f3(*),f4(*),f5(*),
     1          ys(*),rpar(*),ipar(*)
c
c***first executable statement  defehl
      ch=h/4.
      do 230 k=1,neq
  230   ys(k)=y(k)+ch*yp(k)
      call f(t+ch,ys,f1,rpar,ipar)
c
      ch=3.*h/32.
      do 240 k=1,neq
  240   ys(k)=y(k)+ch*(yp(k)+3.*f1(k))
      call f(t+3.*h/8.,ys,f2,rpar,ipar)
c
      ch=h/2197.
      do 250 k=1,neq
  250   ys(k)=y(k)+ch*(1932.*yp(k)+(7296.*f2(k)-7200.*f1(k)))
      call f(t+12.*h/13.,ys,f3,rpar,ipar)
c
      ch=h/4104.
      do 260 k=1,neq
  260   ys(k)=y(k)+ch*((8341.*yp(k)-845.*f3(k))+
     1                            (29440.*f2(k)-32832.*f1(k)))
      call f(t+h,ys,f4,rpar,ipar)
c
      ch=h/20520.
      do 270 k=1,neq
  270   ys(k)=y(k)+ch*((-6080.*yp(k)+(9295.*f3(k)-5643.*f4(k)))+
     1                             (41040.*f1(k)-28352.*f2(k)))
      call f(t+h/2.,ys,f5,rpar,ipar)
c
c     compute approximate solution at t+h
c
      ch=h/7618050.
      do 290 k=1,neq
  290   ys(k)=y(k)+ch*((902880.*yp(k)+(3855735.*f3(k)-1371249.*f4(k)))+
     1                (3953664.*f2(k)+277020.*f5(k)))
c
      return
      end
