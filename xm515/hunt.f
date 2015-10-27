*deck @(#)hunt.f	5.3 5/5/95
      subroutine hunt(xx,n,x,jlo)
c***begin prologue     hunt.f
c***date written       950502  
c***revision date      5/5/95
c
c***keywords           bisection
c***author             press, et al., numerical recipes.
c***source             @(#)hunt.f	5.3 5/5/95
c***purpose            search with correlated values. 
c***description
c   given an array xx of length n, and given a value of x, returns
c   a value jlo such that x is between xx(jlo) and xx(jlo+1).
c   xx must be monotonic, either increasing or decreasing.
c   jlo=0 or jlo=n is returned to indicate that x is out of range.
c   jlo on input is taken as the initial guess for jlo on output.
c
c***references
c
c***routines called
c
c***end prologue       hunt.f
      implicit none
c     --- input variables -----
      integer n
      real*8 x
c     --- input arrays (unmodified) ---
      real*8 xx(n)
c     --- input arrays (scratch) ---
c     --- output arrays ---
c     --- output variables ---
      integer jlo
c     --- scratch arrays ---
c     --- local variables ---
      logical ascnd
      integer jhi,jm,inc
      integer inp,iout
c
      common/io/inp,iout
c
c
c     --- see if it is an ascending list.
      ascnd=xx(n).gt.xx(1)
      if(jlo.le.0.or.jlo.gt.n) then
c        --- the input guess is not useful. go immediately to bisection.
         jlo=0
         jhi=n+1
         goto 3
      endif
c
c     --- set the hunting increment, and bracket x.
      inc=1
c     --- the operator .eqv. is true when its two logical operands
c         are either both true or both false.  the function (l1.eqv.l2)
c         is equivalent to (l1.and.l2).or.(.not.(l1.or.l2)).
      if(x.ge.xx(jlo).eqv.ascnd) then
c        --- hunt up
    1    jhi=jlo+inc
         if(jhi.gt.n) then
c           --- done hunting, we're off the end of the table.
            jhi=n+1
         else if(x.ge.xx(jhi).eqv.ascnd) then
c           --- continue hunting.
            jlo=jhi
c              --- double the increment and try again.
            inc=inc+inc
            goto 1
         endif
c        --- done hunting, value bracketed.
      else 
c        --- hunt down
         jhi=jlo
    2    jlo=jhi-inc
         if(jlo.lt.1) then
            jlo=0
         else if(x.lt.xx(jlo).eqv.ascnd) then
            jhi=jlo
            inc=inc+inc
            goto 2
         endif
      endif
c
c     --- hunt is done. begin the final bisection pass.
    3 if(jhi-jlo.ne.1) then
         jm=(jhi+jlo)/2
         if(x.gt.xx(jm).eqv.ascnd) then
            jlo=jm
         else
            jhi=jm
         endif
         goto 3
      endif
c
c
      return
      end
