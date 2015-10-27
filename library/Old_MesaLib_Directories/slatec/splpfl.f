*deck splpfl
      subroutine splpfl (mrelas, nvars, ienter, ileave, ibasis, ind,
     +   ibb, theta, dirnrm, rprnrm, csc, ww, bl, bu, erp, rprim,
     +   primal, finite, zerolv)
c***begin prologue  splpfl
c***subsidiary
c***purpose  subsidiary to splp
c***library   slatec
c***type      single precision (splpfl-s, dplpfl-d)
c***author  (unknown)
c***description
c
c     the editing required to convert this subroutine from single to
c     double precision involves the following character string changes.
c
c     use an editing command (change) /string-1/(to)string-2/.
c     /real (12 blanks)/double precision/.
c
c     this subprogram is part of the splp( ) package.
c     it implements the procedure (choose variable to leave basis).
c     revised 811130-1045
c     revised yymmdd-hhmm
c
c***see also  splp
c***routines called  (none)
c***revision history  (yymmdd)
c   811215  date written
c   890531  changed all specific intrinsics to generic.  (wrb)
c   890605  removed unreferenced labels.  (wrb)
c   891214  prologue converted to version 4.0 format.  (bab)
c   900328  added type section.  (wrb)
c***end prologue  splpfl
      integer ibasis(*),ind(*),ibb(*)
      real             csc(*),ww(*),bl(*),bu(*),erp(*),rprim(*),
     * primal(*),bound,dirnrm,ratio,rprnrm,theta,zero
      logical finite,zerolv
c***first executable statement  splpfl
      zero=0.e0
c
c     see if the entering variable is restricting the step length
c     because of an upper bound.
      finite=.false.
      j=ibasis(ienter)
      if (.not.(ind(j).eq.3)) go to 20002
      theta=bu(j)-bl(j)
      if(j.le.nvars)theta=theta/csc(j)
      finite=.true.
      ileave=ienter
c
c     now use the basic variables to possibly restrict the step
c     length even further.
20002 i=1
      n20005=mrelas
      go to 20006
20005 i=i+1
20006 if ((n20005-i).lt.0) go to 20007
      j=ibasis(i)
c
c     if this is a free variable, do not use it to
c     restrict the step length.
      if (.not.(ind(j).eq.4)) go to 20009
      go to 20005
c
c     if direction component is about zero, ignore it for computing
c     the step length.
20009 if (.not.(abs(ww(i)).le.dirnrm*erp(i))) go to 20012
      go to 20005
20012 if (.not.(ww(i).gt.zero)) go to 20015
c
c     if rprim(i) is essentially zero, set ratio to zero and exit loop.
      if (.not.(abs(rprim(i)).le.rprnrm*erp(i))) go to 20018
      theta=zero
      ileave=i
      finite=.true.
      go to 20008
c
c     the value of rprim(i) will decrease only to its lower bound or
c     only to its upper bound.  if it decreases to its
c     upper bound, then rprim(i) has already been translated
c     to its upper bound and nothing needs to be done to ibb(j).
20018 if (.not.(rprim(i).gt.zero)) go to 10001
      ratio=rprim(i)/ww(i)
      if (.not.(.not.finite)) go to 20021
      ileave=i
      theta=ratio
      finite=.true.
      go to 20022
20021 if (.not.(ratio.lt.theta)) go to 10002
      ileave=i
      theta=ratio
10002 continue
20022 continue
      go to 20019
c
c     the value rprim(i).lt.zero will not restrict the step.
10001 continue
c
c     the direction component is negative, therefore the variable will
c     increase.
20019 go to 20016
c
c     if the variable is less than its lower bound, it can
c     increase only to its lower bound.
20015 if (.not.(primal(i+nvars).lt.zero)) go to 20024
      ratio=rprim(i)/ww(i)
      if (ratio.lt.zero) ratio=zero
      if (.not.(.not.finite)) go to 20027
      ileave=i
      theta=ratio
      finite=.true.
      go to 20028
20027 if (.not.(ratio.lt.theta)) go to 10003
      ileave=i
      theta=ratio
10003 continue
20028 continue
c
c     if the basic variable is feasible and is not at its upper bound,
c     then it can increase to its upper bound.
      go to 20025
20024 if (.not.(ind(j).eq.3 .and. primal(i+nvars).eq.zero)) go to 10004
      bound=bu(j)-bl(j)
      if(j.le.nvars) bound=bound/csc(j)
      ratio=(bound-rprim(i))/(-ww(i))
      if (.not.(.not.finite)) go to 20030
      ileave=-i
      theta=ratio
      finite=.true.
      go to 20031
20030 if (.not.(ratio.lt.theta)) go to 10005
      ileave=-i
      theta=ratio
10005 continue
20031 continue
      continue
10004 continue
20025 continue
20016 go to 20005
20007 continue
c
c     if step length is finite, see if step length is about zero.
20008 if (.not.(finite)) go to 20033
      zerolv=.true.
      i=1
      n20036=mrelas
      go to 20037
20036 i=i+1
20037 if ((n20036-i).lt.0) go to 20038
      zerolv=zerolv.and. abs(theta*ww(i)).le.erp(i)*rprnrm
      if (.not.(.not. zerolv)) go to 20040
      go to 20039
20040 go to 20036
20038 continue
20039 continue
20033 continue
      return
      end
