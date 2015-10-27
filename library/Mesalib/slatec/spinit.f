*deck spinit
      subroutine spinit (mrelas, nvars, costs, bl, bu, ind, primal,
     +   info, amat, csc, costsc, colnrm, xlamda, anorm, rhs, rhsnrm,
     +   ibasis, ibb, imat, lopt)
c***begin prologue  spinit
c***subsidiary
c***purpose  subsidiary to splp
c***library   slatec
c***type      single precision (spinit-s, dpinit-d)
c***author  (unknown)
c***description
c
c     the editing required to convert this subroutine from single to
c     double precision involves the following character string changes.
c
c     use an editing command (change) /string-1/(to)string-2/.
c     /real (12 blanks)/double precision/,/scopy/dcopy/
c     revised 810519-0900
c     revised yymmdd-hhmm
c
c     initialization subroutine for splp(*) package.
c
c***see also  splp
c***routines called  pnnzrs, sasum, scopy
c***revision history  (yymmdd)
c   811215  date written
c   890531  changed all specific intrinsics to generic.  (wrb)
c   890605  removed unreferenced labels.  (wrb)
c   891214  prologue converted to version 4.0 format.  (bab)
c   900328  added type section.  (wrb)
c***end prologue  spinit
      real             aij,amat(*),anorm,bl(*),bu(*),cmax,
     * colnrm(*),costs(*),costsc,csc(*),csum,one,primal(*),
     * rhs(*),rhsnrm,scalr,testsc,xlamda,zero
      integer ibasis(*),ibb(*),imat(*),ind(*)
      logical contin,usrbas,colscp,cstscp,minprb,lopt(8)
c
c***first executable statement  spinit
      zero=0.
      one=1.
      contin=lopt(1)
      usrbas=lopt(2)
      colscp=lopt(5)
      cstscp=lopt(6)
      minprb=lopt(7)
c
c     scale data. normalize bounds. form column check sums.
      go to 30001
c
c     initialize active basis matrix.
20002 continue
      go to 30002
20003 return
c
c     procedure (scale data. normalize bounds. form column check sums)
c
c     do column scaling if not provided by the user.
30001 if (.not.(.not. colscp)) go to 20004
      j=1
      n20007=nvars
      go to 20008
20007 j=j+1
20008 if ((n20007-j).lt.0) go to 20009
      cmax=zero
      i=0
20011 call pnnzrs(i,aij,iplace,amat,imat,j)
      if (.not.(i.eq.0)) go to 20013
      go to 20012
20013 continue
      cmax=max(cmax,abs(aij))
      go to 20011
20012 if (.not.(cmax.eq.zero)) go to 20016
      csc(j)=one
      go to 20017
20016 csc(j)=one/cmax
20017 continue
      go to 20007
20009 continue
c
c     form check sums of columns. compute matrix norm of scaled matrix.
20004 anorm = zero
      j=1
      n20019=nvars
      go to 20020
20019 j=j+1
20020 if ((n20019-j).lt.0) go to 20021
      primal(j)=zero
      csum = zero
      i=0
20023 call pnnzrs(i,aij,iplace,amat,imat,j)
      if (.not.(i.le.0)) go to 20025
      go to 20024
20025 continue
      primal(j)=primal(j)+aij
      csum = csum+abs(aij)
      go to 20023
20024 if (ind(j).eq.2) csc(j)=-csc(j)
      primal(j)=primal(j)*csc(j)
      colnrm(j)=abs(csc(j)*csum)
      anorm = max(anorm,colnrm(j))
      go to 20019
c
c     if the user has not provided cost vector scaling then scale it
c     using the max. norm of the transformed cost vector, if nonzero.
20021 testsc=zero
      j=1
      n20028=nvars
      go to 20029
20028 j=j+1
20029 if ((n20028-j).lt.0) go to 20030
      testsc=max(testsc,abs(csc(j)*costs(j)))
      go to 20028
20030 if (.not.(.not.cstscp)) go to 20032
      if (.not.(testsc.gt.zero)) go to 20035
      costsc=one/testsc
      go to 20036
20035 costsc=one
20036 continue
      continue
20032 xlamda=(costsc+costsc)*testsc
      if (xlamda.eq.zero) xlamda=one
c
c     if maximization problem, then change sign of costsc and lamda
c     =weight for penalty-feasibility method.
      if (.not.(.not.minprb)) go to 20038
      costsc=-costsc
20038 go to 20002
c:ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     procedure (initialize rhs(*),ibasis(*), and ibb(*))
c
c     initially set right-hand side vector to zero.
30002 call scopy(mrelas,zero,0,rhs,1)
c
c     translate rhs according to classification of independent variables
      j=1
      n20041=nvars
      go to 20042
20041 j=j+1
20042 if ((n20041-j).lt.0) go to 20043
      if (.not.(ind(j).eq.1)) go to 20045
      scalr=-bl(j)
      go to 20046
20045 if (.not.(ind(j).eq.2)) go to 10001
      scalr=-bu(j)
      go to 20046
10001 if (.not.(ind(j).eq.3)) go to 10002
      scalr=-bl(j)
      go to 20046
10002 if (.not.(ind(j).eq.4)) go to 10003
      scalr=zero
10003 continue
20046 continue
      if (.not.(scalr.ne.zero)) go to 20048
      i=0
20051 call pnnzrs(i,aij,iplace,amat,imat,j)
      if (.not.(i.le.0)) go to 20053
      go to 20052
20053 continue
      rhs(i)=scalr*aij+rhs(i)
      go to 20051
20052 continue
20048 continue
      go to 20041
c
c     translate rhs according to classification of dependent variables.
20043 i=nvars+1
      n20056=nvars+mrelas
      go to 20057
20056 i=i+1
20057 if ((n20056-i).lt.0) go to 20058
      if (.not.(ind(i).eq.1)) go to 20060
      scalr=bl(i)
      go to 20061
20060 if (.not.(ind(i).eq.2)) go to 10004
      scalr=bu(i)
      go to 20061
10004 if (.not.(ind(i).eq.3)) go to 10005
      scalr=bl(i)
      go to 20061
10005 if (.not.(ind(i).eq.4)) go to 10006
      scalr=zero
10006 continue
20061 continue
      rhs(i-nvars)=rhs(i-nvars)+scalr
      go to 20056
20058 rhsnrm=sasum(mrelas,rhs,1)
c
c     if this is not a continuation or the user has not provided the
c     initial basis, then the initial basis is comprised of the
c     dependent variables.
      if (.not.(.not.(contin .or. usrbas))) go to 20063
      j=1
      n20066=mrelas
      go to 20067
20066 j=j+1
20067 if ((n20066-j).lt.0) go to 20068
      ibasis(j)=nvars+j
      go to 20066
20068 continue
c
c     define the array ibb(*)
20063 j=1
      n20070=nvars+mrelas
      go to 20071
20070 j=j+1
20071 if ((n20070-j).lt.0) go to 20072
      ibb(j)=1
      go to 20070
20072 j=1
      n20074=mrelas
      go to 20075
20074 j=j+1
20075 if ((n20074-j).lt.0) go to 20076
      ibb(ibasis(j))=-1
      go to 20074
c
c     define the rest of ibasis(*)
20076 ip=mrelas
      j=1
      n20078=nvars+mrelas
      go to 20079
20078 j=j+1
20079 if ((n20078-j).lt.0) go to 20080
      if (.not.(ibb(j).gt.0)) go to 20082
      ip=ip+1
      ibasis(ip)=j
20082 go to 20078
20080 go to 20003
      end
