*deck spopt
      subroutine spopt (prgopt, mrelas, nvars, info, csc, ibasis, ropt,
     +   intopt, lopt)
c***begin prologue  spopt
c***subsidiary
c***purpose  subsidiary to splp
c***library   slatec
c***type      single precision (spopt-s, dpopt-d)
c***author  (unknown)
c***description
c
c     the editing required to convert this subroutine from single to
c     double precision involves the following character string changes.
c
c     use an editing command (change) /string-1/(to)string-2/.
c     /real (12 blanks)/double precision/,
c     /r1mach/d1mach/,/e0/d0/
c
c     revised 821122-1045
c     revised yymmdd-hhmm
c
c     this subroutine processes the option vector, prgopt(*),
c     and validates any modified data.
c
c***see also  splp
c***routines called  r1mach, xermsg
c***revision history  (yymmdd)
c   811215  date written
c   890531  changed all specific intrinsics to generic.  (wrb)
c   890605  removed unreferenced labels.  (wrb)
c   891214  prologue converted to version 4.0 format.  (bab)
c   900315  calls to xerror changed to calls to xermsg.  (thj)
c   900328  added type section.  (wrb)
c***end prologue  spopt
      real             abig,asmall,costsc,csc(*),eps,one,prgopt(*),
     * ropt(07),tolls,tune,zero,r1mach,tolabs
      integer ibasis(*),intopt(08)
      logical contin,usrbas,sizeup,savedt,colscp,cstscp,minprb,
     * stpedg,lopt(8)
c
c***first executable statement  spopt
      iopt=1
      zero=0.e0
      one=1.e0
      go to 30001
20002 continue
      go to 30002
c
20003 lopt(1)=contin
      lopt(2)=usrbas
      lopt(3)=sizeup
      lopt(4)=savedt
      lopt(5)=colscp
      lopt(6)=cstscp
      lopt(7)=minprb
      lopt(8)=stpedg
c
      intopt(1)=idg
      intopt(2)=ipagef
      intopt(3)=isave
      intopt(4)=mxitlp
      intopt(5)=kprint
      intopt(6)=itbrc
      intopt(7)=npp
      intopt(8)=lprg
c
      ropt(1)=eps
      ropt(2)=asmall
      ropt(3)=abig
      ropt(4)=costsc
      ropt(5)=tolls
      ropt(6)=tune
      ropt(7)=tolabs
      return
c
c
c     procedure (initialize parameters and process user options)
30001 contin = .false.
      usrbas = .false.
      sizeup = .false.
      savedt = .false.
      colscp = .false.
      cstscp = .false.
      minprb = .true.
      stpedg = .true.
c
c     get the machine rel. floating point accuracy value from the
c     library subprogram, r1mach( ).
      eps=r1mach(4)
      tolls=r1mach(4)
      tune=one
      tolabs=zero
c
c     define nominal file numbers for matrix pages and data saving.
      ipagef=1
      isave=2
      itbrc=10
      mxitlp=3*(nvars+mrelas)
      kprint=0
      idg=-4
      npp=nvars
      lprg=0
c
      last = 1
      iadbig=10000
      ictmax=1000
      ictopt= 0
20004 next=prgopt(last)
      if (.not.(next.le.0 .or. next.gt.iadbig)) go to 20006
c
c     the checks for small or large values of next are to prevent
c     working with undefined data.
      nerr=14
      call xermsg ('slatec', 'spopt',
     +   'in splp, the user option array has undefined data.', nerr,
     +   iopt)
      info=-nerr
      return
20006 if (.not.(next.eq.1)) go to 10001
      go to 20005
10001 if (.not.(ictopt.gt.ictmax)) go to 10002
      nerr=15
      call xermsg ('slatec', 'spopt',
     +   'in splp, option array processing is cycling.', nerr, iopt)
      info=-nerr
      return
10002 continue
      key = prgopt(last+1)
c
c     if key = 50, this is to be a maximization problem
c     instead of a minimization problem.
      if (.not.(key.eq.50)) go to 20010
      minprb = prgopt(last+2).eq.zero
      lds=3
      go to 20009
20010 continue
c
c     if key = 51, the level of output is being modified.
c     kprint = 0, no output
c            = 1, summary output
c            = 2, lots of output
c            = 3, even more output
      if (.not.(key.eq.51)) go to 20013
      kprint=prgopt(last+2)
      lds=3
      go to 20009
20013 continue
c
c     if key = 52, redefine the format and precision used
c     in the output.
      if (.not.(key.eq.52)) go to 20016
      if (prgopt(last+2).ne.zero) idg=prgopt(last+3)
      lds=4
      go to 20009
20016 continue
c
c     if key = 53, the allotted space for the sparse matrix
c     storage and/or sparse equation solving has been changed.
c     (processed in splp(). this is to compute the length of prgopt(*).)
      if (.not.(key.eq.53)) go to 20019
      lds=5
      go to 20009
20019 continue
c
c     if key = 54, redefine the file number where the pages
c     for the sparse matrix are stored.
      if (.not.(key.eq.54)) go to 20022
      if(prgopt(last+2).ne.zero) ipagef = prgopt(last+3)
      lds=4
      go to 20009
20022 continue
c
c     if key = 55,  a continuation for a problem may be requested.
      if (.not.(key .eq. 55)) go to 20025
      contin = prgopt(last+2).ne.zero
      lds=3
      go to 20009
20025 continue
c
c     if key = 56, redefine the file number where the saved data
c     will be stored.
      if (.not.(key.eq.56)) go to 20028
      if(prgopt(last+2).ne.zero) isave = prgopt(last+3)
      lds=4
      go to 20009
20028 continue
c
c     if key = 57, save data (on external file)  at mxitlp iterations or
c     the optimum, whichever comes first.
      if (.not.(key.eq.57)) go to 20031
      savedt=prgopt(last+2).ne.zero
      lds=3
      go to 20009
20031 continue
c
c     if key = 58,  see if problem is to run only a given
c     number of iterations.
      if (.not.(key.eq.58)) go to 20034
      if (prgopt(last+2).ne.zero) mxitlp = prgopt(last+3)
      lds=4
      go to 20009
20034 continue
c
c     if key = 59,  see if user provides the basis indices.
      if (.not.(key .eq. 59)) go to 20037
      usrbas = prgopt(last+2) .ne. zero
      if (.not.(usrbas)) go to 20040
      i=1
      n20043=mrelas
      go to 20044
20043 i=i+1
20044 if ((n20043-i).lt.0) go to 20045
      ibasis(i) = prgopt(last+2+i)
      go to 20043
20045 continue
20040 continue
      lds=mrelas+3
      go to 20009
20037 continue
c
c     if key = 60,  see if user has provided scaling of columns.
      if (.not.(key .eq. 60)) go to 20047
      colscp = prgopt(last+2).ne.zero
      if (.not.(colscp)) go to 20050
      j=1
      n20053=nvars
      go to 20054
20053 j=j+1
20054 if ((n20053-j).lt.0) go to 20055
      csc(j)=abs(prgopt(last+2+j))
      go to 20053
20055 continue
20050 continue
      lds=nvars+3
      go to 20009
20047 continue
c
c     if key = 61,  see if user has provided scaling of costs.
      if (.not.(key .eq. 61)) go to 20057
      cstscp = prgopt(last+2).ne.zero
      if (cstscp) costsc = prgopt(last+3)
      lds=4
      go to 20009
20057 continue
c
c     if key = 62,  see if size parameters are provided with the data.
c     these will be checked against the matrix element sizes later.
      if (.not.(key .eq. 62)) go to 20060
      sizeup = prgopt(last+2).ne.zero
      if (.not.(sizeup)) go to 20063
      asmall = prgopt(last+3)
      abig = prgopt(last+4)
20063 continue
      lds=5
      go to 20009
20060 continue
c
c     if key = 63, see if tolerance for linear system residual error is
c     provided.
      if (.not.(key .eq. 63)) go to 20066
      if (prgopt(last+2).ne.zero) tolls = max(eps,prgopt(last+3))
      lds=4
      go to 20009
20066 continue
c
c     if key = 64,  see if minimum reduced cost or steepest edge
c     descent is to be used for selecting variables to enter basis.
      if (.not.(key.eq.64)) go to 20069
      stpedg = prgopt(last+2).eq.zero
      lds=3
      go to 20009
20069 continue
c
c     if key = 65, set the number of iterations between recalculating
c     the error in the primal solution.
      if (.not.(key.eq.65)) go to 20072
      if (prgopt(last+2).ne.zero) itbrc=max(one,prgopt(last+3))
      lds=4
      go to 20009
20072 continue
c
c     if key = 66, set the number of negative reduced costs to be found
c     in the partial pricing strategy.
      if (.not.(key.eq.66)) go to 20075
      if (.not.(prgopt(last+2).ne.zero)) go to 20078
      npp=max(prgopt(last+3),one)
      npp=min(npp,nvars)
20078 continue
      lds=4
      go to 20009
20075 continue
c     if key = 67, change the tuning parameter to apply to the error
c     estimates for the primal and dual systems.
      if (.not.(key.eq.67)) go to 20081
      if (.not.(prgopt(last+2).ne.zero)) go to 20084
      tune=abs(prgopt(last+3))
20084 continue
      lds=4
      go to 20009
20081 continue
      if (.not.(key.eq.68)) go to 20087
      lds=6
      go to 20009
20087 continue
c
c     reset the absolute tolerance to be used on the feasibility
c     decision provided the relative error test failed.
      if (.not.(key.eq.69)) go to 20090
      if(prgopt(last+2).ne.zero)tolabs=prgopt(last+3)
      lds=4
      go to 20009
20090 continue
      continue
c
20009 ictopt = ictopt+1
      last = next
      lprg=lprg+lds
      go to 20004
20005 continue
      go to 20002
c
c     procedure (validate optionally modified data)
c
c     if user has defined the basis, check for validity of indices.
30002 if (.not.(usrbas)) go to 20093
      i=1
      n20096=mrelas
      go to 20097
20096 i=i+1
20097 if ((n20096-i).lt.0) go to 20098
      itest=ibasis(i)
      if (.not.(itest.le.0 .or.itest.gt.(nvars+mrelas))) go to 20100
      nerr=16
      call xermsg ('slatec', 'spopt',
     +   'in splp, an index of user-supplied basis is out of range.',
     +   nerr, iopt)
      info=-nerr
      return
20100 continue
      go to 20096
20098 continue
20093 continue
c
c     if user has provided size parameters, make sure they are ordered
c     and positive.
      if (.not.(sizeup)) go to 20103
      if (.not.(asmall.le.zero .or. abig.lt.asmall)) go to 20106
      nerr=17
      call xermsg ('slatec', 'spopt',
     +   'in splp, size parameters for matrix must be smallest and ' //
     +   'largest magnitudes of nonzero entries.', nerr, iopt)
      info=-nerr
      return
20106 continue
20103 continue
c
c     the number of iterations of rev. simplex steps must be positive.
      if (.not.(mxitlp.le.0)) go to 20109
      nerr=18
      call xermsg ('slatec', 'spopt',
     +   'in splp, the number of revised simplex steps between ' //
     +   'check-points must be positive.', nerr, iopt)
      info=-nerr
      return
20109 continue
c
c     check that save and page file numbers are defined and not equal.
      if (.not.(isave.le.0.or.ipagef.le.0.or.(isave.eq.ipagef))) go to 2
     *0112
      nerr=19
      call xermsg ('slatec', 'spopt',
     +   'in splp, file numbers for saved data and matrix pages ' //
     +   'must be positive and not equal.', nerr, iopt)
      info=-nerr
      return
20112 continue
      continue
      go to 20003
      end
