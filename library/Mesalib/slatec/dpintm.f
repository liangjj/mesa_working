*deck dpintm
      subroutine dpintm (m, n, sx, ix, lmx, ipagef)
c***begin prologue  dpintm
c***subsidiary
c***purpose  subsidiary to dsplp
c***library   slatec
c***type      double precision (pinitm-s, dpintm-d)
c***author  hanson, r. j., (snla)
c           wisniewski, j. a., (snla)
c***description
c
c     dpintm limits the type of storage to a sequential scheme.
c     the matrix is stored by columns.
c     sparse matrix initialization subroutine.
c
c            m=number of rows of the matrix.
c            n=number of columns of the matrix.
c  sx(*),ix(*)=the work arrays which are used to store the sparse
c              matrix.  these arrays are automatically maintained by
c              the package for the user.
c          lmx=length of the work array sx(*).
c              lmx must be at least n+7 where
c              for greatest efficiency lmx should be at least n+nz+6
c              where nz is the maximum number of nonzeroes to be
c              stored in the matrix.  values of lmx between n+7 and
c              n+nz+6 will cause demand paging to occur.
c              this is implemented by the package.
c              ix(*) must be dimensioned at least lmx
c      ipagef=unit number where demand pages will be stored.
c
c     this subroutine is a modification of the subroutine linitm,
c     sandia labs. rept. sand78-0785.
c     modifications by k.l. hiebert and r.j. hanson
c     revised 811130-1000
c     revised yymmdd-hhmm
c
c***see also  dsplp
c***routines called  xermsg
c***revision history  (yymmdd)
c   811215  date written
c   890831  modified array declarations.  (wrb)
c   891214  prologue converted to version 4.0 format.  (bab)
c   900315  calls to xerror changed to calls to xermsg.  (thj)
c   900328  added type section.  (wrb)
c   910403  updated author and description sections.  (wrb)
c***end prologue  dpintm
      double precision sx(*),zero,one
      dimension ix(*)
      save zero, one
      data zero,one /0.d0,1.d0/
c***first executable statement  dpintm
      iopt=1
c
c     check for input errors.
c
      if (.not.(m.le.0 .or. n.le.0)) go to 20002
      nerr=55
      call xermsg ('slatec', 'dpintm',
     +   'matrix dimension m or n .le. 0', nerr, iopt)
c
c     verify if value of lmx is large enough.
c
20002 if (.not.(lmx.lt.n+7)) go to 20005
      nerr=55
      call xermsg ('slatec', 'dpintm',
     +   'the value of lmx is too small', nerr, iopt)
c
c     initialize data structure independent values.
c
20005 sx(1)=zero
      sx(2)=zero
      sx(3)=ipagef
      ix(1)=lmx
      ix(2)=m
      ix(3)=n
      ix(4)=0
      sx(lmx-1)=zero
      sx(lmx)=-one
      ix(lmx-1)=-1
      lp4=n+4
c
c     initialize data structure dependent values.
c
      i=4
      n20008=lp4
      go to 20009
20008 i=i+1
20009 if ((n20008-i).lt.0) go to 20010
      sx(i)=zero
      go to 20008
20010 i=5
      n20012=lp4
      go to 20013
20012 i=i+1
20013 if ((n20012-i).lt.0) go to 20014
      ix(i)=lp4
      go to 20012
20014 sx(n+5)=zero
      ix(n+5)=0
      ix(lmx)=0
c
c     initialization complete.
c
      return
      end
