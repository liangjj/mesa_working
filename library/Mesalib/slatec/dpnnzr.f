*deck dpnnzr
      subroutine dpnnzr (i, xval, iplace, sx, ix, ircx)
c***begin prologue  dpnnzr
c***subsidiary
c***purpose  subsidiary to dsplp
c***library   slatec
c***type      double precision (pnnzrs-s, dpnnzr-d)
c***author  hanson, r. j., (snla)
c           wisniewski, j. a., (snla)
c***description
c
c     dpnnzr limits the type of storage to a sequential scheme.
c     sparse matrix non zero retrieval subroutine.
c
c     subroutine dpnnzr() gets the next nonzero value in row or column
c     +/- ircx with an index greater than the value of i.
c
c             i absolute value of this subscript is to be exceeded
c               in the search for the next nonzero value. a negative
c               or zero value of i causes the search to start at
c               the beginning of the vector.  a positive value
c               of i causes the search to continue from the last place
c               accessed. on output, the argument i
c               contains the value of the subscript found.  an output
c               value of i equal to zero indicates that all components
c               with an index greater than the input value of i are
c               zero.
c          xval value of the nonzero element found.  on output,
c               xval=0. whenever i=0.
c     iplace pointer information which is maintained by the package.
c   sx(*),ix(*) the work arrays which are used to store the sparse
c               matrix.  these array contents are automatically
c               maintained by the package for the user.
c          ircx points to the vector of the matrix being scanned.  a
c               negative value of ircx indicates that row -ircx is to be
c               scanned.  a positive value of ircx indicates that
c               column ircx is to be scanned.  a zero value of ircx is
c               an error.
c
c     this subroutine is a modification of the subroutine lnnzrs,
c     sandia labs. rept. sand78-0785.
c     modifications by k.l. hiebert and r.j. hanson
c     revised 811130-1000
c     revised yymmdd-hhmm
c
c***see also  dsplp
c***routines called  idloc, xermsg
c***revision history  (yymmdd)
c   811215  date written
c   890531  changed all specific intrinsics to generic.  (wrb)
c   890605  removed unreferenced labels.  (wrb)
c   890606  changed references from iploc to idloc.  (wrb)
c   891214  prologue converted to version 4.0 format.  (bab)
c   900315  calls to xerror changed to calls to xermsg.  (thj)
c   900328  added type section.  (wrb)
c   910403  updated author and description sections.  (wrb)
c***end prologue  dpnnzr
      dimension ix(*)
      double precision xval,sx(*),zero
      save zero
      data zero /0.d0/
c***first executable statement  dpnnzr
      iopt=1
c
c     check validity of row/col. index.
c
      if (.not.(ircx .eq.0)) go to 20002
      nerr=55
      call xermsg ('slatec', 'dpnnzr', 'ircx=0', nerr, iopt)
c
c     lmx is the length of the in-memory storage area.
c
20002 lmx = ix(1)
      if (.not.(ircx.lt.0)) go to 20005
c
c     check subscripts of the row. the row number must be .le. m and
c     the index must be .le. n.
c
      if (.not.(ix(2).lt.-ircx .or. ix(3).lt.abs(i))) go to 20008
      nerr=55
      call xermsg ('slatec', 'dpnnzr',
     +   'subscripts for array element to be accessed were out of ' //
     +   'bounds.', nerr, iopt)
20008 l=ix(3)
      go to 20006
c
c     check subscripts of the column. the col. number must be .le. n and
c     the index must be .le. m.
c
20005 if (.not.(ircx.gt.ix(3) .or. abs(i).gt.ix(2))) go to 20011
      nerr=55
      call xermsg ('slatec', 'dpnnzr',
     +   'subscripts for array element to be accessed were out of ' //
     +   'bounds', nerr, iopt)
20011 l=ix(2)
c
c     here l is the largest possible subscript within the vector.
c
20006 j=abs(ircx)
      ll=ix(3)+4
      lpg = lmx - ll
      if (.not.(ircx.gt.0)) go to 20014
c
c     searching for the next nonzero in a column.
c
c     initialize starting locations..
      if (.not.(i.le.0)) go to 20017
      if (.not.(j.eq.1)) go to 20020
      iplace=ll+1
      go to 20021
20020 iplace=ix(j+3)+1
20021 continue
c
c     the case i.le.0 signals that the scan for the entry
c     is to begin at the start of the vector.
c
20017 i = abs(i)
      if (.not.(j.eq.1)) go to 20023
      istart = ll+1
      go to 20024
20023 istart=ix(j+3)+1
20024 iend = ix(j+4)
c
c     validate iplace. set to start of vector if out of range.
c
      if (.not.(istart.gt.iplace .or. iplace.gt.iend)) go to 20026
      if (.not.(j.eq.1)) go to 20029
      iplace=ll+1
      go to 20030
20029 iplace=ix(j+3)+1
20030 continue
c
c     scan through several pages, if necessary, to find matrix entry.
c
20026 ipl = idloc(iplace,sx,ix)
c
c     fix up iplace and ipl if they point to paging data.
c     this is necessary because there is control information at the
c     end of each page.
c
      idiff = lmx - ipl
      if (.not.(idiff.le.1.and.ix(lmx-1).gt.0)) go to 20032
c
c     update the relative address in a new page.
c
      iplace = iplace + idiff + 1
      ipl = idloc(iplace,sx,ix)
20032 np = abs(ix(lmx-1))
      go to 20036
20035 if (ilast.eq.iend) go to 20037
20036 ilast = min(iend,np*lpg+ll-2)
c
c     the virtual end of the data for this page is ilast.
c
      il = idloc(ilast,sx,ix)
      il = min(il,lmx-2)
c
c     the relative end of data for this page is il.
c     search for a nonzero value with an index .gt. i on the present
c     page.
c
20038 if (.not.(.not.(ipl.ge.il.or.(ix(ipl).gt.i.and.sx(ipl).ne.zero))))
     * go to 20039
      ipl=ipl+1
      go to 20038
c
c     test if we have found the next nonzero.
c
20039 if (.not.(ix(ipl).gt.i .and. sx(ipl).ne.zero .and. ipl.le.il)) go
     *to 20040
      i = ix(ipl)
      xval = sx(ipl)
      iplace = (np-1)*lpg + ipl
      return
c
c     update to scan the next page.
20040 ipl = ll + 1
      np = np + 1
      go to 20035
c
c     no data was found. end of vector encountered.
c
20037 i = 0
      xval = zero
      il = il + 1
      if(il.eq.lmx-1) il = il + 2
c
c     if a new item would be inserted, iplace points to the place
c     to put it.
c
      iplace = (np-1)*lpg + il
      return
c
c     search a row for the next nonzero.
c     find element j=abs(ircx) in rows abs(i)+1,...,l.
c
20014 i=abs(i)
c
c     check for end of vector.
c
      if (.not.(i.eq.l)) go to 20043
      i=0
      xval=zero
      return
20043 i1 = i+1
      ii=i1
      n20046=l
      go to 20047
20046 ii=ii+1
20047 if ((n20046-ii).lt.0) go to 20048
c
c     initialize ipploc for orthogonal scan.
c     look for j as a subscript in rows ii, ii=i+1,...,l.
c
      if (.not.(ii.eq.1)) go to 20050
      ipploc = ll + 1
      go to 20051
20050 ipploc = ix(ii+3) + 1
20051 iend = ix(ii+4)
c
c     scan through several pages, if necessary, to find matrix entry.
c
      ipl = idloc(ipploc,sx,ix)
c
c     fix up ipploc and ipl to point to matrix data.
c
      idiff = lmx - ipl
      if (.not.(idiff.le.1.and.ix(lmx-1).gt.0)) go to 20053
      ipploc = ipploc + idiff + 1
      ipl = idloc(ipploc,sx,ix)
20053 np = abs(ix(lmx-1))
      go to 20057
20056 if (ilast.eq.iend) go to 20058
20057 ilast = min(iend,np*lpg+ll-2)
      il = idloc(ilast,sx,ix)
      il = min(il,lmx-2)
20059 if (.not.(.not.(ipl.ge.il .or. ix(ipl).ge.j))) go to 20060
      ipl=ipl+1
      go to 20059
c
c     test if we have found the next nonzero.
c
20060 if (.not.(ix(ipl).eq.j .and. sx(ipl).ne.zero .and. ipl.le.il)) go
     *to 20061
      i = ii
      xval = sx(ipl)
      return
20061 if(ix(ipl).ge.j) ilast = iend
      ipl = ll + 1
      np = np + 1
      go to 20056
20058 go to 20046
c
c     orthogonal scan failed. the value j was not a subscript
c     in any row.
c
20048 i=0
      xval=zero
      return
      end
