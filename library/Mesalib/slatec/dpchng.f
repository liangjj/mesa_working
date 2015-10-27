*deck dpchng
      subroutine dpchng (ii, xval, iplace, sx, ix, ircx)
c***begin prologue  dpchng
c***subsidiary
c***purpose  subsidiary to dsplp
c***library   slatec
c***type      double precision (pchngs-s, dpchng-d)
c***author  hanson, r. j., (snla)
c           wisniewski, j. a., (snla)
c***description
c
c     subroutine dpchng changes element ii in vector +/- ircx to the
c     value xval.
c     dpchng limits the type of storage to a sequential scheme.
c     sparse matrix element alteration subroutine.
c
c            ii the absolute value of this integer is the subscript for
c               the element to be changed.
c          xval new value of the matrix element being changed.
c     iplace pointer information which is maintained by the package.
c   sx(*),ix(*) the work arrays which are used to store the sparse
c               matrix. these arrays are automatically maintained by the
c               package for the user.
c          ircx points to the vector of the matrix being updated.
c               a negative value of ircx indicates that row -ircx is
c               being updated.  a positive value of ircx indicates that
c               column ircx is being updated.  a zero value of ircx is
c               an error.
c
c     since data items are kept sorted in the sequential data structure,
c     changing a matrix element can require the movement of all the data
c     items in the matrix. for this reason, it is suggested that data
c     items be added a col. at a time, in ascending col. sequence.
c     furthermore, since deleting items from the data structure may also
c     require moving large amounts of data, zero elements are explicitly
c     stored in the matrix.
c
c     this subroutine is a modification of the subroutine lchngs,
c     sandia labs. rept. sand78-0785.
c     modifications by k.l. hiebert and r.j. hanson
c     revised 811130-1000
c     revised yymmdd-hhmm
c
c***see also  dsplp
c***routines called  dprwpg, idloc, xermsg
c***revision history  (yymmdd)
c   811215  date written
c   890531  changed all specific intrinsics to generic.  (wrb)
c   890606  changed references from iploc to idloc.  (wrb)
c   891214  prologue converted to version 4.0 format.  (bab)
c   900315  calls to xerror changed to calls to xermsg.  (thj)
c   900328  added type section.  (wrb)
c   910403  updated author and description sections.  (wrb)
c***end prologue  dpchng
      dimension ix(*)
      integer idloc
      double precision sx(*),xval,zero,one,sxlast,sxval
      save zero, one
      data zero,one /0.d0,1.d0/
c***first executable statement  dpchng
      iopt=1
c
c     determine null-cases..
      if(ii.eq.0) return
c
c     check validity of row/col. index.
c
      if (.not.(ircx.eq.0)) go to 20002
      nerr=55
      call xermsg ('slatec', 'dpchng', 'ircx=0', nerr, iopt)
20002 lmx = ix(1)
c
c     lmx is the length of the in-memory storage area.
c
      if (.not.(ircx.lt.0)) go to 20005
c
c     check subscripts of the row. the row number must be .le. m and
c     the index must be .le. n.
c
      if (.not.(ix(2).lt.-ircx .or. ix(3).lt.abs(ii))) go to 20008
      nerr=55
      call xermsg ('slatec', 'dpchng',
     +   'subscripts for array element to be accessed were out of ' //
     +   'bounds', nerr, iopt)
20008 go to 20006
c
c     check subscripts of the column. the col. number must be .le. n and
c     the index must be .le. m.
c
20005 if (.not.(ix(3).lt.ircx .or. ix(2).lt.abs(ii))) go to 20011
      nerr=55
      call xermsg ('slatec', 'dpchng',
     +   'subscripts for array element to be accessed were out of ' //
     +   'bounds', nerr, iopt)
20011 continue
c
c     set i to be the element of row/column j to be changed.
c
20006 if (.not.(ircx.gt.0)) go to 20014
      i = abs(ii)
      j = abs(ircx)
      go to 20015
20014 i = abs(ircx)
      j = abs(ii)
c
c     the integer ll points to the start of the matrix element data.
c
20015 ll=ix(3)+4
      ii = abs(ii)
      lpg = lmx - ll
c
c     set iplace to start our scan for the element at the beginning
c     of the vector.
c
      if (.not.(j.eq.1)) go to 20017
      iplace=ll+1
      go to 20018
20017 iplace=ix(j+3)+1
c
c     iend points to the last element of the vector to be scanned.
c
20018 iend = ix(j+4)
c
c     scan through several pages, if necessary, to find matrix element.
c
      ipl = idloc(iplace,sx,ix)
      np = abs(ix(lmx-1))
      go to 20021
20020 if (ilast.eq.iend) go to 20022
c
c     the virtual end of data for this page is ilast.
c
20021 ilast = min(iend,np*lpg+ll-2)
c
c     the relative end of data for this page is il.
c     search for a matrix value with an index .ge. i on the present
c     page.
c
      il = idloc(ilast,sx,ix)
      il = min(il,lmx-2)
20023 if (.not.(.not.(ipl.ge.il .or. ix(ipl).ge.i))) go to 20024
      ipl=ipl+1
      go to 20023
c
c     set iplace and store data item if found.
c
20024 if (.not.(ix(ipl).eq.i .and. ipl.le.il)) go to 20025
      sx(ipl) = xval
      sx(lmx) = one
      return
c
c     exit from loop if item was found.
c
20025 if(ix(ipl).gt.i .and. ipl.le.il) ilast = iend
      if (.not.(ilast.ne.iend)) go to 20028
      ipl = ll + 1
      np = np + 1
20028 go to 20020
c
c     insert new data item into location at iplace(ipl).
c
20022 if (.not.(ipl.gt.il.or.(ipl.eq.il.and.i.gt.ix(ipl)))) go to 20031
      ipl = il + 1
      if(ipl.eq.lmx-1) ipl = ipl + 2
20031 iplace = (np-1)*lpg + ipl
c
c     go to a new page, if necessary, to insert the item.
c
      if (.not.(ipl.le.lmx .or. ix(lmx-1).ge.0)) go to 20034
      ipl=idloc(iplace,sx,ix)
20034 iend = ix(ll)
      np = abs(ix(lmx-1))
      sxval = xval
c
c     loop through all subsequent pages of the matrix moving data down.
c     this is necessary to make room for the new matrix element and
c     keep the entries sorted.
c
      go to 20038
20037 if (ix(lmx-1).le.0) go to 20039
20038 ilast = min(iend,np*lpg+ll-2)
      il = idloc(ilast,sx,ix)
      il = min(il,lmx-2)
      sxlast = sx(il)
      ixlast = ix(il)
      istart = ipl + 1
      if (.not.(istart.le.il)) go to 20040
      k = istart + il
      do 50 jj=istart,il
      sx(k-jj) = sx(k-jj-1)
      ix(k-jj) = ix(k-jj-1)
50    continue
      sx(lmx) = one
20040 if (.not.(ipl.le.lmx)) go to 20043
      sx(ipl) = sxval
      ix(ipl) = i
      sxval = sxlast
      i = ixlast
      sx(lmx) = one
      if (.not.(ix(lmx-1).gt.0)) go to 20046
      ipl = ll + 1
      np = np + 1
20046 continue
20043 go to 20037
20039 np = abs(ix(lmx-1))
c
c     determine if a new page is to be created for the last element
c     moved down.
c
      il = il + 1
      if (.not.(il.eq.lmx-1)) go to 20049
c
c     create a new page.
c
      ix(lmx-1) = np
c
c     write the old page.
c
      sx(lmx) = zero
      key = 2
      call dprwpg(key,np,lpg,sx,ix)
      sx(lmx) = one
c
c     store last element moved down in a new page.
c
      ipl = ll + 1
      np = np + 1
      ix(lmx-1) = -np
      sx(ipl) = sxval
      ix(ipl) = i
      go to 20050
c
c     last element moved remained on the old page.
c
20049 if (.not.(ipl.ne.il)) go to 20052
      sx(il) = sxval
      ix(il) = i
      sx(lmx) = one
20052 continue
c
c     increment pointers to last element in vectors j,j+1,... .
c
20050 jstart = j + 4
      jj=jstart
      n20055=ll
      go to 20056
20055 jj=jj+1
20056 if ((n20055-jj).lt.0) go to 20057
      ix(jj) = ix(jj) + 1
      if(mod(ix(jj)-ll,lpg).eq.lpg-1) ix(jj) = ix(jj) + 2
      go to 20055
c
c     iplace points to the inserted data item.
c
20057 ipl=idloc(iplace,sx,ix)
      return
      end
