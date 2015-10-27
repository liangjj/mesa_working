*deck dprwvr
      subroutine dprwvr (key, ipage, lpg, sx, ix)
c***begin prologue  dprwvr
c***subsidiary
c***purpose  subsidiary to dsplp
c***library   slatec
c***type      double precision (prwvir-s, dprwvr-d)
c***author  hanson, r. j., (snla)
c           wisniewski, j. a., (snla)
c***description
c
c     dprwvr limits the type of storage to a sequential sparse matrix
c     storage scheme.  the page storage is on random access disk.
c     dprwvr is part of the sparse lp package, dsplp.
c
c     key       is a flag which indicates whether a read or write
c               operation is to be performed. a value of key=1 indicates
c               a read. a value of key=2 indicates a write.
c     ipage     is the page of matrix mn we are accessing.
c     lpg       is the length of the page.
c   sx(*),ix(*) is the matrix data.
c
c     this subroutine is a modification of the subroutine lrwvir,
c     sandia labs. rept. sand78-0785.
c     modifications by k.l. hiebert and r.j. hanson
c
c***see also  dsplp
c***routines called  dreadp, dwritp, sopenm
c***revision history  (yymmdd)
c   811215  date written
c   891009  removed unreferenced variables.  (wrb)
c   891214  prologue converted to version 4.0 format.  (bab)
c   900328  added type section.  (wrb)
c   910403  updated author and description sections.  (wrb)
c***end prologue  dprwvr
      dimension ix(*)
      double precision sx(*),zero,one
      logical first
      save zero, one
      data zero,one/0.d0,1.d0/
c***first executable statement  dprwvr
c
c     compute starting address of page.
c
      ipagef=sx(3)
      istart = ix(3) + 5
c
c     open random access file number ipagef, if first page write.
c
      first=sx(4).eq.zero
      if (.not.(first)) go to 20002
      call sopenm(ipagef,lpg)
      sx(4)=one
c
c     perform either a read or a write.
c
20002 iaddr = 2*ipage - 1
      if (.not.(key.eq.1)) go to 20005
      call dreadp(ipagef,ix(istart),sx(istart),lpg,iaddr)
      go to 20006
20005 if (.not.(key.eq.2)) go to 10001
      call dwritp(ipagef,ix(istart),sx(istart),lpg,iaddr)
10001 continue
20006 return
      end
