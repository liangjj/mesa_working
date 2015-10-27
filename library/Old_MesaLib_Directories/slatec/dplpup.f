*deck dplpup
      subroutine dplpup (dusrmt, mrelas, nvars, prgopt, dattrv, bl, bu,
     +   ind, info, amat, imat, sizeup, asmall, abig)
c***begin prologue  dplpup
c***subsidiary
c***purpose  subsidiary to dsplp
c***library   slatec
c***type      double precision (splpup-s, dplpup-d)
c***author  (unknown)
c***description
c
c     the editing required to convert this subroutine from single to
c     double precision involves the following character string changes.
c
c     use an editing command (change) /string-1/(to)string-2/.
c     /real (12 blanks)/double precision/.
c
c     revised 810613-1130
c     revised yymmdd-hhmm
c
c     this subroutine collects information about the bounds and matrix
c     from the user.  it is part of the dsplp( ) package.
c
c***see also  dsplp
c***routines called  dpchng, dpnnzr, xermsg
c***revision history  (yymmdd)
c   811215  date written
c   890531  changed all specific intrinsics to generic.  (wrb)
c   890605  corrected references to xerrwv.  (wrb)
c   890605  removed unreferenced labels.  (wrb)
c   891009  removed unreferenced variables.  (wrb)
c   891214  prologue converted to version 4.0 format.  (bab)
c   900315  calls to xerror changed to calls to xermsg.  (thj)
c   900328  added type section.  (wrb)
c   900510  convert xerrwv calls to xermsg calls, changed do-it-yourself
c           do loops to do loops.  (rwc)
c   900602  get rid of assigned gotos.  (rwc)
c***end prologue  dplpup
      double precision abig,aij,amat(*),amn,amx,asmall,bl(*),
     * bu(*),dattrv(*),prgopt(*),xval,zero
      integer iflag(10),imat(*),ind(*)
      logical sizeup,first
      character*8 xern1, xern2
      character*16 xern3, xern4
c
c***first executable statement  dplpup
      zero = 0.d0
c
c     check user-supplied bounds
c
c     check that ind(*) values are 1,2,3 or 4.
c     also check consistency of upper and lower bounds.
c
      do 10 j=1,nvars
         if (ind(j).lt.1 .or. ind(j).gt.4) then
            write (xern1, '(i8)') j
            call xermsg ('slatec', 'dplpup',
     *         'in dsplp, independent variable = ' // xern1 //
     *         ' is not defined.', 10, 1)
            info = -10
            return
         endif
c
         if (ind(j).eq.3) then
            if (bl(j).gt.bu(j)) then
               write (xern1, '(i8)') j
               write (xern3, '(1pe15.6)') bl(j)
               write (xern4, '(1pe15.6)') bu(j)
               call xermsg ('slatec', 'dplpup',
     *            'in dsplp, lower bound = ' // xern3 //
     *            ' and upper bound = ' // xern4 //
     *            ' for independent variable = ' // xern1 //
     *            ' are not consistent.', 11, 1)
               return
            endif
         endif
   10 continue
c
      do 20 i=nvars+1,nvars+mrelas
         if (ind(i).lt.1 .or. ind(i).gt.4) then
            write (xern1, '(i8)') i-nvars
            call xermsg ('slatec', 'dplpup',
     *         'in dsplp, dependent variable = ' // xern1 //
     *         ' is not defined.', 12, 1)
            info = -12
            return
         endif
c
         if (ind(i).eq.3) then
            if (bl(i).gt.bu(i)) then
               write (xern1, '(i8)') i
               write (xern3, '(1pe15.6)') bl(i)
               write (xern4, '(1pe15.6)') bu(i)
               call xermsg ('slatec', 'dplpup',
     *            'in dsplp, lower bound = ' // xern3 //
     *            ' and upper bound = ' // xern4 //
     *            ' for dependant variable = ' // xern1 //
     *            ' are not consistent.',13,1)
               info = -13
               return
            endif
         endif
   20 continue
c
c     get updates or data for matrix from the user
c
c     get the elements of the matrix from the user.  it will be stored
c     by columns using the sparse storage codes of rj hanson and
c     ja wisniewski.
c
      iflag(1) = 1
c
c     keep accepting elements until the user is finished giving them.
c     limit this loop to 2*nvars*mrelas iterations.
c
      itmax = 2*nvars*mrelas+1
      itcnt = 0
      first = .true.
c
c     check on the iteration count.
c
   30 itcnt = itcnt+1
      if (itcnt.gt.itmax) then
         call xermsg ('slatec', 'dplpup',
     +      'in dsplp, more than 2*nvars*mrelas iterations defining ' //
     +      'or updating matrix data.', 7, 1)
         info = -7
         return
      endif
c
      aij = zero
      call dusrmt(i,j,aij,indcat,prgopt,dattrv,iflag)
      if (iflag(1).eq.1) then
         iflag(1) = 2
         go to 30
      endif
c
c     check to see that the subscripts i and j are valid.
c
      if (i.lt.1 .or. i.gt.mrelas .or. j.lt.1 .or. j.gt.nvars) then
c
c        check on size of matrix data
c        record the largest and smallest(in magnitude) nonzero elements.
c
         if (iflag(1).eq.3) then
            if (sizeup .and. abs(aij).ne.zero) then
               if (first) then
                  amx = abs(aij)
                  amn = abs(aij)
                  first = .false.
               elseif (abs(aij).gt.amx) then
                  amx = abs(aij)
               elseif (abs(aij).lt.amn) then
                  amn = abs(aij)
               endif
            endif
            go to 40
         endif
c
         write (xern1, '(i8)') i
         write (xern2, '(i8)') j
         call xermsg ('slatec', 'dplpup',
     *      'in dsplp, row index = ' // xern1 // ' or column index = '
     *      // xern2 // ' is out of range.', 8, 1)
         info = -8
         return
      endif
c
c     if indcat=0 then set a(i,j)=aij.
c     if indcat=1 then accumulate element, a(i,j)=a(i,j)+aij.
c
      if (indcat.eq.0) then
         call dpchng(i,aij,iplace,amat,imat,j)
      elseif (indcat.eq.1) then
         index = -(i-1)
         call dpnnzr(index,xval,iplace,amat,imat,j)
         if (index.eq.i) aij=aij+xval
         call dpchng(i,aij,iplace,amat,imat,j)
      else
         write (xern1, '(i8)') indcat
         call xermsg ('slatec', 'dplpup',
     *      'in dsplp, indication flag = ' // xern1 //
     *      ' for matrix data must be either 0 or 1.', 9, 1)
         info = -9
         return
      endif
c
c     check on size of matrix data
c     record the largest and smallest(in magnitude) nonzero elements.
c
      if (sizeup .and. abs(aij).ne.zero) then
         if (first) then
            amx = abs(aij)
            amn = abs(aij)
            first = .false.
         elseif (abs(aij).gt.amx) then
            amx = abs(aij)
         elseif (abs(aij).lt.amn) then
            amn = abs(aij)
         endif
      endif
      if (iflag(1).ne.3) go to 30
c
   40 if (sizeup .and. .not. first) then
         if (amn.lt.asmall .or. amx.gt.abig) then
            call xermsg ('slatec', 'dplpup',
     +         'in dsplp, a matrix element''s size is out of the ' //
     +         'specified range.', 22, 1)
            info = -22
            return
         endif
      endif
      return
      end
