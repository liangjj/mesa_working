*deck dplpdm
      subroutine dplpdm (mrelas, nvars, lmx, lbm, nredc, info, iopt,
     +   ibasis, imat, ibrc, ipr, iwr, ind, ibb, anorm, eps, uu, gg,
     +   amat, basmat, csc, wr, singlr, redbas)
c***begin prologue  dplpdm
c***subsidiary
c***purpose  subsidiary to dsplp
c***library   slatec
c***type      double precision (splpdm-s, dplpdm-d)
c***author  (unknown)
c***description
c
c     this subprogram is from the dsplp( ) package.  it performs the
c     task of defining the entries of the basis matrix and
c     decomposing it using the la05 package.
c     it is the main part of the procedure (decompose basis matrix).
c
c***see also  dsplp
c***routines called  dasum, dpnnzr, la05ad, xermsg
c***common blocks    la05dd
c***revision history  (yymmdd)
c   811215  date written
c   890605  added dasum to list of double precision variables.
c   890605  removed unreferenced labels.  (wrb)
c   891009  removed unreferenced variable.  (wrb)
c   891214  prologue converted to version 4.0 format.  (bab)
c   900315  calls to xerror changed to calls to xermsg.  (thj)
c   900328  added type section.  (wrb)
c   900510  convert xerrwv calls to xermsg calls, convert do-it-yourself
c           do loops to do loops.  (rwc)
c***end prologue  dplpdm
      integer ibasis(*),imat(*),ibrc(lbm,2),ipr(*),iwr(*),ind(*),ibb(*)
      double precision aij,amat(*),basmat(*),csc(*),wr(*),anorm,dasum,
     * eps,gg,one,small,uu,zero
      logical singlr,redbas
      character*16 xern3
c
c     common block used by la05 () package..
      common /la05dd/ small,lp,lenl,lenu,ncp,lrow,lcol
c
c***first executable statement  dplpdm
      zero = 0.d0
      one = 1.d0
c
c     define basis matrix by columns for sparse matrix equation solver.
c     the la05ad() subprogram requires the nonzero entries of the matrix
c     together with the row and column indices.
c
      nzbm = 0
c
c     define dependent variable columns. these are
c     cols. of the identity matrix and implicitly generated.
c
      do 20 k = 1,mrelas
         j = ibasis(k)
         if (j.gt.nvars) then
            nzbm = nzbm+1
            if (ind(j).eq.2) then
               basmat(nzbm) = one
            else
               basmat(nzbm) = -one
            endif
            ibrc(nzbm,1) = j-nvars
            ibrc(nzbm,2) = k
         else
c
c           define the indep. variable cols.  this requires retrieving
c           the cols. from the sparse matrix data structure.
c
            i = 0
   10       call dpnnzr(i,aij,iplace,amat,imat,j)
            if (i.gt.0) then
               nzbm = nzbm+1
               basmat(nzbm) = aij*csc(j)
               ibrc(nzbm,1) = i
               ibrc(nzbm,2) = k
               go to 10
            endif
         endif
   20 continue
c
      singlr = .false.
c
c     recompute matrix norm using crude norm  =  sum of magnitudes.
c
      anorm = dasum(nzbm,basmat,1)
      small = eps*anorm
c
c     get an l-u factorization of the basis matrix.
c
      nredc = nredc+1
      redbas = .true.
      call la05ad(basmat,ibrc,nzbm,lbm,mrelas,ipr,iwr,wr,gg,uu)
c
c     check return value of error flag, gg.
c
      if (gg.ge.zero) return
      if (gg.eq.(-7.)) then
         call xermsg ('slatec', 'dplpdm',
     *      'in dsplp, short on storage for la05ad.  ' //
     *      'use prgopt(*) to give more.', 28, iopt)
         info = -28
      elseif (gg.eq.(-5.)) then
         singlr = .true.
      else
         write (xern3, '(1pe15.6)') gg
         call xermsg ('slatec', 'dplpdm',
     *      'in dsplp, la05ad returned error flag = ' // xern3,
     *      27, iopt)
         info = -27
      endif
      return
      end
