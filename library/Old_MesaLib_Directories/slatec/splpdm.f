*deck splpdm
      subroutine splpdm (mrelas, nvars, lmx, lbm, nredc, info, iopt,
     +   ibasis, imat, ibrc, ipr, iwr, ind, ibb, anorm, eps, uu, gg,
     +   amat, basmat, csc, wr, singlr, redbas)
c***begin prologue  splpdm
c***subsidiary
c***purpose  subsidiary to splp
c***library   slatec
c***type      single precision (splpdm-s, dplpdm-d)
c***author  (unknown)
c***description
c
c     this subprogram is from the splp( ) package.  it performs the
c     task of defining the entries of the basis matrix and
c     decomposing it using the la05 package.
c     it is the main part of the procedure (decompose basis matrix).
c
c***see also  splp
c***routines called  la05as, pnnzrs, sasum, xermsg
c***common blocks    la05ds
c***revision history  (yymmdd)
c   811215  date written
c   890605  corrected references to xerrwv.  (wrb)
c   890605  removed unreferenced labels.  (wrb)
c   891214  prologue converted to version 4.0 format.  (bab)
c   900315  calls to xerror changed to calls to xermsg.  (thj)
c   900328  added type section.  (wrb)
c   900510  convert xerrwv calls to xermsg calls, changed do-it-yourself
c           do loops to do loops.  (rwc)
c***end prologue  splpdm
      integer ibasis(*),imat(*),ibrc(lbm,2),ipr(*),iwr(*),ind(*),ibb(*)
      real             amat(*),basmat(*),csc(*),wr(*),anorm,eps,gg,
     * one,small,uu,zero
      logical singlr,redbas
      character*16 xern3
c
c     common block used by la05 () package..
      common /la05ds/ small,lp,lenl,lenu,ncp,lrow,lcol
c
c***first executable statement  splpdm
      zero = 0.e0
      one = 1.e0
c
c     define basis matrix by columns for sparse matrix equation solver.
c     the la05as() subprogram requires the nonzero entries of the matrix
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
   10       call pnnzrs(i,aij,iplace,amat,imat,j)
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
      anorm = sasum(nzbm,basmat,1)
      small = eps*anorm
c
c     get an l-u factorization of the basis matrix.
c
      nredc = nredc+1
      redbas = .true.
      call la05as(basmat,ibrc,nzbm,lbm,mrelas,ipr,iwr,wr,gg,uu)
c
c     check return value of error flag, gg.
c
      if (gg.ge.zero) return
      if (gg.eq.(-7.)) then
         call xermsg ('slatec', 'splpdm',
     *      'in splp, short on storage for la05as.  ' //
     *      'use prgopt(*) to give more.', 28, iopt)
         info = -28
      elseif (gg.eq.(-5.)) then
         singlr = .true.
      else
         write (xern3, '(1pe15.6)') gg
         call xermsg ('slatec', 'splpdm',
     *      'in splp, la05as returned error flag = ' // xern3,
     *      27, iopt)
         info = -27
      endif
      return
      end
