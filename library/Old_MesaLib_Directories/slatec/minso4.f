*deck minso4
      subroutine minso4 (usol, idmn, zn, zm, pertb)
c***begin prologue  minso4
c***subsidiary
c***purpose  subsidiary to sepx4
c***library   slatec
c***type      single precision (minso4-s)
c***author  (unknown)
c***description
c
c     this subroutine orthogonalizes the array usol with respect to
c     the constant array in a weighted least squares norm.
c
c     entry at minso4 occurs when the final solution is
c     to be minimized with respect to the weighted
c     least squares norm.
c
c***see also  sepx4
c***routines called  (none)
c***common blocks    spl4
c***revision history  (yymmdd)
c   801001  date written
c   891214  prologue converted to version 4.0 format.  (bab)
c   900402  added type section.  (wrb)
c***end prologue  minso4
c
      common /spl4/   kswx       ,kswy       ,k          ,l          ,
     1                ait        ,bit        ,cit        ,dit        ,
     2                mit        ,nit        ,is         ,ms         ,
     3                js         ,ns         ,dlx        ,dly        ,
     4                tdlx3      ,tdly3      ,dlx4       ,dly4
      dimension       usol(idmn,*)           ,zn(*)      ,zm(*)
c***first executable statement  minso4
      istr = 1
      ifnl = k
      jstr = 1
      jfnl = l
c
c     compute weighted inner products
c
      ute = 0.0
      ete = 0.0
      do  20 i=is,ms
         ii = i-is+1
         do  10 j=js,ns
            jj = j-js+1
            ete = ete+zm(ii)*zn(jj)
            ute = ute+usol(i,j)*zm(ii)*zn(jj)
   10    continue
   20 continue
c
c     set perturbation parameter
c
      pertrb = ute/ete
c
c     subtract off constant pertrb
c
      do  40 i=istr,ifnl
         do  30 j=jstr,jfnl
            usol(i,j) = usol(i,j)-pertrb
   30    continue
   40 continue
      return
      end
