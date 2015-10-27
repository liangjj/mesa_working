*deck dplpce
      subroutine dplpce (mrelas, nvars, lmx, lbm, itlp, itbrc, ibasis,
     +   imat, ibrc, ipr, iwr, ind, ibb, erdnrm, eps, tune, gg, amat,
     +   basmat, csc, wr, ww, primal, erd, erp, singlr, redbas)
c***begin prologue  dplpce
c***subsidiary
c***purpose  subsidiary to dsplp
c***library   slatec
c***type      double precision (splpce-s, dplpce-d)
c***author  (unknown)
c***description
c
c     the editing required to convert this subroutine from single to
c     double precision involves the following character string changes.
c
c     use an editing command (change) /string-1/(to)string-2/.
c     /real (12 blanks)/double precision/,
c     /sasum/dasum/,/dcopy/,dcopy/.
c
c     revised 811219-1630
c     revised yymmdd-hhmm
c
c     this subprogram is from the dsplp( ) package.  it calculates
c     the approximate error in the primal and dual systems.  it is
c     the main part of the procedure (compute error in dual and primal
c     systems).
c
c***see also  dsplp
c***routines called  dasum, dcopy, dprwpg, idloc, la05bd
c***revision history  (yymmdd)
c   811215  date written
c   890531  changed all specific intrinsics to generic.  (wrb)
c   890605  removed unreferenced labels.  (wrb)
c   890606  changed references from iploc to idloc.  (wrb)
c   891214  prologue converted to version 4.0 format.  (bab)
c   900328  added type section.  (wrb)
c***end prologue  dplpce
      integer ibasis(*),imat(*),ibrc(lbm,2),ipr(*),iwr(*),ind(*),ibb(*)
      double precision amat(*),basmat(*),csc(*),wr(*),ww(*),primal(*),
     * erd(*),erp(*),eps,erdnrm,factor,gg,one,zero,ten,tune
      double precision dasum
      logical singlr,redbas,trans,pagepl
c***first executable statement  dplpce
      zero=0.d0
      one=1.d0
      ten=10.d0
      lpg=lmx-(nvars+4)
      singlr=.false.
      factor=0.01
c
c     copy colsums in ww(*), and solve transposed system.
      i=1
      n20002=mrelas
      go to 20003
20002 i=i+1
20003 if ((n20002-i).lt.0) go to 20004
      j=ibasis(i)
      if (.not.(j.le.nvars)) go to 20006
      ww(i) = primal(j)
      go to 20007
20006 if (.not.(ind(j).eq.2)) go to 20009
      ww(i)=one
      go to 20010
20009 ww(i)=-one
20010 continue
20007 continue
      go to 20002
c
c     perturb right-side in units of last bits to better reflect
c     errors in the check sum solns.
20004 i=1
      n20012=mrelas
      go to 20013
20012 i=i+1
20013 if ((n20012-i).lt.0) go to 20014
      ww(i)=ww(i)+ten*eps*ww(i)
      go to 20012
20014 trans = .true.
      call la05bd(basmat,ibrc,lbm,mrelas,ipr,iwr,wr,gg,ww,trans)
      i=1
      n20016=mrelas
      go to 20017
20016 i=i+1
20017 if ((n20016-i).lt.0) go to 20018
      erd(i)=max(abs(ww(i)-one),eps)*tune
c
c     system becomes singular when accuracy of solution is .gt. factor.
c     this value (factor) might need to be changed.
      singlr=singlr.or.(erd(i).ge.factor)
      go to 20016
20018 erdnrm=dasum(mrelas,erd,1)
c
c     recalculate row check sums every itbrc iterations or when
c     a redecomposition has occurred.
      if (.not.(mod(itlp,itbrc).eq.0 .or. redbas)) go to 20020
c
c     compute row sums, store in ww(*), solve primal system.
      ww(1)=zero
      call dcopy(mrelas,ww,0,ww,1)
      pagepl=.true.
      j=1
      n20023=nvars
      go to 20024
20023 j=j+1
20024 if ((n20023-j).lt.0) go to 20025
      if (.not.(ibb(j).ge.zero)) go to 20027
c
c     the variable is non-basic.
      pagepl=.true.
      go to 20023
20027 if (.not.(j.eq.1)) go to 20030
      ilow=nvars+5
      go to 20031
20030 ilow=imat(j+3)+1
20031 if (.not.(pagepl)) go to 20033
      il1=idloc(ilow,amat,imat)
      if (.not.(il1.ge.lmx-1)) go to 20036
      ilow=ilow+2
      il1=idloc(ilow,amat,imat)
20036 continue
      ipage=abs(imat(lmx-1))
      go to 20034
20033 il1=ihi+1
20034 ihi=imat(j+4)-(ilow-il1)
20039 iu1=min(lmx-2,ihi)
      if (.not.(il1.gt.iu1)) go to 20041
      go to 20040
20041 continue
      do 20 i=il1,iu1
      ww(imat(i))=ww(imat(i))+amat(i)*csc(j)
20    continue
      if (.not.(ihi.le.lmx-2)) go to 20044
      go to 20040
20044 continue
      ipage=ipage+1
      key=1
      call dprwpg(key,ipage,lpg,amat,imat)
      il1=nvars+5
      ihi=ihi-lpg
      go to 20039
20040 pagepl=ihi.eq.(lmx-2)
      go to 20023
20025 l=1
      n20047=mrelas
      go to 20048
20047 l=l+1
20048 if ((n20047-l).lt.0) go to 20049
      j=ibasis(l)
      if (.not.(j.gt.nvars)) go to 20051
      i=j-nvars
      if (.not.(ind(j).eq.2)) go to 20054
      ww(i)=ww(i)+one
      go to 20055
20054 ww(i)=ww(i)-one
20055 continue
      continue
20051 continue
      go to 20047
c
c     perturb right-side in units of last bit positions.
20049 i=1
      n20057=mrelas
      go to 20058
20057 i=i+1
20058 if ((n20057-i).lt.0) go to 20059
      ww(i)=ww(i)+ten*eps*ww(i)
      go to 20057
20059 trans = .false.
      call la05bd(basmat,ibrc,lbm,mrelas,ipr,iwr,wr,gg,ww,trans)
      i=1
      n20061=mrelas
      go to 20062
20061 i=i+1
20062 if ((n20061-i).lt.0) go to 20063
      erp(i)=max(abs(ww(i)-one),eps)*tune
c
c     system becomes singular when accuracy of solution is .gt. factor.
c     this value (factor) might need to be changed.
      singlr=singlr.or.(erp(i).ge.factor)
      go to 20061
20063 continue
c
20020 return
      end
