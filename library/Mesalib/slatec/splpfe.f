*deck splpfe
      subroutine splpfe (mrelas, nvars, lmx, lbm, ienter, ibasis, imat,
     +   ibrc, ipr, iwr, ind, ibb, erdnrm, eps, gg, dulnrm, dirnrm,
     +   amat, basmat, csc, wr, ww, bl, bu, rz, rg, colnrm, duals,
     +   found)
c***begin prologue  splpfe
c***subsidiary
c***purpose  subsidiary to splp
c***library   slatec
c***type      single precision (splpfe-s, dplpfe-d)
c***author  (unknown)
c***description
c
c     the editing required to convert this subroutine from single to
c     double precision involves the following character string changes.
c
c     use an editing command (change) /string-1/(to)string-2/.
c     /real (12 blanks)/double precision/,/sasum/dasum/,
c     /scopy/dcopy/.
c
c     this subprogram is part of the splp( ) package.
c     it implements the procedure (find variable to enter basis
c     and get search direction).
c     revised 811130-1100
c     revised yymmdd-hhmm
c
c***see also  splp
c***routines called  iploc, la05bs, prwpge, sasum, scopy
c***revision history  (yymmdd)
c   811215  date written
c   890531  changed all specific intrinsics to generic.  (wrb)
c   890605  removed unreferenced labels.  (wrb)
c   891214  prologue converted to version 4.0 format.  (bab)
c   900328  added type section.  (wrb)
c***end prologue  splpfe
      integer ibasis(*),imat(*),ibrc(lbm,2),ipr(*),iwr(*),ind(*),ibb(*)
      real             amat(*),basmat(*),csc(*),wr(*),ww(*),bl(*),bu(*),
     * rz(*),rg(*),colnrm(*),duals(*),cnorm,dirnrm,dulnrm,eps,erdnrm,gg,
     * one,ratio,rcost,rmax,zero
      logical found,trans
c***first executable statement  splpfe
      lpg=lmx-(nvars+4)
      zero=0.e0
      one=1.e0
      rmax=zero
      found=.false.
      i=mrelas+1
      n20002=mrelas+nvars
      go to 20003
20002 i=i+1
20003 if ((n20002-i).lt.0) go to 20004
      j=ibasis(i)
c
c     if j=ibasis(i) .lt. 0 then the variable left at a zero level
c     and is not considered as a candidate to enter.
      if (.not.(j.gt.0)) go to 20006
c
c     do not consider variables corresponding to unbounded step lengths.
      if (.not.(ibb(j).eq.0)) go to 20009
      go to 20002
20009 continue
c
c     if a variable corresponds to an equation(ind=3 and bl=bu),
c     then do not consider it as a candidate to enter.
      if (.not.(ind(j).eq.3)) go to 20012
      if (.not.((bu(j)-bl(j)).le.eps*(abs(bl(j))+abs(bu(j))))) go to 200
     *15
      go to 20002
20015 continue
20012 continue
      rcost=rz(j)
c
c     if variable is at upper bound it can only decrease.  this
c     accounts for the possible change of sign.
      if(mod(ibb(j),2).eq.0) rcost=-rcost
c
c     if the variable is free, use the negative magnitude of the
c     reduced cost for that variable.
      if(ind(j).eq.4) rcost=-abs(rcost)
      cnorm=one
      if(j.le.nvars)cnorm=colnrm(j)
c
c     test for negativity of reduced costs.
      if (.not.(rcost+erdnrm*dulnrm*cnorm.lt.zero)) go to 20018
      found=.true.
      ratio=rcost**2/rg(j)
      if (.not.(ratio.gt.rmax)) go to 20021
      rmax=ratio
      ienter=i
20021 continue
20018 continue
20006 go to 20002
c
c     use col. chosen to compute search direction.
20004 if (.not.(found)) go to 20024
      j=ibasis(ienter)
      ww(1)=zero
      call scopy(mrelas,ww,0,ww,1)
      if (.not.(j.le.nvars)) go to 20027
      if (.not.(j.eq.1)) go to 20030
      ilow=nvars+5
      go to 20031
20030 ilow=imat(j+3)+1
20031 continue
      il1=iploc(ilow,amat,imat)
      if (.not.(il1.ge.lmx-1)) go to 20033
      ilow=ilow+2
      il1=iploc(ilow,amat,imat)
20033 continue
      ipage=abs(imat(lmx-1))
      ihi=imat(j+4)-(ilow-il1)
20036 iu1=min(lmx-2,ihi)
      if (.not.(il1.gt.iu1)) go to 20038
      go to 20037
20038 continue
      do 30 i=il1,iu1
      ww(imat(i))=amat(i)*csc(j)
30    continue
      if (.not.(ihi.le.lmx-2)) go to 20041
      go to 20037
20041 continue
      ipage=ipage+1
      key=1
      call prwpge(key,ipage,lpg,amat,imat)
      il1=nvars+5
      ihi=ihi-lpg
      go to 20036
20037 go to 20028
20027 if (.not.(ind(j).eq.2)) go to 20044
      ww(j-nvars)=one
      go to 20045
20044 ww(j-nvars)=-one
20045 continue
      continue
c
c     compute search direction.
20028 trans=.false.
      call la05bs(basmat,ibrc,lbm,mrelas,ipr,iwr,wr,gg,ww,trans)
c
c     the search direction requires the following sign change if either
c     variable entering is at its upper bound or is free and has
c     positive reduced cost.
      if (.not.(mod(ibb(j),2).eq.0.or.(ind(j).eq.4 .and. rz(j).gt.zero))
     *) go to 20047
      i=1
      n20050=mrelas
      go to 20051
20050 i=i+1
20051 if ((n20050-i).lt.0) go to 20052
      ww(i)=-ww(i)
      go to 20050
20052 continue
20047 dirnrm=sasum(mrelas,ww,1)
c
c     copy contents of wr(*) to duals(*) for use in
c     add-drop (exchange) step, la05cs( ).
      call scopy(mrelas,wr,1,duals,1)
20024 return
      end
