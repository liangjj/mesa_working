*deck spincw
      subroutine spincw (mrelas, nvars, lmx, lbm, npp, jstrt, ibasis,
     +   imat, ibrc, ipr, iwr, ind, ibb, costsc, gg, erdnrm, dulnrm,
     +   amat, basmat, csc, wr, ww, rz, rg, costs, colnrm, duals,
     +   stpedg)
c***begin prologue  spincw
c***subsidiary
c***purpose  subsidiary to splp
c***library   slatec
c***type      single precision (spincw-s, dpincw-d)
c***author  (unknown)
c***description
c
c     the editing required to convert this subroutine from single to
c     double precision involves the following character string changes.
c
c     use an editing command (change) /string-1/(to)string-2/,
c     real (12 blanks)/double precision/,/scopy/dcopy/,/sdot/ddot/.
c
c     this subprogram is part of the splp( ) package.
c     it implements the procedure (initialize reduced costs and
c     steepest edge weights).
c
c***see also  splp
c***routines called  iploc, la05bs, prwpge, scopy, sdot
c***revision history  (yymmdd)
c   811215  date written
c   890531  changed all specific intrinsics to generic.  (wrb)
c   890605  removed unreferenced labels.  (wrb)
c   891214  prologue converted to version 4.0 format.  (bab)
c   900328  added type section.  (wrb)
c***end prologue  spincw
      integer ibasis(*),imat(*),ibrc(lbm,2),ipr(*),iwr(*),ind(*),ibb(*)
      real             amat(*),basmat(*),csc(*),wr(*),ww(*),rz(*),rg(*),
     * costs(*),colnrm(*),duals(*),costsc,erdnrm,dulnrm,gg,one,rzj,
     * scalr,zero,rcost
      logical stpedg,pagepl,trans
c***first executable statement  spincw
      lpg=lmx-(nvars+4)
      zero=0.
      one=1.
c
c     form reduced costs, rz(*), and steepest edge weights, rg(*).
      pagepl=.true.
      rz(1)=zero
      call scopy(nvars+mrelas,rz,0,rz,1)
      rg(1)=one
      call scopy(nvars+mrelas,rg,0,rg,1)
      nnegrc=0
      j=jstrt
20002 if (.not.(ibb(j).le.0)) go to 20004
      pagepl=.true.
      go to 20005
c
c     these are nonbasic independent variables. the cols. are in sparse
c     matrix format.
20004 if (.not.(j.le.nvars)) go to 20007
      rzj=costsc*costs(j)
      ww(1)=zero
      call scopy(mrelas,ww,0,ww,1)
      if (.not.(j.eq.1)) go to 20010
      ilow=nvars+5
      go to 20011
20010 ilow=imat(j+3)+1
20011 continue
      if (.not.(pagepl)) go to 20013
      il1=iploc(ilow,amat,imat)
      if (.not.(il1.ge.lmx-1)) go to 20016
      ilow=ilow+2
      il1=iploc(ilow,amat,imat)
20016 continue
      ipage=abs(imat(lmx-1))
      go to 20014
20013 il1=ihi+1
20014 continue
      ihi=imat(j+4)-(ilow-il1)
20019 iu1=min(lmx-2,ihi)
      if (.not.(il1.gt.iu1)) go to 20021
      go to 20020
20021 continue
      do 60 i=il1,iu1
      rzj=rzj-amat(i)*duals(imat(i))
      ww(imat(i))=amat(i)*csc(j)
60    continue
      if (.not.(ihi.le.lmx-2)) go to 20024
      go to 20020
20024 continue
      ipage=ipage+1
      key=1
      call prwpge(key,ipage,lpg,amat,imat)
      il1=nvars+5
      ihi=ihi-lpg
      go to 20019
20020 pagepl=ihi.eq.(lmx-2)
      rz(j)=rzj*csc(j)
      if (.not.(stpedg)) go to 20027
      trans=.false.
      call la05bs(basmat,ibrc,lbm,mrelas,ipr,iwr,wr,gg,ww,trans)
      rg(j)=sdot(mrelas,ww,1,ww,1)+one
20027 continue
c
c     these are nonbasic dependent variables. the cols. are implicitly
c     defined.
      go to 20008
20007 pagepl=.true.
      ww(1)=zero
      call scopy(mrelas,ww,0,ww,1)
      scalr=-one
      if (ind(j).eq.2) scalr=one
      i=j-nvars
      rz(j)=-scalr*duals(i)
      ww(i)=scalr
      if (.not.(stpedg)) go to 20030
      trans=.false.
      call la05bs(basmat,ibrc,lbm,mrelas,ipr,iwr,wr,gg,ww,trans)
      rg(j)=sdot(mrelas,ww,1,ww,1)+one
20030 continue
      continue
20008 continue
c
20005 rcost=rz(j)
      if (mod(ibb(j),2).eq.0) rcost=-rcost
      if (ind(j).eq.4) rcost=-abs(rcost)
      cnorm=one
      if (j.le.nvars) cnorm=colnrm(j)
      if (rcost+erdnrm*dulnrm*cnorm.lt.zero) nnegrc=nnegrc+1
      j=mod(j,mrelas+nvars)+1
      if (.not.(nnegrc.ge.npp .or. j.eq.jstrt)) go to 20033
      go to 20003
20033 go to 20002
20003 jstrt=j
      return
      end
