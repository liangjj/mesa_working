*deck dplpmu
      subroutine dplpmu (mrelas, nvars, lmx, lbm, nredc, info, ienter,
     +   ileave, iopt, npp, jstrt, ibasis, imat, ibrc, ipr, iwr, ind,
     +   ibb, anorm, eps, uu, gg, rprnrm, erdnrm, dulnrm, theta, costsc,
     +   xlamda, rhsnrm, amat, basmat, csc, wr, rprim, ww, bu, bl, rhs,
     +   erd, erp, rz, rg, colnrm, costs, primal, duals, singlr, redbas,
     +   zerolv, stpedg)
c***begin prologue  dplpmu
c***subsidiary
c***purpose  subsidiary to dsplp
c***library   slatec
c***type      double precision (splpmu-s, dplpmu-d)
c***author  (unknown)
c***description
c
c     the editing required to convert this subroutine from single to
c     double precision involves the following character string changes.
c
c     use an editing command (change) /string-1/(to)string-2/.
c     /real (12 blanks)/double precision/,
c     /sasum/dasum/,/scopy/dcopy/,/sdot/ddot/,
c     /.e0/.d0/
c
c     this subprogram is from the dsplp( ) package.  it performs the
c     tasks of updating the primal solution, edge weights, reduced
c     costs, and matrix decomposition.
c     it is the main part of the procedure (make move and update).
c
c     revised 821122-1100
c     revised yymmdd
c
c***see also  dsplp
c***routines called  dasum, dcopy, ddot, dplpdm, dpnnzr, dprwpg, idloc,
c                    la05bd, la05cd, xermsg
c***revision history  (yymmdd)
c   811215  date written
c   890531  changed all specific intrinsics to generic.  (wrb)
c   890605  removed unreferenced labels.  (wrb)
c   890606  changed references from iploc to idloc.  (wrb)
c   890606  removed unused common block la05dd.  (wrb)
c   891214  prologue converted to version 4.0 format.  (bab)
c   900315  calls to xerror changed to calls to xermsg.  (thj)
c   900328  added type section.  (wrb)
c***end prologue  dplpmu
      integer ibasis(*),imat(*),ibrc(lbm,2),ipr(*),iwr(*),ind(*),ibb(*)
      double precision aij,alpha,anorm,costsc,erdnrm,dulnrm,eps,gamma,
     * gg,gq,one,rprnrm,rzj,scalr,theta,two,uu,wp,xlamda,rhsnrm,
     * zero,amat(*),basmat(*),csc(*),wr(*),rprim(*),ww(*),bu(*),bl(*),
     * rhs(*),erd(*),erp(*),rz(*),rg(*),costs(*),primal(*),duals(*),
     * colnrm(*),rcost,dasum,ddot,cnorm
      logical singlr,redbas,pagepl,trans,zerolv,stpedg
c
c***first executable statement  dplpmu
      zero=0.d0
      one=1.d0
      two=2.d0
      lpg=lmx-(nvars+4)
c
c     update the primal solution with a multiple of the search
c     direction.
      i=1
      n20002=mrelas
      go to 20003
20002 i=i+1
20003 if ((n20002-i).lt.0) go to 20004
      rprim(i)=rprim(i)-theta*ww(i)
      go to 20002
c
c     if ejected variable is leaving at an upper bound,  then
c     translate right hand side.
20004 if (.not.(ileave.lt.0)) go to 20006
      ibas=ibasis(abs(ileave))
      scalr=rprim(abs(ileave))
      assign 20009 to npr001
      go to 30001
20009 ibb(ibas)=abs(ibb(ibas))+1
c
c     if entering variable is restricted to its upper bound, translate
c     right hand side.  if the variable decreased from its upper
c     bound, a sign change is required in the translation.
20006 if (.not.(ienter.eq.ileave)) go to 20010
      ibas=ibasis(ienter)
      scalr=theta
      if (mod(ibb(ibas),2).eq.0) scalr=-scalr
      assign 20013 to npr001
      go to 30001
20013 ibb(ibas)=ibb(ibas)+1
      go to 20011
20010 ibas=ibasis(ienter)
c
c     if entering variable is decreasing from its upper bound,
c     complement its primal value.
      if (.not.(ind(ibas).eq.3.and.mod(ibb(ibas),2).eq.0)) go to 20014
      scalr=-(bu(ibas)-bl(ibas))
      if (ibas.le.nvars) scalr=scalr/csc(ibas)
      assign 20017 to npr001
      go to 30001
20017 theta=-scalr-theta
      ibb(ibas)=ibb(ibas)+1
20014 continue
      rprim(abs(ileave))=theta
      ibb(ibas)=-abs(ibb(ibas))
      i=ibasis(abs(ileave))
      ibb(i)=abs(ibb(i))
      if(primal(abs(ileave)+nvars).gt.zero) ibb(i)=ibb(i)+1
c
c     interchange column pointers to note exchange of columns.
20011 ibas=ibasis(ienter)
      ibasis(ienter)=ibasis(abs(ileave))
      ibasis(abs(ileave))=ibas
c
c     if variable was exchanged at a zero level, mark it so that
c     it can't be brought back in.  this is to help prevent cycling.
      if(zerolv) ibasis(ienter)=-abs(ibasis(ienter))
      rprnrm=max(rprnrm,dasum(mrelas,rprim,1))
      k=1
      n20018=mrelas
      go to 20019
20018 k=k+1
20019 if ((n20018-k).lt.0) go to 20020
c
c     see if variables that were classified as infeasible have now
c     become feasible.  this may required translating upper bounded
c     variables.
      if (.not.(primal(k+nvars).ne.zero .and.
     *          abs(rprim(k)).le.rprnrm*erp(k))) go to 20022
      if (.not.(primal(k+nvars).gt.zero)) go to 20025
      ibas=ibasis(k)
      scalr=-(bu(ibas)-bl(ibas))
      if(ibas.le.nvars)scalr=scalr/csc(ibas)
      assign 20028 to npr001
      go to 30001
20028 rprim(k)=-scalr
      rprnrm=rprnrm-scalr
20025 primal(k+nvars)=zero
20022 continue
      go to 20018
c
c     update reduced costs, edge weights, and matrix decomposition.
20020 if (.not.(ienter.ne.ileave)) go to 20029
c
c     the incoming variable is always classified as feasible.
      primal(abs(ileave)+nvars)=zero
c
      wp=ww(abs(ileave))
      gq=ddot(mrelas,ww,1,ww,1)+one
c
c     compute inverse (transpose) times search direction.
      trans=.true.
      call la05bd(basmat,ibrc,lbm,mrelas,ipr,iwr,wr,gg,ww,trans)
c
c     update the matrix decomposition.  col. abs(ileave) is leaving.
c     the array duals(*) contains intermediate results for the
c     incoming column.
      call la05cd(basmat,ibrc,lbm,mrelas,ipr,iwr,duals,gg,uu,
     * abs(ileave))
      redbas=.false.
      if (.not.(gg.lt.zero)) go to 20032
c
c     redecompose basis matrix when an error return from
c     la05cd( ) is noted.  this will probably be due to
c     space being exhausted, gg=-7.
      call dplpdm(
     *mrelas,nvars,lmx,lbm,nredc,info,iopt,
     *ibasis,imat,ibrc,ipr,iwr,ind,ibb,
     *anorm,eps,uu,gg,
     *amat,basmat,csc,wr,
     *singlr,redbas)
      if (.not.(singlr)) go to 20035
      nerr=26
      call xermsg ('slatec', 'dplpmu',
     +   'in dsplp, moved to a singular point. this should not happen.',
     +   nerr, iopt)
      info=-nerr
      return
20035 continue
      go to 30002
20038 continue
20032 continue
c
c     if steepest edge pricing is used, update reduced costs
c     and edge weights.
      if (.not.(stpedg)) go to 20039
c
c     compute col. abs(ileave) of the new inverse (transpose) matrix
c     here abs(ileave) points to the ejected column.
c     use erd(*) for temp. storage.
      call dcopy(mrelas,zero,0,erd,1)
      erd(abs(ileave))=one
      trans=.true.
      call la05bd(basmat,ibrc,lbm,mrelas,ipr,iwr,wr,gg,erd,trans)
c
c     compute updated dual variables in duals(*).
      assign 20042 to npr003
      go to 30003
c
c     compute the dot product of col. j of the new inverse (transpose)
c     with each non-basic column.  also compute the dot product of the
c     inverse (transpose) of non-updated matrix (times) the
c     search direction with each non-basic column.
c     recompute reduced costs.
20042 pagepl=.true.
      call dcopy(nvars+mrelas,zero,0,rz,1)
      nnegrc=0
      j=jstrt
20043 if (.not.(ibb(j).le.0)) go to 20045
      pagepl=.true.
      rg(j)=one
      go to 20046
c
c     nonbasic independent variables (column in sparse matrix storage)
20045 if (.not.(j.le.nvars)) go to 20048
      rzj=costs(j)*costsc
      alpha=zero
      gamma=zero
c
c     compute the dot product of the sparse matrix nonbasic columns
c     with three vectors involved in the updating step.
      if (.not.(j.eq.1)) go to 20051
      ilow=nvars+5
      go to 20052
20051 ilow=imat(j+3)+1
20052 if (.not.(pagepl)) go to 20054
      il1=idloc(ilow,amat,imat)
      if (.not.(il1.ge.lmx-1)) go to 20057
      ilow=ilow+2
      il1=idloc(ilow,amat,imat)
20057 continue
      ipage=abs(imat(lmx-1))
      go to 20055
20054 il1=ihi+1
20055 ihi=imat(j+4)-(ilow-il1)
20060 iu1=min(lmx-2,ihi)
      if (.not.(il1.gt.iu1)) go to 20062
      go to 20061
20062 continue
      do 10 i=il1,iu1
      rzj=rzj-amat(i)*duals(imat(i))
      alpha=alpha+amat(i)*erd(imat(i))
      gamma=gamma+amat(i)*ww(imat(i))
10    continue
      if (.not.(ihi.le.lmx-2)) go to 20065
      go to 20061
20065 continue
      ipage=ipage+1
      key=1
      call dprwpg(key,ipage,lpg,amat,imat)
      il1=nvars+5
      ihi=ihi-lpg
      go to 20060
20061 pagepl=ihi.eq.(lmx-2)
      rz(j)=rzj*csc(j)
      alpha=alpha*csc(j)
      gamma=gamma*csc(j)
      rg(j)=max(rg(j)-two*alpha*gamma+alpha**2*gq,one+alpha**2)
c
c     nonbasic dependent variables (columns defined implicitly)
      go to 20049
20048 pagepl=.true.
      scalr=-one
      if(ind(j).eq.2) scalr=one
      i=j-nvars
      alpha=scalr*erd(i)
      rz(j)=-scalr*duals(i)
      gamma=scalr*ww(i)
      rg(j)=max(rg(j)-two*alpha*gamma+alpha**2*gq,one+alpha**2)
20049 continue
20046 continue
c
      rcost=rz(j)
      if (mod(ibb(j),2).eq.0) rcost=-rcost
      if (.not.(ind(j).eq.3)) go to 20068
      if(bu(j).eq.bl(j)) rcost=zero
20068 continue
      if (ind(j).eq.4) rcost=-abs(rcost)
      cnorm=one
      if (j.le.nvars) cnorm=colnrm(j)
      if (rcost+erdnrm*dulnrm*cnorm.lt.zero) nnegrc=nnegrc+1
      j=mod(j,mrelas+nvars)+1
      if (.not.(nnegrc.ge.npp .or. j.eq.jstrt)) go to 20071
      go to 20044
20071 continue
      go to 20043
20044 jstrt=j
c
c     update the edge weight for the ejected variable.
      rg(abs(ibasis(ienter)))= gq/wp**2
c
c     if minimum reduced cost (dantzig) pricing is used,
c     calculate the new reduced costs.
      go to 20040
c
c     compute the updated duals in duals(*).
20039 assign 20074 to npr003
      go to 30003
20074 call dcopy(nvars+mrelas,zero,0,rz,1)
      nnegrc=0
      j=jstrt
      pagepl=.true.
c
20075 if (.not.(ibb(j).le.0)) go to 20077
      pagepl=.true.
      go to 20078
c
c     nonbasic independent variable (column in sparse matrix storage)
20077 if (.not.(j.le.nvars)) go to 20080
      rz(j)=costs(j)*costsc
      if (.not.(j.eq.1)) go to 20083
      ilow=nvars+5
      go to 20084
20083 ilow=imat(j+3)+1
20084 continue
      if (.not.(pagepl)) go to 20086
      il1=idloc(ilow,amat,imat)
      if (.not.(il1.ge.lmx-1)) go to 20089
      ilow=ilow+2
      il1=idloc(ilow,amat,imat)
20089 continue
      ipage=abs(imat(lmx-1))
      go to 20087
20086 il1=ihi+1
20087 continue
      ihi=imat(j+4)-(ilow-il1)
20092 iu1=min(lmx-2,ihi)
      if (.not.(iu1.ge.il1 .and.mod(iu1-il1,2).eq.0)) go to 20094
      rz(j)=rz(j)-amat(il1)*duals(imat(il1))
      il1=il1+1
20094 continue
      if (.not.(il1.gt.iu1)) go to 20097
      go to 20093
20097 continue
c
c     unroll the dot product loop to a depth of two.  (this is done
c     for increased efficiency).
      do 40 i=il1,iu1,2
      rz(j)=rz(j)-amat(i)*duals(imat(i))-amat(i+1)*duals(imat(i+1))
40    continue
      if (.not.(ihi.le.lmx-2)) go to 20100
      go to 20093
20100 continue
      ipage=ipage+1
      key=1
      call dprwpg(key,ipage,lpg,amat,imat)
      il1=nvars+5
      ihi=ihi-lpg
      go to 20092
20093 pagepl=ihi.eq.(lmx-2)
      rz(j)=rz(j)*csc(j)
c
c     nonbasic dependent variables (columns defined implicitly)
      go to 20081
20080 pagepl=.true.
      scalr=-one
      if(ind(j).eq.2) scalr=one
      i=j-nvars
      rz(j)=-scalr*duals(i)
20081 continue
20078 continue
c
      rcost=rz(j)
      if (mod(ibb(j),2).eq.0) rcost=-rcost
      if (.not.(ind(j).eq.3)) go to 20103
      if(bu(j).eq.bl(j)) rcost=zero
20103 continue
      if (ind(j).eq.4) rcost=-abs(rcost)
      cnorm=one
      if (j.le.nvars) cnorm=colnrm(j)
      if (rcost+erdnrm*dulnrm*cnorm.lt.zero) nnegrc=nnegrc+1
      j=mod(j,mrelas+nvars)+1
      if (.not.(nnegrc.ge.npp .or. j.eq.jstrt)) go to 20106
      go to 20076
20106 continue
      go to 20075
20076 jstrt=j
20040 continue
      go to 20030
c
c     this is necessary only for printing of intermediate results.
20029 assign 20109 to npr003
      go to 30003
20109 continue
20030 return
c     procedure (translate right hand side)
c
c     perform the translation on the right-hand side.
30001 if (.not.(ibas.le.nvars)) go to 20110
      i=0
20113 call dpnnzr(i,aij,iplace,amat,imat,ibas)
      if (.not.(i.le.0)) go to 20115
      go to 20114
20115 continue
      rhs(i)=rhs(i)-scalr*aij*csc(ibas)
      go to 20113
20114 go to 20111
20110 i=ibas-nvars
      if (.not.(ind(ibas).eq.2)) go to 20118
      rhs(i)=rhs(i)-scalr
      go to 20119
20118 rhs(i)=rhs(i)+scalr
20119 continue
20111 continue
      rhsnrm=max(rhsnrm,dasum(mrelas,rhs,1))
      go to npr001, (20009,20013,20017,20028)
c     procedure (compute new primal)
c
c     copy rhs into ww(*), solve system.
30002 call dcopy(mrelas,rhs,1,ww,1)
      trans = .false.
      call la05bd(basmat,ibrc,lbm,mrelas,ipr,iwr,wr,gg,ww,trans)
      call dcopy(mrelas,ww,1,rprim,1)
      rprnrm=dasum(mrelas,rprim,1)
      go to 20038
c     procedure (compute new duals)
c
c     solve for dual variables. first copy costs into duals(*).
30003 i=1
      n20121=mrelas
      go to 20122
20121 i=i+1
20122 if ((n20121-i).lt.0) go to 20123
      j=ibasis(i)
      if (.not.(j.le.nvars)) go to 20125
      duals(i)=costsc*costs(j)*csc(j) + xlamda*primal(i+nvars)
      go to 20126
20125 duals(i)=xlamda*primal(i+nvars)
20126 continue
      go to 20121
c
20123 trans=.true.
      call la05bd(basmat,ibrc,lbm,mrelas,ipr,iwr,wr,gg,duals,trans)
      dulnrm=dasum(mrelas,duals,1)
      go to npr003, (20042,20074,20109)
      end
