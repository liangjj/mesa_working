*deck splpmn
      subroutine splpmn (usrmat, mrelas, nvars, costs, prgopt, dattrv,
     +   bl, bu, ind, info, primal, duals, amat, csc, colnrm, erd, erp,
     +   basmat, wr, rz, rg, rprim, rhs, ww, lmx, lbm, ibasis, ibb,
     +   imat, ibrc, ipr, iwr)
c***begin prologue  splpmn
c***subsidiary
c***purpose  subsidiary to splp
c***library   slatec
c***type      single precision (splpmn-s, dplpmn-d)
c***author  (unknown)
c***description
c
c     marvel option(s).. output=yes/no to eliminate printed output.
c     this does not apply to the calls to the error processor.
c
c     main subroutine for splp package.
c
c***see also  splp
c***routines called  ivout, la05bs, pinitm, pnnzrs, prwpge, sasum,
c                    sclosm, scopy, sdot, spincw, spinit, splpce,
c                    splpdm, splpfe, splpfl, splpmu, splpup, spopt,
c                    svout, xermsg
c***common blocks    la05ds
c***revision history  (yymmdd)
c   811215  date written
c   890531  changed all specific intrinsics to generic.  (wrb)
c   890605  corrected references to xerrwv.  (wrb)
c   890605  removed unreferenced labels.  (wrb)
c   891009  removed unreferenced variable.  (wrb)
c   891214  prologue converted to version 4.0 format.  (bab)
c   900315  calls to xerror changed to calls to xermsg.  (thj)
c   900328  added type section.  (wrb)
c   900510  convert xerrwv calls to xermsg calls.  (rwc)
c***end prologue  splpmn
      real             abig,aij,amat(*),anorm,asmall,basmat(*),
     * bl(*),bu(*),colnrm(*),costs(*),costsc,csc(*),dattrv(*),
     * dirnrm,duals(*),dulnrm,eps,tune,erd(*),erdnrm,erp(*),factor,gg,
     * one,prgopt(*),primal(*),resnrm,rg(*),rhs(*),rhsnrm,ropt(07),
     * rprim(*),rprnrm,rz(*),rzj,scalr,scosts,size,small,theta,
     * tolls,upbnd,uu,wr(*),ww(*),xlamda,xval,zero,rdum(01),tolabs
c
      integer ibasis(*),ibb(*),ibrc(lbm,2),imat(*),ind(*),
     * ipr(*),iwr(*),intopt(08),idum(01)
c
c     array local variables
c     name(length)          description
c
c     costs(nvars)          cost coefficients
c     prgopt( )             option vector
c     dattrv( )             data transfer vector
c     primal(nvars+mrelas)  as output it is primal solution of lp.
c                           internally, the first nvars positions hold
c                           the column check sums.  the next mrelas
c                           positions hold the classification for the
c                           basic variables  -1 violates lower
c                           bound, 0 feasible, +1 violates upper bound
c     duals(mrelas+nvars)   dual solution. internally holds r.h. side
c                           as first mrelas entries.
c     amat(lmx)             sparse form of data matrix
c     imat(lmx)             sparse form of data matrix
c     bl(nvars+mrelas)      lower bounds for variables
c     bu(nvars+mrelas)      upper bounds for variables
c     ind(nvars+mrelas)     indicator for variables
c     csc(nvars)            column scaling
c     ibasis(nvars+mrelas)  cols. 1-mrelas are basic, rest are non-basic
c     ibb(nvars+mrelas)     indicator for non-basic vars., polarity of
c                           vars., and potentially infinite vars.
c                           if ibb(j).lt.0, variable j is basic
c                           if ibb(j).gt.0, variable j is non-basic
c                           if ibb(j).eq.0, variable j has to be ignored
c                           because it would cause unbounded soln.
c                           when mod(ibb(j),2).eq.0, variable is at its
c                           upper bound, otherwise it is at its lower
c                           bound
c     colnrm(nvars)         norm of columns
c     erd(mrelas)           errors in dual variables
c     erp(mrelas)           errors in primal variables
c     basmat(lbm)           basis matrix for harwell sparse code
c     ibrc(lbm,2)           row and column pointers for basmat(*)
c     ipr(2*mrelas)         work array for harwell sparse code
c     iwr(8*mrelas)         work array for harwell sparse code
c     wr(mrelas)            work array for harwell sparse code
c     rz(nvars+mrelas)      reduced costs
c     rprim(mrelas)         internal primal solution
c     rg(nvars+mrelas)      column weights
c     ww(mrelas)            work array
c     rhs(mrelas)           holds translated right hand side
c
c     scalar local variables
c     name       type         description
c
c     lmx        integer      length of amat(*)
c     lpg        integer      length of page for amat(*)
c     eps        real         machine precision
c     tune       real         parameter to scale error estimates
c     tolls      real         relative tolerance for small residuals
c     tolabs     real         absolute tolerance for small residuals.
c                             used if relative error test fails.
c                             in constraint equations
c     factor     real         .01--determines if basis is singular
c                             or component is feasible.  may need to
c                             be increased to 1.e0 on short word
c                             length machines.
c     asmall     real         lower bound for non-zero magn. in amat(*)
c     abig       real         upper bound for non-zero magn. in amat(*)
c     mxitlp     integer      maximum number of iterations for lp
c     itlp       integer      iteration counter for total lp iters
c     costsc     real         costs(*) scaling
c     scosts     real         temp loc. for costsc.
c     xlamda     real         weight parameter for pen. method.
c     anorm      real         norm of data matrix amat(*)
c     rprnrm     real         norm of the solution
c     dulnrm     real         norm of the duals
c     erdnrm     real         norm of error in dual variables
c     dirnrm     real         norm of the direction vector
c     rhsnrm     real         norm of translated right hand side vector
c     resnrm     real         norm of residual vector for checking
c                             feasibility
c     nzbm       integer      number of non-zeros in basmat(*)
c     lbm        integer      length of basmat(*)
c     small      real         eps*anorm  used in harwell sparse code
c     lp         integer      used in harwell la05*() pack as output
c                             file number. set=i1mach(4) now.
c     uu         real         0.1--used in harwell sparse code
c                             for relative pivoting tolerance.
c     gg         real         output info flag in harwell sparse code
c     iplace     integer      integer used by sparse matrix codes
c     ienter     integer      next column to enter basis
c     nredc      integer      no. of full redecompositions
c     kprint     integer      level of output, =0-3
c     idg        integer      format and precision of output
c     itbrc      integer      no. of iters. between recalculating
c                             the error in the primal solution.
c     npp        integer      no. of negative reduced costs required
c                             in partial pricing
c     jstrt      integer      starting place for partial pricing.
c
      logical colscp,savedt,contin,cstscp,unbnd,
     *        feas,finite,found,minprb,redbas,
     *        singlr,sizeup,stpedg,trans,usrbas,zerolv,lopt(08)
      character*8 xern1, xern2
      equivalence (contin,lopt(1)),(usrbas,lopt(2)),
     *  (sizeup,lopt(3)),(savedt,lopt(4)),(colscp,lopt(5)),
     *  (cstscp,lopt(6)),(minprb,lopt(7)),(stpedg,lopt(8)),
     *  (idg,intopt(1)),(ipagef,intopt(2)),(isave,intopt(3)),
     *  (mxitlp,intopt(4)),(kprint,intopt(5)),(itbrc,intopt(6)),
     *  (npp,intopt(7)),(lprg,intopt(8)),(eps,ropt(1)),(asmall,ropt(2)),
     *  (abig,ropt(3)),(costsc,ropt(4)),(tolls,ropt(5)),(tune,ropt(6)),
     *   (tolabs,ropt(7))
c
c     common block used by la05 () package..
      common /la05ds/ small,lp,lenl,lenu,ncp,lrow,lcol
      external usrmat
c
c     set lp=0 so no error messages will print within la05 () package.
c***first executable statement  splpmn
      lp=0
c
c     the values zero and one.
      zero=0.e0
      one=1.e0
      factor=0.01e0
      lpg=lmx-(nvars+4)
      iopt=1
      info=0
      unbnd=.false.
      jstrt=1
c
c     process user options in prgopt(*).
c     check that any user-given changes are well-defined.
      call spopt(prgopt,mrelas,nvars,info,csc,ibasis,ropt,intopt,lopt)
      if (.not.(info.lt.0)) go to 20002
      go to 30001
20002 if (.not.(contin)) go to 20003
      go to 30002
20006 go to 20004
c
c     initialize sparse data matrix, amat(*) and imat(*).
20003 call pinitm(mrelas,nvars,amat,imat,lmx,ipagef)
c
c     update matrix data and check bounds for consistency.
20004 call splpup(usrmat,mrelas,nvars,prgopt,dattrv,
     *     bl,bu,ind,info,amat,imat,sizeup,asmall,abig)
      if (.not.(info.lt.0)) go to 20007
      go to 30001
c
c++  code for output=yes is active
20007 if (.not.(kprint.ge.1)) go to 20008
      go to 30003
20011 continue
c++  code for output=no is inactive
c++  end
c
c     initialization. scale data, normalize bounds, form column
c     check sums, and form initial basis matrix.
20008 call spinit(mrelas,nvars,costs,bl,bu,ind,primal,info,
     * amat,csc,costsc,colnrm,xlamda,anorm,rhs,rhsnrm,
     * ibasis,ibb,imat,lopt)
      if (.not.(info.lt.0)) go to 20012
      go to 30001
c
20012 nredc=0
      assign 20013 to npr004
      go to 30004
20013 if (.not.(singlr)) go to 20014
      nerr=23
      call xermsg ('slatec', 'splpmn',
     +   'in splp,  a singular initial basis was encountered.', nerr,
     +   iopt)
      info=-nerr
      go to 30001
20014 assign 20018 to npr005
      go to 30005
20018 assign 20019 to npr006
      go to 30006
20019 assign 20020 to npr007
      go to 30007
20020 if (.not.(usrbas)) go to 20021
      assign 20024 to npr008
      go to 30008
20024 if (.not.(.not.feas)) go to 20025
      nerr=24
      call xermsg ('slatec', 'splpmn',
     +   'in splp, an infeasible initial basis was encountered.', nerr,
     +   iopt)
      info=-nerr
      go to 30001
20025 continue
20021 itlp=0
c
c     lamda has been set to a constant, perform penalty method.
      assign 20029 to npr009
      go to 30009
20029 assign 20030 to npr010
      go to 30010
20030 assign 20031 to npr006
      go to 30006
20031 assign 20032 to npr008
      go to 30008
20032 if (.not.(.not.feas)) go to 20033
c
c     set lamda to infinity by setting costsc to zero (save the value of
c     costsc) and perform standard phase-1.
      if(kprint.ge.2)call ivout(0,idum,'('' enter standard phase-1'')',
     *idg)
      scosts=costsc
      costsc=zero
      assign 20036 to npr007
      go to 30007
20036 assign 20037 to npr009
      go to 30009
20037 assign 20038 to npr010
      go to 30010
20038 assign 20039 to npr006
      go to 30006
20039 assign 20040 to npr008
      go to 30008
20040 if (.not.(feas)) go to 20041
c
c     set lamda to zero, costsc=scosts, perform standard phase-2.
      if(kprint.gt.1)call ivout(0,idum,'('' enter standard phase-2'')',
     *idg)
      xlamda=zero
      costsc=scosts
      assign 20044 to npr009
      go to 30009
20044 continue
20041 go to 20034
c     check if any basic variables are still classified as
c     infeasible.  if any are, then this may not yet be an
c     optimal point.  therefore set lamda to zero and try
c     to perform more simplex steps.
20033 i=1
      n20046=mrelas
      go to 20047
20046 i=i+1
20047 if ((n20046-i).lt.0) go to 20048
      if (primal(i+nvars).ne.zero) go to 20045
      go to 20046
20048 go to 20035
20045 xlamda=zero
      assign 20050 to npr009
      go to 30009
20050 continue
20034 continue
c
20035 assign 20051 to npr011
      go to 30011
20051 if (.not.(feas.and.(.not.unbnd))) go to 20052
      info=1
      go to 20053
20052 if (.not.((.not.feas).and.(.not.unbnd))) go to 10001
      nerr=1
      call xermsg ('slatec', 'splpmn',
     +   'in splp, the problem appears to be infeasible', nerr, iopt)
      info=-nerr
      go to 20053
10001 if (.not.(feas .and. unbnd)) go to 10002
      nerr=2
      call xermsg ('slatec', 'splpmn',
     +   'in splp, the problem appears to have no finite solution.',
     +   nerr, iopt)
      info=-nerr
      go to 20053
10002 if (.not.((.not.feas).and.unbnd)) go to 10003
      nerr=3
      call xermsg ('slatec', 'splpmn',
     +   'in splp, the problem appears to be infeasible and to have ' //
     +   'no finite solution.', nerr, iopt)
      info=-nerr
10003 continue
20053 continue
c
      if (.not.(info.eq.(-1) .or. info.eq.(-3))) go to 20055
      size=sasum(nvars,primal,1)*anorm
      size=size/sasum(nvars,csc,1)
      size=size+sasum(mrelas,primal(nvars+1),1)
      i=1
      n20058=nvars+mrelas
      go to 20059
20058 i=i+1
20059 if ((n20058-i).lt.0) go to 20060
      nx0066=ind(i)
      if (nx0066.lt.1.or.nx0066.gt.4) go to 20066
      go to (20062,20063,20064,20065), nx0066
20062 if (.not.(size+abs(primal(i)-bl(i))*factor.eq.size)) go to 20068
      go to 20058
20068 if (.not.(primal(i).gt.bl(i))) go to 10004
      go to 20058
10004 ind(i)=-4
      go to 20067
20063 if (.not.(size+abs(primal(i)-bu(i))*factor.eq.size)) go to 20071
      go to 20058
20071 if (.not.(primal(i).lt.bu(i))) go to 10005
      go to 20058
10005 ind(i)=-4
      go to 20067
20064 if (.not.(size+abs(primal(i)-bl(i))*factor.eq.size)) go to 20074
      go to 20058
20074 if (.not.(primal(i).lt.bl(i))) go to 10006
      ind(i)=-4
      go to 20075
10006 if (.not.(size+abs(primal(i)-bu(i))*factor.eq.size)) go to 10007
      go to 20058
10007 if (.not.(primal(i).gt.bu(i))) go to 10008
      ind(i)=-4
      go to 20075
10008 go to 20058
20075 go to 20067
20065 go to 20058
20066 continue
20067 go to 20058
20060 continue
20055 continue
c
      if (.not.(info.eq.(-2) .or. info.eq.(-3))) go to 20077
      j=1
      n20080=nvars
      go to 20081
20080 j=j+1
20081 if ((n20080-j).lt.0) go to 20082
      if (.not.(ibb(j).eq.0)) go to 20084
      nx0091=ind(j)
      if (nx0091.lt.1.or.nx0091.gt.4) go to 20091
      go to (20087,20088,20089,20090), nx0091
20087 bu(j)=bl(j)
      ind(j)=-3
      go to 20092
20088 bl(j)=bu(j)
      ind(j)=-3
      go to 20092
20089 go to 20080
20090 bl(j)=zero
      bu(j)=zero
      ind(j)=-3
20091 continue
20092 continue
20084 go to 20080
20082 continue
20077 continue
c++  code for output=yes is active
      if (.not.(kprint.ge.1)) go to 20093
      assign 20096 to npr012
      go to 30012
20096 continue
20093 continue
c++  code for output=no is inactive
c++  end
      go to 30001
c ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     procedure (compute right hand side)
30010 rhs(1)=zero
      call scopy(mrelas,rhs,0,rhs,1)
      j=1
      n20098=nvars+mrelas
      go to 20099
20098 j=j+1
20099 if ((n20098-j).lt.0) go to 20100
      nx0106=ind(j)
      if (nx0106.lt.1.or.nx0106.gt.4) go to 20106
      go to (20102,20103,20104,20105), nx0106
20102 scalr=-bl(j)
      go to 20107
20103 scalr=-bu(j)
      go to 20107
20104 scalr=-bl(j)
      go to 20107
20105 scalr=zero
20106 continue
20107 if (.not.(scalr.ne.zero)) go to 20108
      if (.not.(j.le.nvars)) go to 20111
      i=0
20114 call pnnzrs(i,aij,iplace,amat,imat,j)
      if (.not.(i.le.0)) go to 20116
      go to 20115
20116 rhs(i)=rhs(i)+aij*scalr
      go to 20114
20115 go to 20112
20111 rhs(j-nvars)=rhs(j-nvars)-scalr
20112 continue
20108 go to 20098
20100 j=1
      n20119=nvars+mrelas
      go to 20120
20119 j=j+1
20120 if ((n20119-j).lt.0) go to 20121
      scalr=zero
      if(ind(j).eq.3.and.mod(ibb(j),2).eq.0) scalr=bu(j)-bl(j)
      if (.not.(scalr.ne.zero)) go to 20123
      if (.not.(j.le.nvars)) go to 20126
      i=0
20129 call pnnzrs(i,aij,iplace,amat,imat,j)
      if (.not.(i.le.0)) go to 20131
      go to 20130
20131 rhs(i)=rhs(i)-aij*scalr
      go to 20129
20130 go to 20127
20126 rhs(j-nvars)=rhs(j-nvars)+scalr
20127 continue
20123 go to 20119
20121 continue
      go to npr010, (20030,20038)
c     procedure (perform simplex steps)
30009 assign 20134 to npr013
      go to 30013
20134 assign 20135 to npr014
      go to 30014
20135 if (.not.(kprint.gt.2)) go to 20136
      call svout(mrelas,duals,'('' basic (internal) dual soln.'')',idg)
      call svout(nvars+mrelas,rz,'('' reduced costs'')',idg)
20136 continue
20139 assign 20141 to npr015
      go to 30015
20141 if (.not.(.not. found)) go to 20142
      go to 30016
20145 continue
20142 if (.not.(found)) go to 20146
      if (kprint.ge.3) call svout(mrelas,ww,'('' search direction'')',
     *idg)
      go to 30017
20149 if (.not.(finite)) go to 20150
      go to 30018
20153 assign 20154 to npr005
      go to 30005
20154 go to 20151
20150 unbnd=.true.
      ibb(ibasis(ienter))=0
20151 go to 20147
20146 go to 20140
20147 itlp=itlp+1
      go to 30019
20155 go to 20139
20140 continue
      go to npr009, (20029,20037,20044,20050)
c ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     procedure (retrieve saved data from file isave)
30002 lpr=nvars+4
      rewind isave
      read(isave) (amat(i),i=1,lpr),(imat(i),i=1,lpr)
      key=2
      ipage=1
      go to 20157
20156 if (np.lt.0) go to 20158
20157 lpr1=lpr+1
      read(isave) (amat(i),i=lpr1,lmx),(imat(i),i=lpr1,lmx)
      call prwpge(key,ipage,lpg,amat,imat)
      np=imat(lmx-1)
      ipage=ipage+1
      go to 20156
20158 nparm=nvars+mrelas
      read(isave) (ibasis(i),i=1,nparm)
      rewind isave
      go to 20006
c ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     procedure (save data on file isave)
c
c     some pages may not be written yet.
30020 if (.not.(amat(lmx).eq.one)) go to 20159
      amat(lmx)=zero
      key=2
      ipage=abs(imat(lmx-1))
      call prwpge(key,ipage,lpg,amat,imat)
c
c     force page file to be opened on restarts.
20159 key=amat(4)
      amat(4)=zero
      lpr=nvars+4
      write(isave) (amat(i),i=1,lpr),(imat(i),i=1,lpr)
      amat(4)=key
      ipage=1
      key=1
      go to 20163
20162 if (np.lt.0) go to 20164
20163 call prwpge(key,ipage,lpg,amat,imat)
      lpr1=lpr+1
      write(isave) (amat(i),i=lpr1,lmx),(imat(i),i=lpr1,lmx)
      np=imat(lmx-1)
      ipage=ipage+1
      go to 20162
20164 nparm=nvars+mrelas
      write(isave) (ibasis(i),i=1,nparm)
      endfile isave
c
c     close file, ipagef, where pages are stored. this is needed so that
c     the pages may be restored at a continuation of splp().
      go to 20317
c ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     procedure (decompose basis matrix)
c++  code for output=yes is active
30004 if (.not.(kprint.ge.2)) go to 20165
      call ivout(mrelas,ibasis,
     *'('' subscripts of basic variables during redecomposition'')',
     *idg)
c++  code for output=no is inactive
c++  end
c
c     set relative pivoting factor for use in la05 () package.
20165 uu=0.1
      call splpdm(
     *mrelas,nvars,lmx,lbm,nredc,info,iopt,
     *ibasis,imat,ibrc,ipr,iwr,ind,ibb,
     *anorm,eps,uu,gg,
     *amat,basmat,csc,wr,
     *singlr,redbas)
      if (.not.(info.lt.0)) go to 20168
      go to 30001
20168 continue
      go to npr004, (20013,20204,20242)
c ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     procedure (classify variables)
c
c     define the classification of the basic variables
c     -1 violates lower bound, 0 feasible, +1 violates upper bound.
c     (this info is stored in primal(nvars+1)-primal(nvars+mrelas))
c     translate variable to its upper bound, if .gt. upper bound
30007 primal(nvars+1)=zero
      call scopy(mrelas,primal(nvars+1),0,primal(nvars+1),1)
      i=1
      n20172=mrelas
      go to 20173
20172 i=i+1
20173 if ((n20172-i).lt.0) go to 20174
      j=ibasis(i)
      if (.not.(ind(j).ne.4)) go to 20176
      if (.not.(rprim(i).lt.zero)) go to 20179
      primal(i+nvars)=-one
      go to 20180
20179 if (.not.(ind(j).eq.3)) go to 10009
      upbnd=bu(j)-bl(j)
      if (j.le.nvars) upbnd=upbnd/csc(j)
      if (.not.(rprim(i).gt.upbnd)) go to 20182
      rprim(i)=rprim(i)-upbnd
      if (.not.(j.le.nvars)) go to 20185
      k=0
20188 call pnnzrs(k,aij,iplace,amat,imat,j)
      if (.not.(k.le.0)) go to 20190
      go to 20189
20190 rhs(k)=rhs(k)-upbnd*aij*csc(j)
      go to 20188
20189 go to 20186
20185 rhs(j-nvars)=rhs(j-nvars)+upbnd
20186 primal(i+nvars)=one
20182 continue
      continue
10009 continue
20180 continue
20176 go to 20172
20174 continue
      go to npr007, (20020,20036)
c cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     procedure (compute error in dual and primal systems)
30005 ntries=1
      go to 20195
20194 ntries=ntries+1
20195 if ((2-ntries).lt.0) go to 20196
      call splpce(
     *mrelas,nvars,lmx,lbm,itlp,itbrc,
     *ibasis,imat,ibrc,ipr,iwr,ind,ibb,
     *erdnrm,eps,tune,gg,
     *amat,basmat,csc,wr,ww,primal,erd,erp,
     *singlr,redbas)
      if (.not.(.not. singlr)) go to 20198
c++  code for output=yes is active
      if (.not.(kprint.ge.3)) go to 20201
      call svout(mrelas,erp,'('' est. error in primal comps.'')',idg)
      call svout(mrelas,erd,'('' est. error in dual comps.'')',idg)
20201 continue
c++  code for output=no is inactive
c++  end
      go to 20193
20198 if (ntries.eq.2) go to 20197
      assign 20204 to npr004
      go to 30004
20204 continue
      go to 20194
20196 continue
20197 nerr=26
      call xermsg ('slatec', 'splpmn',
     +   'in splp, moved to a singular point.  this should not happen.',
     +   nerr, iopt)
      info=-nerr
      go to 30001
20193 continue
      go to npr005, (20018,20154,20243)
c ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     procedure (check feasibility)
c
c     see if nearby feasible point satisfies the constraint
c     equations.
c
c     copy rhs into ww(*), then update ww(*).
30008 call scopy(mrelas,rhs,1,ww,1)
      j=1
      n20206=mrelas
      go to 20207
20206 j=j+1
20207 if ((n20206-j).lt.0) go to 20208
      ibas=ibasis(j)
      xval=rprim(j)
c
c     all variables bounded below have zero as that bound.
      if (ind(ibas).le.3) xval=max(zero,xval)
c
c     if the variable has an upper bound, compute that bound.
      if (.not.(ind(ibas).eq.3)) go to 20210
      upbnd=bu(ibas)-bl(ibas)
      if (ibas.le.nvars) upbnd=upbnd/csc(ibas)
      xval=min(upbnd,xval)
20210 continue
c
c     subtract xval times column vector from right-hand side in ww(*)
      if (.not.(xval.ne.zero)) go to 20213
      if (.not.(ibas.le.nvars)) go to 20216
      i=0
20219 call pnnzrs(i,aij,iplace,amat,imat,ibas)
      if (.not.(i.le.0)) go to 20221
      go to 20220
20221 ww(i)=ww(i)-xval*aij*csc(ibas)
      go to 20219
20220 go to 20217
20216 if (.not.(ind(ibas).eq.2)) go to 20224
      ww(ibas-nvars)=ww(ibas-nvars)-xval
      go to 20225
20224 ww(ibas-nvars)=ww(ibas-nvars)+xval
20225 continue
20217 continue
20213 continue
      go to 20206
c
c   compute norm of difference and check for feasibility.
20208 resnrm=sasum(mrelas,ww,1)
      feas=resnrm.le.tolls*(rprnrm*anorm+rhsnrm)
c
c     try an absolute error test if the relative test fails.
      if(.not. feas)feas=resnrm.le.tolabs
      if (.not.(feas)) go to 20227
      primal(nvars+1)=zero
      call scopy(mrelas,primal(nvars+1),0,primal(nvars+1),1)
20227 continue
      go to npr008, (20024,20032,20040)
c ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     procedure (initialize reduced costs and steepest edge weights)
30014 call spincw(
     *mrelas,nvars,lmx,lbm,npp,jstrt,
     *ibasis,imat,ibrc,ipr,iwr,ind,ibb,
     *costsc,gg,erdnrm,dulnrm,
     *amat,basmat,csc,wr,ww,rz,rg,costs,colnrm,duals,
     *stpedg)
c
      go to npr014, (20135,20246)
c ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     procedure (check and return with excess iterations)
30019 if (.not.(itlp.gt.mxitlp)) go to 20230
      nerr=25
      assign 20233 to npr011
      go to 30011
c++  code for output=yes is active
20233 if (.not.(kprint.ge.1)) go to 20234
      assign 20237 to npr012
      go to 30012
20237 continue
20234 continue
c++  code for output=no is inactive
c++  end
      idum(1)=0
      if(savedt) idum(1)=isave
      write (xern1, '(i8)') mxitlp
      write (xern2, '(i8)') idum(1)
      call xermsg ('slatec', 'splpmn',
     *   'in splp, max iterations = ' // xern1 //
     *   ' taken.  up-to-date results saved on file no. ' // xern2 //
     *   '.  if file no. = 0, no save.', nerr, iopt)
      info=-nerr
      go to 30001
20230 continue
      go to 20155
c ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     procedure (redecompose basis matrix and try again)
30016 if (.not.(.not.redbas)) go to 20239
      assign 20242 to npr004
      go to 30004
20242 assign 20243 to npr005
      go to 30005
20243 assign 20244 to npr006
      go to 30006
20244 assign 20245 to npr013
      go to 30013
20245 assign 20246 to npr014
      go to 30014
20246 continue
c
c     erase non-cycling markers near completion.
20239 i=mrelas+1
      n20247=mrelas+nvars
      go to 20248
20247 i=i+1
20248 if ((n20247-i).lt.0) go to 20249
      ibasis(i)=abs(ibasis(i))
      go to 20247
20249 assign 20251 to npr015
      go to 30015
20251 continue
      go to 20145
c ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     procedure (compute new primal)
c
c     copy rhs into ww(*), solve system.
30006 call scopy(mrelas,rhs,1,ww,1)
      trans = .false.
      call la05bs(basmat,ibrc,lbm,mrelas,ipr,iwr,wr,gg,ww,trans)
      call scopy(mrelas,ww,1,rprim,1)
      rprnrm=sasum(mrelas,rprim,1)
      go to npr006, (20019,20031,20039,20244,20275)
c ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     procedure (compute new duals)
c
c     solve for dual variables. first copy costs into duals(*).
30013 i=1
      n20252=mrelas
      go to 20253
20252 i=i+1
20253 if ((n20252-i).lt.0) go to 20254
      j=ibasis(i)
      if (.not.(j.le.nvars)) go to 20256
      duals(i)=costsc*costs(j)*csc(j) + xlamda*primal(i+nvars)
      go to 20257
20256 duals(i)=xlamda*primal(i+nvars)
20257 continue
      go to 20252
c
20254 trans=.true.
      call la05bs(basmat,ibrc,lbm,mrelas,ipr,iwr,wr,gg,duals,trans)
      dulnrm=sasum(mrelas,duals,1)
      go to npr013, (20134,20245,20267)
c ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     procedure (find variable to enter basis and get search direction)
30015 call splpfe(
     *mrelas,nvars,lmx,lbm,ienter,
     *ibasis,imat,ibrc,ipr,iwr,ind,ibb,
     *erdnrm,eps,gg,dulnrm,dirnrm,
     *amat,basmat,csc,wr,ww,bl,bu,rz,rg,colnrm,duals,
     *found)
      go to npr015, (20141,20251)
c ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     procedure (choose variable to leave basis)
30017 call splpfl(
     *mrelas,nvars,ienter,ileave,
     *ibasis,ind,ibb,
     *theta,dirnrm,rprnrm,
     *csc,ww,bl,bu,erp,rprim,primal,
     *finite,zerolv)
      go to 20149
c ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     procedure (make move and update)
30018 call splpmu(
     *mrelas,nvars,lmx,lbm,nredc,info,ienter,ileave,iopt,npp,jstrt,
     *ibasis,imat,ibrc,ipr,iwr,ind,ibb,
     *anorm,eps,uu,gg,rprnrm,erdnrm,dulnrm,theta,costsc,xlamda,rhsnrm,
     *amat,basmat,csc,wr,rprim,ww,bu,bl,rhs,erd,erp,rz,rg,colnrm,costs,
     *primal,duals,singlr,redbas,zerolv,stpedg)
      if (.not.(info.eq.(-26))) go to 20259
      go to 30001
c++  code for output=yes is active
20259 if (.not.(kprint.ge.2)) go to 20263
      go to 30021
20266 continue
c++  code for output=no is inactive
c++  end
20263 continue
      go to 20153
c ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     procedure(rescale and rearrange variables)
c
c     rescale the dual variables.
30011 assign 20267 to npr013
      go to 30013
20267 if (.not.(costsc.ne.zero)) go to 20268
      i=1
      n20271=mrelas
      go to 20272
20271 i=i+1
20272 if ((n20271-i).lt.0) go to 20273
      duals(i)=duals(i)/costsc
      go to 20271
20273 continue
20268 assign 20275 to npr006
      go to 30006
c
c     reapply column scaling to primal.
20275 i=1
      n20276=mrelas
      go to 20277
20276 i=i+1
20277 if ((n20276-i).lt.0) go to 20278
      j=ibasis(i)
      if (.not.(j.le.nvars)) go to 20280
      scalr=csc(j)
      if(ind(j).eq.2)scalr=-scalr
      rprim(i)=rprim(i)*scalr
20280 go to 20276
c
c     replace translated basic variables into array primal(*)
20278 primal(1)=zero
      call scopy(nvars+mrelas,primal,0,primal,1)
      j=1
      n20283=nvars+mrelas
      go to 20284
20283 j=j+1
20284 if ((n20283-j).lt.0) go to 20285
      ibas=abs(ibasis(j))
      xval=zero
      if (j.le.mrelas) xval=rprim(j)
      if (ind(ibas).eq.1) xval=xval+bl(ibas)
      if (ind(ibas).eq.2) xval=bu(ibas)-xval
      if (.not.(ind(ibas).eq.3)) go to 20287
      if (mod(ibb(ibas),2).eq.0) xval=bu(ibas)-bl(ibas)-xval
      xval = xval+bl(ibas)
20287 primal(ibas)=xval
      go to 20283
c
c     compute duals for independent variables with bounds.
c     other entries are zero.
20285 j=1
      n20290=nvars
      go to 20291
20290 j=j+1
20291 if ((n20290-j).lt.0) go to 20292
      rzj=zero
      if (.not.(ibb(j).gt.zero .and. ind(j).ne.4)) go to 20294
      rzj=costs(j)
      i=0
20297 call pnnzrs(i,aij,iplace,amat,imat,j)
      if (.not.(i.le.0)) go to 20299
      go to 20298
20299 continue
      rzj=rzj-aij*duals(i)
      go to 20297
20298 continue
20294 duals(mrelas+j)=rzj
      go to 20290
20292 continue
      go to npr011, (20051,20233)
c ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c++  code for output=yes is active
c     procedure (print prologue)
30003 idum(1)=mrelas
      call ivout(1,idum,'(''1num. of dependent vars., mrelas'')',idg)
      idum(1)=nvars
      call ivout(1,idum,'('' num. of independent vars., nvars'')',idg)
      call ivout(1,idum,'('' dimension of costs(*)='')',idg)
      idum(1)=nvars+mrelas
      call ivout(1,idum, '('' dimensions of bl(*),bu(*),ind(*)''
     */'' primal(*),duals(*) ='')',idg)
      call ivout(1,idum,'('' dimension of ibasis(*)='')',idg)
      idum(1)=lprg+1
      call ivout(1,idum,'('' dimension of prgopt(*)='')',idg)
      call ivout(0,idum,
     * '('' 1-nvars=independent variable indices.''/
     * '' (nvars+1)-(nvars+mrelas)=dependent variable indices.''/
     * '' constraint indicators are 1-4 and mean'')',idg)
      call ivout(0,idum,
     * '('' 1=variable has only lower bound.''/
     * '' 2=variable has only upper bound.''/
     * '' 3=variable has both bounds.''/
     * '' 4=variable has no bounds, it is free.'')',idg)
      call svout(nvars,costs,'('' array of costs'')',idg)
      call ivout(nvars+mrelas,ind,
     * '('' constraint indicators'')',idg)
      call svout(nvars+mrelas,bl,
     *'('' lower bounds for variables  (ignore unused entries.)'')',idg)
      call svout(nvars+mrelas,bu,
     *'('' upper bounds for variables  (ignore unused entries.)'')',idg)
      if (.not.(kprint.ge.2)) go to 20302
      call ivout(0,idum,
     * '(''0non-basic indices that are negative show variables''
     * '' exchanged at a zero''/'' step length'')',idg)
      call ivout(0,idum,
     * '('' when col. no. leaving=col. no. entering, the entering ''
     * ''variable moved''/'' to its bound.  it remains non-basic.''/
     * '' when col. no. of basis exchanged is negative, the leaving''/
     * '' variable is at its upper bound.'')',idg)
20302 continue
      go to 20011
c ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     procedure (print summary)
30012 idum(1)=info
      call ivout(1,idum,'('' the output value of info is'')',idg)
      if (.not.(minprb)) go to 20305
      call ivout(0,idum,'('' this is a minimization problem.'')',idg)
      go to 20306
20305 call ivout(0,idum,'('' this is a maximization problem.'')',idg)
20306 if (.not.(stpedg)) go to 20308
      call ivout(0,idum,'('' steepest edge pricing was used.'')',idg)
      go to 20309
20308 call ivout(0,idum,'('' minimum reduced cost pricing was used.'')',
     * idg)
20309 rdum(1)=sdot(nvars,costs,1,primal,1)
      call svout(1,rdum,
     * '('' output value of the objective function'')',idg)
      call svout(nvars+mrelas,primal,
     * '('' the output independent and dependent variables'')',idg)
      call svout(mrelas+nvars,duals,
     * '('' the output dual variables'')',idg)
      call ivout(nvars+mrelas,ibasis,
     * '('' variable indices in positions 1-mrelas are basic.'')',idg)
      idum(1)=itlp
      call ivout(1,idum,'('' no. of iterations'')',idg)
      idum(1)=nredc
      call ivout(1,idum,'('' no. of full redecomps'')',idg)
      go to npr012, (20096,20237)
c ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     procedure (print iteration summary)
30021 idum(1)=itlp+1
      call ivout(1,idum,'(''0iteration number'')',idg)
      idum(1)=ibasis(abs(ileave))
      call ivout(1,idum,
     * '('' index of variable entering the basis'')',idg)
      idum(1)=ileave
      call ivout(1,idum,'('' column of the basis exchanged'')',idg)
      idum(1)=ibasis(ienter)
      call ivout(1,idum,
     * '('' index of variable leaving the basis'')',idg)
      rdum(1)=theta
      call svout(1,rdum,'('' length of the exchange step'')',idg)
      if (.not.(kprint.ge.3)) go to 20311
      call svout(mrelas,rprim,'('' basic (internal) primal soln.'')',
     * idg)
      call ivout(nvars+mrelas,ibasis,
     * '('' variable indices in positions 1-mrelas are basic.'')',idg)
      call ivout(nvars+mrelas,ibb,'('' ibb array'')',idg)
      call svout(mrelas,rhs,'('' translated rhs'')',idg)
      call svout(mrelas,duals,'('' basic (internal) dual soln.'')',idg)
20311 continue
      go to 20266
c++  code for output=no is inactive
c++  end
c ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     procedure (return to user)
30001 if (.not.(savedt)) go to 20314
      go to 30020
20317 continue
20314 if(imat(lmx-1).ne.(-1)) call sclosm(ipagef)
c
c     this test is there only to avoid diagnostics on some fortran
c     compilers.
      return
      end
