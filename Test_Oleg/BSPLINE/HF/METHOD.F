*     ------------------------------------------------------------------
*               m e t h o d
*     ------------------------------------------------------------------
*
*       uses m1, m2, or m3 to solve the radial equation. if the input
*   data indicated meth(i) = 3, then this  solution  is  returned  to
*   de.  otherwise,  the routine searches for an acceptable  solution
*   which  is  both  positive  near  the  origin and has the required
*   number  of nodes.
*
      SUBROUTINE method(i,first,ed,pde,m)
      IMPLICIT DOUBLE PRECISION(a-h,o-z)
      PARAMETER(nt=100,kx=8,nwd=10)
*
      COMMON /coeff/coef(100),ijptr(5,5)
      character config*50,el*3,atom*6,term*6
      COMMON /label/config,el(nwd),atom,term
*
      COMMON /ORBITS/p(nt,nwd), n(nwd), l(nwd), max(nwd)
      COMMON /ORBSCF/ sum(nwd), s(nwd), dpm(nwd), acc(nwd), meth(nwd),
     :             iord(nwd), e(nwd,nwd)
*
      COMMON /PARAM/ z,h,t(nt),fine,mass,ks,ns,nwf,nclosd,ncfg,ib,ic,
     :               id,nscf,no,kmax,lmax,noff
      LOGICAL fail,omit,rel,all,trace
      COMMON /test/fail,omit,rel,all,trace
      DOUBLE PRECISION pde(nt),h1(nt,kx),bpv(nt),w(nt),yrm(nt,kx),xv(nt)
      LOGICAL first
      SAVE em,eu,fu,fm
*
*  *****  'first' must be 'true' the first time solve is called for
*  *****  potential and exchange to be computed
*  *****  'eu' is the upper bound of the energy parameter
*  *****  'em' is the minimum value of the energy parameter
*
      fail = .false.
      mit = 0
      if (first) then
        first = .false.
        em = 0.d0
        eu = ((z - dmin1(d5*s(i),d2*s(i)))/n(i))**2
        fu = eu
 	call potl(iel,rel,yrm)
	call xch(iel,rel,xv)
	call hlm(z,l(iel),nt,kx,ks,ns,db2,rm1,rm2,h1)
	do 10 m = 1,ks
	  do 10 ii = 1,ns
	    h1(ii,m) = h1(ii,m) - yrm(ii,m)
10      continue
	call bxv(nt,ks,ns,h1,p(1,iel),bpv)
	do 12 ii = 1,ns
	  bpv(ii) = bpv(ii) - xv(ii)
12      continue
	ed = ddot(ns,p(1,iel),1,bpv,1)
	do 14 j = 1,nwf
	  eij = e(iel,j)
	  if( eij .ne. 0 .and. iej .ne. 1.d-10 .and. iel.ne.j) then
	    call bxvpw(nt,ks,nc,eij,sb,p(1,j),xv,w)
	  end if
14	continue
      end if
17    call solve(i,first,ed,h1,xv,pde,m)
*
*  *****  if mk equals 3, omit the node checking
*
      if (mk .le. 2) then
*
*  *****  count the number of nodes
*
        mn = m
        nc = nodec(mn)
        if (trace) write (6,99) el(i),nc,mn,pde(mn),ed,eu,em
99      format(2x,a3,' nc =',i3,' mn =',i3,' pde(mn) =',
     :         d10.2,' ed =',d10.2,' eu =',d10.2,' em =',d10.2)
*
*  *****  if node count is off by no more than 1 and deltae is still
*  *****  quite large, apply the deltae correction
*
        if (iabs(nc-node) .eq. 0) then
*
*         .. we have the correct number of nodes
          if (pde(mn) .lt. 0.d0 ) then
            do 52 j = 1,no
              pde(j) = - pde(j)
52          continue
	  end if
          return
	else if (nc - node .lt. 0) then
*
*  *****  the solution has too few nodes
*
          if (pde(mn) .gt. 0.d0) then
            del = 1.d0 - ed/eu
            eu = ed
            if ( del .lt. .05d0) fu = fu*((l(i)+1+nc)/fn)**2.5
            if (del  .ge. .05d0) fu = ed*((l(i)+1+nc)/fn)**2.5
            if (fu .lt. em) fu = d5*(eu + em)
            if (dabs(fu - ed) .lt. 0.001d0) then
              write (6,28) mk,el(i),nc,ed,em,eu
	      fail = .true.
	      return
	    else
	      mit = mit+1
              ed = fu
            end if
	  else
*
*           ..try a new value of ed which must lie within the upper and
*             lower bound
*
            edp = ed
            ed = ed*((l(i)+1+nc)/fn)**2.5
            if (ed .ge. eu ) ed = d5*(eu + edp)
            if (ed .le. em ) ed = d5*(em + edp)
33          mit = mit + 1
	  end if
	else 
*
*         ..the solution has too many nodes
*
          if (pde(mn) .gt. 0.d0) then
            del = 1.d0 - em/ed
            em = ed
            if (del .lt. 0.05d0) fm = fm*((l(i)+1+nc)/fn)**2.5
            if (del .ge. 0.05d0) fm = ed*((l(i)+1+nc)/fn)**2.5
            if (fm .gt. eu) fm = d5*(eu + em)
            if (dabs(fm - ed) .lt. 0.001d0) then
               write (6,28) mk,el(i),nc,ed,em,eu
	       fail = .true.
	       return
	    end if
            ed = fm
          end if
	end if
      end if
      if ( eu .le. em ) write(6,30) em,eu,ed
30    format(6x,48hwarning: difficulty with node counting procedure/
     :       6x,42hlower bound on ed greater than upper bound/
     :       6x,5hel = ,f10.6,7h  eu = ,f10.6,7h  ed = ,f10.6)
      if ( mit .gt. 3*n(i) .or. eu-em .lt. fn**(-3)) then
         write (6,28) mk,el(i),nc,ed,em,eu
	 fail = .true.
	 return
      else
        go to 17
      end if
*
28    format(10x,6hmethod,i2,38h unable to solve equation for electron,
     :   a3/10x,5hnc = ,i3,3x,5hed = ,f10.6,3x,5hel = ,
     :   f10.6,3x,5heu = ,f10.6)
      end
