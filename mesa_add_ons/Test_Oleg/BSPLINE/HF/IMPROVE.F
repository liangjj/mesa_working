*     ------------------------------------------------------------------
*    		I M P R O V E
*     ------------------------------------------------------------------
*
*     If initial estimates were screened hyrogenics, solve equation
*     without exchange using matrix methods, for lowest eigenfunction
*     satisfying orthogonality.
*
*
      SUBROUTINE improve(ind)
      IMPLICIT DOUBLE PRECISION (a-h,o-z)
      PARAMETER (nt=100,kx=8,nwd=10)
*
      INTEGER ouf
      COMMON /INOUT/ iuf,ouf,iscw
*
      COMMON /ORBITS/p(nt,nwd), n(nwd), l(nwd), max(nwd)
      COMMON /ORBSCF/sum(nwd), s(nwd), dpm(nwd), acc(nwd), meth(nwd),
     :               iord(nwd), e(nwd,nwd)
*
      COMMON /PARAM/ z,h,t(nt),fine,mass,ks,ns,nwf,nclosd,ncfg,ib,ic,
     :               id,nscf,no,kmax,lmax,noff
*
      COMMON /SPLINE/gr(nt,kx),grm(nt,kx),grw(nt,kx),
     :        bsp(nt,kx,kx),bspd(nt,kx,kx,2),db1(nt,kx),db2(nt,kx),
     :        sb(nt,kx),r1(nt,kx),rm1(nt,kx),rm2(nt,kx),
     :        h0(nt,kx), rkm(nt,nt,kx,kx), lh, krk
*
      LOGICAL fail,omit,rel,all,trace
      COMMON /TEST/fail,omit,rel,all,trace
*
      DIMENSION h1(nt,kx),yrm(nt,kx),a(nt,nt),b(nt,nt),
     :          eivr(nt),eivi(nt),den(nt),eigvec(nt,nt), w1(nt),
     :          eigval(nt), w2(nt),ind(nwd)
      Equivalence (w1(1), eivr(1)), (w2(1),eivi(1)), (eigval(1),den(1))
*
      do 300 i = 1,nwf
	if (ind(i) .ne. -2 ) then
	  node = n(i)-l(i) -1
   	  call potl(i,rel,yrm)
	  call   hlm(z,l(i),nt,kx,ks,ns,db2,rm1,rm2,h1)
	  do 10 m = 1,ks
	    do 10 ii = 1,ns
	      h1(ii,m) = h1(ii,m) - yrm(ii,m)
10        continue
*
*         .. form full coulomb matrix
*
	  call coulom(nt,kx,ks,ns,sb,h1,a,b)
*         .. add orthogonality constraints to inner orbitals
	  icount = 0
	  do 20 j = 1,i-1
	    if( e(i,j) .ne. 0) then
	      icount = icount + 1
 	      call bxv(nt,ks,ns,sb,p(1,j),w1)
 	      call bxv(nt,ks,ns,h1,p(1,j),w2)
	      do 30 jj = 1,ns
	        do 30 ii = 1,ns
	          a(ii,jj) = a(ii,jj) - w1(ii)*w2(jj)
30            continue
	    end if
20        continue
*
*	  if (icount .eq. 0) then
*	    call rsg(nt,ns-2,a(2,2),b(2,2),eigval,1,eigvec,w1,w2)
*	    keig = ns-node-2
*	  else 
	    call rgg(nt,ns-2,a(2,2),b(2,2),eivr,eivi,den,1,eigvec,ier)
	    do 40 ii = 1,ns-2
	      eivr(ii) = eivr(ii)/den(ii)
40          continue
	  print * , 'Eigenvalues:'
	  print '(6f12.8)', (eivr(ii),ii=1,ns-2)
	  print '(6f12.8)', (eivi(ii),ii=1,ns-2)
	  print '(6f12.8)', (den(ii),ii=1,ns-2)
	    node = node -icount
	    do 50 ii = 1,node
	      keig = idmax(ns-2,eivr,1)
	      eivr(keig) = 0.d0
50          continue
	    keig = idmax(ns-2,eivr,1)
	    print *, 'keig, eigval', keig, eivr(keig)
*	  end if
	  do 200 ii = 1,ns-2
	    p(ii+1,i) = eigvec(ii,keig)
200       continue
	  pn=1.d0/sqrt(quadr(i,i,0))
	  if (p(l(i)+3,i) .lt. 0.d0) pn = -pn
	  do 210 ii = 1,ns
	    p(ii,i) = pn*p(ii,i)
210       continue
	  Print *, 'Eigenvector: ',keig
	  print '(6f12.8)', (p(ii,i),ii=1,ns)
	end if
300   continue
      end
