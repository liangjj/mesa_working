*====================================================================
	subroutine facdyk(nt,kx,ktx,ks,k,nint,t,bsp,bspd,db2,rm2,
     :                    ipvtd,dyk)
*====================================================================
*
*      	Sets up and factorizes the matrix
*
*           DD(i,j)-k*(k+1)*rm2
*
*       for solving the following differential equation by the
*       spline-Galerkin method:
*   
*	    (d^2-k(k+1)/r^2)yk(r) = g(r), 
*
*       with boundary conditionas yk(0) = 0, yk(rmax) = const,
*       if k=0, and  dy/dr + (k/r)y = 0 otherwise.
*
*	SUBROUTINES called:
* 	    dgbfa (or sgbfa)
*
*	This version dated 6/90
*
*--------------------------------------------------------------------
*
        implicit double precision(a-h,o-z)
        dimension bsp(nt,kx,kx),bspd(nt,kx,kx,2),db2(nt,*),rm2(nt,kx),
     :            dyk(ktx,nt),t(*),ipvtd(nt)
*
*       on entry
*       --------
*
*       nt   the leading dimension of arrays. must be greater than 
*            the number of elements in the knot sequence.
*       kx   maximum order of the spline.
*       ktx  the leading dimension of dyk. must be greater than 
*            3*kx-2.
*       ks   the order of the spline.
*       k    the order of y. 
*       nint the number of intervals starting at t(k).
*	t    the grid knots for B-spline
*       bsp  array of values of the B-splines at Gaussian points
*       bspd array of B', and B" at Gaussian points
*       db2   two dimensional array, stores integral <B_i,B"_j>.
*       rm2  two dimensional array, stores integral <B_i,(1/r^2)B_j>.
*
*       on exit
*       -------
* 
*       dyk  factorized array for the differential with operator
*            (d^2-k(k+1)/r^2). dyk is banded with
*            width 3*ks-2. The first ks-1 rows are zeros, and the
*            following k rows are the elements above the diagnal,
*            the last ks-1 rows are the elements below the diagnal of
*            the original arrays.
*       iptvd  an integer vector of pivot indices.
*
*--------------------------------------------------------------------
*
	ns = nint+ks-1
	fkk = -(1.d0+k)*k
*
*  clear the dyk array
*
  	do 10 i=1,ks-1
	  do 10 j=1,ns
            dyk(i,j)=0.d0
 10	continue 
*
*  set up dyk
*
*       ... lower portion
	do 20 j=1,ks
   	  do 20 i=ks-j+1,ns
	    dyk(3*ks-1-j,i-ks+j)=db2(i,j)+fkk*rm2(i,j)
 20	continue            
*      ... upper portion
	do 30 j=2,ks
	  do 30 i=1,ns-j+1
	    dyk(2*ks-j,i+j-1)=db2(i+j-1,ks-j+1)+fkk*rm2(i+j-1,ks-j+1)
 30	continue
*
*       .. correct the upper portion of the matrix for the asymmetry
*          of d^2/dr^2
	j = 2*ks-2
	dyk(j,ns) = dyk(j,ns) + db2(1,1)-db2(ns,ks-1)
*
* apply the zero boundary condition at the origin and modify the
* last equation for the boundary condition at rmax
*
	do 40 i=1,ks
      	    j=2*ks-i
	    dyk(j,i)=0.d0
	    j=3*ks-1-i
	    dyk(j,ns-ks+i)=0.d0
 40	continue
 	dyk(2*ks-1,1)=1.d0 
*
* apply the boundary condition at the end of the range:  if k <> 0, 
* we will use an explicit condition:
*             dy/dr + (k/r) y = 0 at rmax.
* otherwise it will be a constant determined from integration.
*
	if (k .ne. 0) then
	  nn = ns-ks+2
          dyk(2*ks,ns-1)=bspd(nn,1,ks-1,1)
 	  dyk(2*ks-1,ns)=bspd(nn,1,ks,1) +bsp(nn,1,ks)*k/t(ns+1)
	else
          dyk(2*ks-1,ns) = 1.d0
	end if
*
* factorize dyk
*
	call dgbfa(dyk,ktx,ns,ks-1,ks-1,ipvtd,ier)
	if (ier .ne. 0) stop 'dyk is non_invertible'
* 
*	do 99 m = 1,3*ks-2
*	  print '(6f12.8)', (dyk(m,ii),ii=1,ns)
*99      continue
	return
	end	
