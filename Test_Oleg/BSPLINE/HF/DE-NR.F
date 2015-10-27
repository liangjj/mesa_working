!=======================================================================
  SUBROUTINE de(iel) 
!=======================================================================
!       This routine controls the improvement of the radial functions
!   with the same orbital quantum number.  
!
!   The Newton Raphson iteration is used to solve the
!   non-linear systems of equations.
!
!----------------------------------------------------------------------
!
    USE spline_param
    USE spline_galerkin
    USE atomic_case
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: iel


    REAL(KIND=8), DIMENSION(ns,ks) :: h_one, yrm
    REAL(KIND=8), DIMENSION(ns)    :: bpv

    print *, 'DE:  nwf',nwf, 'iel',iel
    call potl(iel,rel,yrm)
    call xch(iel,rel,xv)
    call hlm(z,l(iel))  
    h_one = hl - yrm
    call bxv(h_one,p(1,iel),bpv)
    bpv = bpv - xv
    eii = ddot(ns,p(1,iel),1,bpv,1)
    print *,'eii', eii
    h_one = h_one - eii*sb
    do j = 1,nwf
      eij = e(iel,j)
      if( eij N> 0 .and. iej <> 1.d-10 .and. iel <> j) then
	call bxvpw(nt,ks,nc,eij,sb,p(1,j),xv,w)
      end if
    end do
  
    !  clear the bandde array
    bandde = 0
*
*  set up bandde
*
*       ... lower portion
	do 30 j=1,ks
   	  do 30 i=ks-j+1,ns
	    bandde(3*ks-1-j,i-ks+j)=h1(i,j)
 30	continue            
*      ... upper portion
	do 40 j=2,ks
	  do 40 i=1,ns-j+1
	    bandde(2*ks-j,i+j-1)= h1(i+j-1,ks-j+1)
 40	continue
*
* apply the zero boundary condition at the origin and at rmax
*
	do 50 i=1,ks
*           .. diagonals
      	    j=2*ks-i
	    bandde(j,i)=0.d0
	    j=3*ks-1-i
	    bandde(j,ns-ks+i)=0.d0
*           .. columns
	    bandde(2*ks -2 +i,1) = 0.d0
	    bandde(ks-1+i,ns) = 0.d0
 50	continue
 	bandde(2*ks-1,1)=1.d0 
        bandde(2*ks-1,ns) = 1.d0

*	print *, ' Band DE'
*	do 99 m = 1,3*ks-2
*	  print '(6f12.8)', (bandde(m,ii),ii=1,ns)
*99      continue
*	print *, ' h1'
*	do 199 m = 1,ks
*	  print '(6f12.8)', (h1(ii,m),ii=1,ns)
*199     continue
*
* factorize bandd
*
	call dgbfa(bandde,ktx,ns,ks-1,ks-1,ipvtd,ier)
	if (ier .ne. 0) then
          print *, 'bandDE is non_invertible'
	  print *, ' Ier = ',ier
	  stop
	end if
	call dgbsl(bandDE,ktx,ns,ks-1,ks-1,ipvtd,xv,0)
* 
	return
	end	
