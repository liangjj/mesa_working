!============================================================================
!  PROGRAM: Parallel_FEDVR_VTV_TDCC_Imaginary_Time.f90
!
!	Copyright Reserved by Suxing Hu @ LANL
!		Written on July 15, 2005
!============================================================================
!	Using MPI to parallelize in space !!!
!============================================================================
!
!	The program is to use imaginary-time relaxation to get the GS of a
!	two-electron quantum systems in 3D (TDCC method)
!
!==========================================================================================================
!==========================================================================================================
!    Driver for multi-dimensional Finite-Element Discrete-Variable-Representation (FEDVR) codes
!==========================================================================================================
!
!    By Suxing Hu at LANL, references from Barry and Nicolai's DVR codes
!
!    Written on:  07/14/2005
!    Revised on:  mm/dd/yy
!==========================================================================================================
!==========================================================================================================

!==========================================================================================================
!==========================================================================================================
!==========================================================================================================
!==========================================================================================================
!==========================================================================================================
!   MAIN PROGRAM STARTING FROM HERE!
!==========================================================================================================
!==========================================================================================================
!==========================================================================================================
!==========================================================================================================
!==========================================================================================================
Program Parallel_FEDVR_All
  USE nrtype
  USE globalmod
  USE pfedvrmod
  USE dvrmod
  USE anglib
  USE keops
  USE integrals
  USE tdcc_itsubroutines
  USE mpi
  !USE la_precision
  USE f95_lapack
  Implicit NONE

  REAL(DP), dimension(:), allocatable :: lbound, rbound, wavenorms, tmpwavenorms, dvr_weights
  REAL(DP), dimension(:,:,:), allocatable :: V_r1r2lambda
  REAL(DP), dimension(:,:), allocatable :: Ttemp, Tinv, Tident
  INTEGER(I4B), dimension(:), allocatable :: ipivot


  REAL(DP) :: pkx, pky, pkz, pnorm, time, x0

  REAL(DP) :: ttx1,ttx2,starttime,endtime,slope1,slope2

  INTEGER(I4B) :: maxregion, maxfunc, maxpoint, npro, info
  INTEGER(I4B) :: tag, tagp, lambda, idim, i1start, i1end, ireg, NR, NReff, ntmp
  INTEGER(I4B) :: i, j, k, m, n, l

  REAL(DP), dimension(:), allocatable :: work
  INTEGER(I4B) :: lwork

  !------Various real working variables--------------

  REAL(DP) :: temp,temp1,temp2,temp3

  REAL(DP) :: ctemp,ctemp1,ctemp2,ctemp3,ctemp4,ctemp5

  !===========for parallel code using=============================

  INTEGER(I4B) ::	 ntemp3, ntemp4, ntemp5, ntemp6, nmax, LLmax,l1max,l2max

  INTEGER(I4B) ::	 MYID, IERR, status(mpi_status_size)

  !==========================================================================================================
  !==========================================================================================================

  CALL MPI_INIT(IERR)

  CALL MPI_COMM_SIZE(MPI_COMM_WORLD,ncount,IERR)

  CALL MPI_COMM_RANK(MPI_COMM_WORLD,MYID,IERR)

  nrcount=ncount-1

  write(6,*)'MYID=',MYID,' Come here after initializ...'

  ALLOCATE (ntarget(nrcount), nst1(nrcount), nst2(nrcount) )

  !==========================================================================================================
  !==========================================================================================================

803 FORMAT       (9999(G20.10E3))
801 FORMAT (I10,  9999(G20.6E3))
802 FORMAT (G13.6,9999(G20.10E3))

  !-----------------------------------------------------------------
  !-----------------------------------------------------------------

  !***************MASTER CPU DOES THE FOLLOWING WORKS***************


  IF(MYID.EQ.0)THEN

     open(9,file='par.dat',status='unknown')
     open(10,file='norm.dat',status='unknown')
     open(11,file='avecr1.dat',status='unknown')
     open(12,file='avecr2.dat',status='unknown')
     open(14,file='avecr.dat',status='unknown')
     open(111,file='avepr1.dat',status='unknown')
     open(112,file='avepr2.dat',status='unknown')
     open(114,file='avepr.dat',status='unknown')
     open(15,file='pop1.dat',status='unknown')
     open(16,file='fieldx.dat',status='unknown')

     open(66,file='eige.dat',status='unknown')

  ENDIF

  !-----------------------------------------------------------------
  !-----------------------------------------------------------------

  !-----------------------------------------------------------------
  !-----------------------------------------------------------------
  ! Parameters are read from file (parameters.inp & regions.inp)
  !-----------------------------------------------------------------
  !-----------------------------------------------------------------

  IF(MYID.EQ.0)THEN

     open(1,file='parameters.inp',status='OLD')
     read(1,*)ntmp ! this used to be ndim, but this was set as a parameter in globalmod now
     read(1,*)intensity_x, wavelx, durationx, phasex
     read(1,*)Z1
     read(1,*)LLmax,nnp,l1max,l2max
     read(1,*)nabx, naby, VRi, alphai
     read(1,*)tstart, tstop, ntim, ntrec, ntprob
     read(1,*)ncpux, ncpuy
     close(1)

     !--------------------------------------------------------

     if( ncount.ne.(ncpux*ncpuy) )then
        WRITE(6,*)'ncount=',ncount
        WRITE(6,*)'ncpux=',ncpux
        WRITE(6,*)'ncpuy=',ncpuy
        WRITE(6,*)'The total number of CPUs can NOT be decomposed to ncpux*ncpuy!!!'
        STOP
     endif

  ENDIF

  !--------------------------------------------------------
  !--------------------------------------------------------

  CALL MPI_BARRIER(MPI_COMM_WORLD, IERR)

  CALL MPI_BCAST(ncpux,1,MPI_INTEGER,0,MPI_COMM_WORLD, IERR)
  CALL MPI_BCAST(ncpuy,1,MPI_INTEGER,0,MPI_COMM_WORLD, IERR)

  CALL MPI_BARRIER(MPI_COMM_WORLD, IERR)

  !==========================================================================================================
  !	Constructing the 2D Cartesian Communicator for "2D decomposition"
  !==========================================================================================================

  cperiods(1) = 0 ! 0 should be .false.
  cperiods(2) = 0 ! 0 should be .false.
  creorder =    1 ! 1 should be .true.
  cdims(1)=ncpux
  cdims(2)=ncpuy

  CALL MPI_Cart_Create (MPI_COMM_WORLD, cndim, cdims, cperiods, creorder, comm2d, IERR)

  nbrup=MPI_PROC_NULL
  nbrdown=MPI_PROC_NULL
  nbrleft=MPI_PROC_NULL
  nbrright=MPI_PROC_NULL

  CALL MPI_Cart_Shift (comm2d, 0, 1, nbrup, nbrdown, IERR)
  CALL MPI_Cart_Shift (comm2d, 1, 1, nbrleft, nbrright, IERR)

  CALL MPI_Cart_Coords (comm2d, MYID, cndim, coords, IERR)

804 FORMAT (A,I4.4,A,9999(I5.4,1x))
  write(6,804)'MYID = ',MYID,'   coords  = ',coords
  write(6,804)'MYID = ',MYID,'   nbrup   = ',nbrup
  write(6,804)'MYID = ',MYID,'   nbrdown = ',nbrdown
  write(6,804)'MYID = ',MYID,'   nbrleft = ',nbrleft
  write(6,804)'MYID = ',MYID,'   nbrright= ',nbrright
  write(6,*)

  !==========================================================================================================

  CALL MPI_BARRIER(comm2d, IERR)

  CALL MPI_BCAST(ndim,1,MPI_INTEGER,0,comm2d, IERR)
  CALL MPI_BCAST(nabx,1,MPI_INTEGER,0,comm2d, IERR)
  CALL MPI_BCAST(naby,1,MPI_INTEGER,0,comm2d, IERR)
  CALL MPI_BCAST(ntim,1,MPI_INTEGER,0,comm2d, IERR)
  CALL MPI_BCAST(ntrec,1,MPI_INTEGER,0,comm2d, IERR)
  CALL MPI_BCAST(ntprob,1,MPI_INTEGER,0,comm2d, IERR)
  CALL MPI_BCAST(LLmax,1,MPI_INTEGER,0,comm2d, IERR)
  CALL MPI_BCAST(nnp,1,MPI_INTEGER,0,comm2d, IERR)
  CALL MPI_BCAST(l1max,1,MPI_INTEGER,0,comm2d, IERR)
  CALL MPI_BCAST(l2max,1,MPI_INTEGER,0,comm2d, IERR)

  CALL MPI_BCAST(intensity_x,1,MPI_DOUBLE_PRECISION,0,comm2d, IERR)
  CALL MPI_BCAST(wavelx,1,MPI_DOUBLE_PRECISION,0,comm2d, IERR)
  CALL MPI_BCAST(durationx,1,MPI_DOUBLE_PRECISION,0,comm2d, IERR)
  CALL MPI_BCAST(phasex,1,MPI_DOUBLE_PRECISION,0,comm2d, IERR)
  CALL MPI_BCAST(tstart,1,MPI_DOUBLE_PRECISION,0,comm2d, IERR)
  CALL MPI_BCAST(tstop,1,MPI_DOUBLE_PRECISION,0,comm2d, IERR)

  CALL MPI_BCAST(Z1,1,MPI_DOUBLE_PRECISION,0,comm2d, IERR)
  CALL MPI_BCAST(VRi,1,MPI_DOUBLE_PRECISION,0,comm2d, IERR)
  CALL MPI_BCAST(alphai,1,MPI_DOUBLE_PRECISION,0,comm2d, IERR)

  CALL MPI_BARRIER(comm2d, IERR)

  !--------------------------------------------------------
  ! Reading setups for the finite-element regions
  !--------------------------------------------------------

  ALLOCATE(  num_reg(ndim) , ngridstart(ndim), ngridstop(ndim), ngstart(ndim), ngstop(ndim) )	
  ALLOCATE(  bdstart(ndim) , bdstop(ndim), bdslope(ndim), numfunc(ndim) )	


  IF(MYID.EQ.0)THEN

     !================================================================
     open(2,file='regions.inp',status='OLD')

     read(2,*)(num_reg(j), j=1, ndim)
     read(2,*)(numfunc(j), j=1, ndim)
     DO j=1, ndim
        read(2,*)bdstart(j), bdstop(j), bdslope(j)
     END DO

     close(2)
     !================================================================

     nmax=1
     DO j=1, ndim
        if(nmax.lt.num_reg(j))nmax=num_reg(j)
     END DO


     !----------testing if the num_reg(n)/ncpuxy=even--------------------

     IF(mod(num_reg(1),ncpux).eq.0)THEN
        nelemx=num_reg(1)/ncpux
        if(mod(nelemx,2).ne.0)then
           write(9,*)
           write(9,*)'Num_reg(1)/ncpux=',num_reg(1)/dfloat(ncpux)
           write(9,*)'The ratio of number of elements along the x-dimension to the number of CPUs along that dimension must be EVEN integer!'
           STOP
        endif
     ELSE
        write(9,*)
        write(9,*)'Num_reg(1)/ncpux=',num_reg(1)/dfloat(ncpux)
        write(9,*)'This ratio must be EVEN integer!'
        STOP
     ENDIF

     IF(mod(num_reg(2),ncpuy).eq.0)THEN
        nelemy=num_reg(2)/ncpuy
        if(mod(nelemy,2).ne.0)then
           write(9,*)
           write(9,*)'Num_reg(2)/ncpuy=',num_reg(2)/dfloat(ncpuy)
           write(9,*)'The ratio of number of elements along the y-dimension to the number of CPUs along that dimension must be EVEN integer!'
           STOP
        endif
     ELSE
        write(9,*)
        write(9,*)'Num_reg(2)/ncpuy=',num_reg(2)/dfloat(ncpuy)
        write(9,*)'This ratio must be EVEN integer!'
        STOP
     ENDIF

  ENDIF

  !-------------------------------------------------------------
  !-------------------------------------------------------------

  CALL MPI_BARRIER(comm2d, IERR)

  CALL MPI_BCAST(nmax,1,MPI_INTEGER,0,comm2d, IERR)
  CALL MPI_BCAST(nelemx,1,MPI_INTEGER,0,comm2d, IERR)
  CALL MPI_BCAST(nelemy,1,MPI_INTEGER,0,comm2d, IERR)
  CALL MPI_BCAST(num_reg,ndim,MPI_INTEGER,0,comm2d, IERR)
  CALL MPI_BCAST(numfunc,ndim,MPI_INTEGER,0,comm2d, IERR)
  CALL MPI_BCAST(bdstart,ndim,MPI_DOUBLE_PRECISION,0,comm2d, IERR)
  CALL MPI_BCAST(bdstop,ndim,MPI_DOUBLE_PRECISION,0,comm2d, IERR)
  CALL MPI_BCAST(bdslope,ndim,MPI_DOUBLE_PRECISION,0,comm2d, IERR)

  CALL MPI_BARRIER(comm2d, IERR)

  !-------------------------------------------------------------
  !-------------------------------------------------------------

  ALLOCATE( ntot(ndim), regstart(ndim), regstop(ndim) )	
  ALLOCATE( num_fun(ndim,nmax) )
  ALLOCATE( bounds(ndim,nmax,1:2) )

  ALLOCATE( mat_reg(ndim,nmax) )
  ALLOCATE( nindex(ndim,nmax,500) )

  bounds(:,:,:)=0.d0
  num_fun(:,:)=0

  !*******************deciding how many regions in each dimension for different CPU********************

  regstart(1)=coords(1)*nelemx+1
  regstop(1)=(coords(1)+1)*nelemx

  regstart(2)=coords(2)*nelemy+1
  regstop(2)=(coords(2)+1)*nelemy

  !*************************************************************

  !-------------------------------------------------------------
  ! Generating the regions
  !-------------------------------------------------------------

  DO j=1, ndim

     temp=dabs(bdstop(j)-bdstart(j))

     i=1
     bounds(j,i,1)=bdstart(j)
     bounds(j,i,2)=bdstart(j)+temp
     num_fun(j,i)=numfunc(j)

     do i=2, num_reg(j)
        bounds(j,i,1)=bounds(j,i-1,2)
        bounds(j,i,2)=bounds(j,i,1)+temp*dexp(bdslope(j)* &
             dfloat(i)/dfloat(num_reg(j)))
        num_fun(j,i)=numfunc(j)
     end do

  END DO


  !--------------------------------------------------------
  ! Check if the regions are continuous
  !--------------------------------------------------------

  IF(MYID.eq.0)THEN

     DO j=1, ndim
        do i=2,num_reg(j)-1
           if(dabs(bounds(j,i,1) - bounds(j,i-1,2)).gt.1.d-10) then
              write(6,*)
              write(6,*)'  The ',j,'th dimension has problem:'
              write(6,*)'  Gap between element',i-1,'and',i,' program halted'
              write(6,*)' bounds(j,i-1,2)=',bounds(j,i-1,2)
              write(6,*)' bounds(j,i,1)=',bounds(j,i,1)
              write(6,*)
              stop
           end if
        end do
     END DO

  ENDIF

  CALL MPI_BARRIER(comm2d, IERR)

  !--------------------------------------------------------
  ! Find total number of grid points
  ! this is the sum of number of functions in each region
  ! minus the number of overlap points
  !--------------------------------------------------------

  do j=1, ndim
     ntot(j)=sum(num_fun(j,:))-(num_reg(j)-1)
     !-----------throw out the first point----------------
     !     ntot(j)=ntot(j)-1
  end do
  ntemp1 = maxval(ntot(:))
  !-------------------------------------------------------------
  !-------------------------------------------------------------

  ALLOCATE ( point_order(ndim,ntemp1), factor(ndim,ntemp1) )

  !-------------------------------------------------------------
  !   Allocating matrix array for calculations in each element!
  !-------------------------------------------------------------

  !*********first checking if the number of basis functions are the same for each dimension******


  do i=2, num_reg(1)
     if( num_fun(1,i).ne.num_fun(1,i-1) )then
        write(6,*)'i=',i
        write(6,*)'num_fun(1,i)=',num_fun(1,i),' num_fun(1,i-1)=',num_fun(1,i-1)
        write(6,*)'!!!! WE NEED KEEP THE NUMBER OF BASIS SAME FOR EACH ELEMENT IN THE SAME DIMENSION!!!!'
        STOP
     endif
  end do

  do i=2, num_reg(2)
     if( num_fun(2,i).ne.num_fun(2,i-1) )then
        write(6,*)'i=',i
        write(6,*)'num_fun(2,i)=',num_fun(2,i),' num_fun(2,i-1)=',num_fun(2,i-1)
        write(6,*)'!!!! WE NEED KEEP THE NUMBER OF BASIS SAME FOR EACH ELEMENT IN THE SAME DIMENSION!!!!'
        STOP
     endif
  end do


  !--------------Allocating matrice-------------------------------------------------


  ALLOCATE ( prob2d(num_fun(1,1),num_fun(2,1)),   probx2d(num_fun(1,1),num_fun(2,1)), &
       proby2d(num_fun(1,1),num_fun(2,1)),  cprob2d(num_fun(1,1),num_fun(2,1)), &
       cprobpx2d(num_fun(1,1),num_fun(2,1)), cprobpy2d(num_fun(1,1),num_fun(2,1)) )


  !-------------------------------------------------------------
  !   Assuming equal number of basis in each element!
  !-------------------------------------------------------------


  do i=1, ndim
     if(i.eq.1)then
        ntemp1=coords(1)*nelemx*(num_fun(i,1)-1)+1
        ngridstart(i)=ntemp1
        ngstart(i)=ntemp1
        do j=regstart(i), regstop(i)
           ntemp1=ntemp1+(num_fun(i,j)-1)
        end do
        ngridstop(i)=ntemp1	
        ngstop(i)=ntemp1-1
     else if(i.eq.2)then
        ntemp2=coords(2)*nelemy*(num_fun(i,1)-1)+1
        ngridstart(i)=ntemp2
        ngstart(i)=ntemp2
        do j=regstart(i), regstop(i)
           ntemp2=ntemp2+(num_fun(i,j)-1)
        end do
        ngridstop(i)=ntemp2	
        ngstop(i)=ntemp2-1	
     endif

     write(6,*)'i=',i,' MYID=',MYID,' ngridstart(i)=', ngridstart(i)
     write(6,*)'i=',i,' MYID=',MYID,' ngridstop(i)=', ngridstop(i)
     write(6,*)'i=',i,' MYID=',MYID,' regstart(i)=', regstart(i),' regstop(i)=', regstop(i)
     write(6,*)

  end do


  !-------------------------------------------------------------
  !-----------throw out the first point for both r1 and r2------
  !-------------------------------------------------------------
  !***********This adjustment is very important for the energy calculation***********

  if(coords(1).eq.0)then
     ngstart(1)=ngstart(1)+1
  endif

  if(coords(2).eq.0)then
     ngstart(2)=ngstart(2)+1
  endif


  !	if(MYID.eq.0)then
  !	  ngridstop(1)=ngridstop(1)-1
  !	  ngstop(1)=ngstop(1)-1

  !	  ngridstop(2)=ngridstop(2)-1
  !	  ngstop(2)=ngstop(2)-1
  !	else
  !	  ngridstart(1)=ngridstart(1)-1
  !	  ngridstop(1)=ngridstop(1)-1
  !	  ngstart(1)=ngstart(1)-1
  !	  ngstop(1)=ngstop(1)-1
  !
  !	  ngridstart(2)=ngridstart(2)-1
  !	  ngridstop(2)=ngridstop(2)-1
  !	  ngstart(2)=ngstart(2)-1
  !	  ngstop(2)=ngstop(2)-1
  !	endif
  !	
  !------------------------------------------------------------------------


  CALL MPI_BARRIER(comm2d, IERR)


  !-------------------------------------------------------------
  !-------------------------------------------------------------



  !----------point_order stores the #-order of point in that element----------

  !-------------------------------------------------------------
  !	In order to keep the number of points in each CPU be
  !	equal, we do not throw this origin point, but just
  !	impose the boundary in the propagator!!
  !-------------------------------------------------------------

  r0BC=.false.


  ! When solving radial equation we throw away the first point
  !   if((symmetry == 'spherical').or.(symmetry == 'cylindrical')) then
  !      ntot=ntot-1
  !      r0BC=.true.
  !   else
  !      mu=0.d0
  !      r0BC=.false.
  !   end if

  ! Centrifugal barrier constant
  !   if(symmetry == 'spherical') then
  !      alpha=mu*(mu+1.d0)
  !   elseif(symmetry == 'cylindrical') then
  !      alpha=mu**2-1.d0/4.d0
  !   else
  !      alpha=0.d0
  !   end if


  !------------------------------------------------------------------
  !   Generating the parameters
  !------------------------------------------------------------------

  efieldx=0.53309178d-8*dsqrt(intensity_x)
  omegax=2.d0*pi*3.d0*2.418884d0/wavelx
  TCx=twopi/omegax
  durationx=durationx/convfs
  phasex=phasex*pi/180.0
  tstart=tstart/convfs
  tstop=tstop/convfs
  dt=(tstop-tstart)/dfloat(ntim)

  !------------------------------------------------------------------

  !------------------------------------------------
  ! Outputing parameters to file
  !------------------------------------------------

  IF(MYID.EQ.0)THEN

     write(9,*)'-----------------------------------------------------------'
     write(9,*)'   The code is to relax an initial trial wave packet to the '
     write(9,*)'   GS of the two-electron system'
     write(9,*)'-----------------------------------------------------------'
     write(9,*)


     write(9,*)'==============================================================='
     write(9,*)'==Total number of CPUs =',ncount
     write(9,*)'==The number of CPUs along x-dimensional decomposition =',ncpux
     write(9,*)'==The number of CPUs along y-dimensional decomposition =',ncpuy
     write(9,*)'==============================================================='


     write(9,*)


     write(9,*)'=================Parameters for the FEDVR program=============='
     write(9,*)
     write(9,*)'The x-axis FCP laser intensity I=',intensity_x,'  W/cm^2'
     write(9,*)'The x-axis FCP laser wavelength wavelx=',wavelx,' (nm)'
     write(9,*)'The x-axis FCP durationx =',durationx*convfs,' (fs)'
     write(9,*)'The x-axis phase between the carrier and the FCP envelope phasex=',phasex,' (degree)'
     write(9,*)
     write(9,*)
     write(9,*)'The core charges Z1=',Z1
     write(9,*)
     write(9,*)'The start time tstart==',tstart*convfs, '  (fs)'
     write(9,*)'The stop time tstop==',tstop*convfs, '  (fs)'
     write(9,*)
     write(9,*)
     write(9,*)'The absorption grids nabr1 =',nabx
     write(9,*)'The absorption grids nabr2 =',naby
     write(9,*)'The absorption strength VRi=',VRi
     write(9,*)'The absorption parameter alphai=',alphai
     write(9,*)
     write(9,*)
     write(9,*)'The order of multipole expansion np= 0 --- ',nnp
     write(9,*)'The total number of PWs for expansion Lmax=',0,' - ',LLmax
     write(9,*)'The number of PWs for expansion l1max= 0 - ',l1max
     write(9,*)'The number of PWs for expansion l2max= 0 - ',l2max
     write(9,*)

     write(9,*)
     write(9,*)'The total time steps ntim=',ntim
     write(9,*)'The recording frequency ntrec=',ntrec
     write(9,*)'The probability recording frequency ntprob=',ntprob
     write(9,*)
     write(9,*)'-----------------------------------------------------'
     write(9,*)'The exponential increasing rate of element size:'
     write(9,*)'-----------------------------------------------------'
     do i=1, ndim
        write(9,*)' The ',i,'-th dimension : bdslope(i)=', bdslope(i)
     end do
     write(9,*)'-----------------------------------------------------'
     write(9,*)

     !------------------------------------------------------------------

     write(9,*)
     write(9,*)'------------------------------------------------'
     write(9,*)'The time step dt=',dt
     write(9,*)'------------------------------------------------'
     write(9,*)
     write(9,*)

     !------------------------------------------------------------------

     write(9,*)'=================Regions of FEDVR being used=============='
     write(9,*)
     write(9,*)'This is a ',ndim,'-dimensional radial problem + partial waves '
     write(9,*)

     DO j=1, ndim
	write(9,*)'============================================================================='
	write(9,*)'     		For the ',j,'th dimension:'
	write(9,*)'============================================================================='
	write(9,*)
	write(9,*)'The total number of elements:: num_reg=',num_reg(j)
	write(9,*)'The total DVR-grid points:: ntot(j)=',ntot(j)
	write(9,*)
	write(9,*)'----------------------------------------------------------------------------'
	write(9,*)'# of region       left-boundary        right-boundary     # of basis used'
	write(9,*)'----------------------------------------------------------------------------'
        do i=1, num_reg(j)
           write(9,888)i,bounds(j,i,1),bounds(j,i,2),num_fun(j,i)
        end do
	write(9,*)
	write(9,*)
     END DO

888  format (I15, G18.6, G18.6, I20)

     call flush(5)


  ENDIF


  CALL MPI_BARRIER(comm2d, IERR)

  !-----------------------------------------------------------------
  !-----------------------------------------------------------------

  IF(MYID.eq.0)THEN

     ALLOCATE ( lbound(ndim), rbound(ndim) )

     !-----------------------------------------------------------------
     ! Endpoint of interval: i.e., space for each dimension
     !-----------------------------------------------------------------

     DO j=1, ndim
        lbound(j)=bounds(j,1,1)
        rbound(j)=bounds(j,num_reg(j),2)
        write(9,*)
        write(9,*)'The ',j,'-th dimension space:',lbound(j),' bohr ----',rbound(j),' bohr'
        write(9,*)
        write(9,*)
        write(9,*)
        write(9,*)
     END DO
     write(9,*)
     call flush(5)

  ENDIF


  !==============================================================================
  !==============================================================================
  !	Preparing Kinetic-energy matrice for later propagation use
  !==============================================================================
  !==============================================================================

  CALL MPI_BARRIER(comm2d, IERR)


  DO idim = 1, ndim
     do ireg = 1, num_reg(idim)

        ntmp = num_fun(idim,ireg)
        ALLOCATE( mat_reg(idim,ireg)%ke_mat(ntmp,ntmp), &
             &    mat_reg(idim,ireg)%eigvec_mat(ntmp,ntmp), &
             &    mat_reg(idim,ireg)%pt(ntmp),   &
             &    mat_reg(idim,ireg)%wt(ntmp),   &
             &    mat_reg(idim,ireg)%fac1(ntmp),   &
             &    mat_reg(idim,ireg)%fac2(ntmp),   &
             &    mat_reg(idim,ireg)%df(ntmp,ntmp),   &
             &    mat_reg(idim,ireg)%ddf(ntmp,ntmp),   &
             &    mat_reg(idim,ireg)%eigval_mat(ntmp))
     end do
  END DO


  CALL MPI_BARRIER(comm2d, IERR)

  !--------------------------------------------------------------------------------
  ! Set up finite element dvr to get points, weights, and Kinetic-energy matrice
  !--------------------------------------------------------------------------------



  write(6,*)'Myid=',MYID,'  Coming here before dvr_setup !'

  CALL dvr_setup

  write(6,*)'Myid=',MYID,'  Coming here after dvr_setup !'

  CALL MPI_BARRIER(comm2d, IERR)

  !--------------------------------------------------------------------------------
  ! Diagonalizing Kinetic-energy matrice for each dimension & region;
  ! Storing the propagators for later propagation use
  !--------------------------------------------------------------------------------


  !**********finding the maximum space to allocate the arrays**************

  maxregion=1
  maxfunc=1
  do j=1, ndim
     if(maxregion.lt.num_reg(j))maxregion=num_reg(j)
     do i=1, num_reg(j)
        if(maxfunc.lt.num_fun(j,i))maxfunc=num_fun(j,i)
     end do
  end do


  IF(MYID.eq.0)THEN
     write(9,*)
     write(9,*)'The maximum number of elements for each dimension:: maxregion=',maxregion
     write(9,*)'The maximum number of functions for each element::  maxfunc=',maxfunc
     write(9,*)
     call flush(5)
  ENDIF

  !*****************************************************************************

  IF( MYID.eq.0 )THEN

     ALLOCATE ( Kin_eigvec(ndim,maxregion,maxfunc,maxfunc), &
          Kin_eigvec_inv(ndim,maxregion,maxfunc,maxfunc), &
          Kin_eigval_dt(ndim,maxregion,maxfunc), &
          Kin_eigval_halfdt(ndim,maxregion,maxfunc), &
          Kin_eigval_quarterdt(ndim,maxregion,maxfunc), &
          Kin_operator_dt(ndim,maxregion,maxfunc,maxfunc), &
          Kin_operator_halfdt(ndim,maxregion,maxfunc,maxfunc), &
          Kin_operator_quarterdt(ndim,maxregion,maxfunc,maxfunc) )
  ELSE

     ALLOCATE ( Kin_operator_dt(ndim,maxregion,maxfunc,maxfunc), &
          Kin_operator_halfdt(ndim,maxregion,maxfunc,maxfunc), &
          Kin_operator_quarterdt(ndim,maxregion,maxfunc,maxfunc) )

  ENDIF

  !*****************************************************************************

  Kin_operator_dt(:,:,:,:)=zero
  Kin_operator_halfdt(:,:,:,:)=zero
  Kin_operator_quarterdt(:,:,:,:)=zero

  CALL MPI_BARRIER(comm2d, IERR)

  !---------------------------------------------------------------------------------
  !	Let only the Master CPU do the diagonalization!
  ! 	If let all CPU do it, will be wrong; but don't know why?
  !---------------------------------------------------------------------------------

  IF( MYID.eq.0 )THEN

     DO j=1, ndim

        DO i=1, num_reg(j)

           !--------------------considering throwing-out the first point for the first element-------------------

           if( i.eq.1 )then

              ntemp1=num_fun(j,1)-1

              ALLOCATE ( bmat_eigval(ntemp1), bmat_eigvec(ntemp1,ntemp1) )
              do m=2, num_fun(j,1)
                 do n=2, num_fun(j,1)
                    bmat_eigvec(m-1,n-1)=mat_reg(j,i)%ke_mat(m,n)
                 end do
              end do

              CALL LA_SYEV(bmat_eigvec,bmat_eigval,JOBZ='V')

              do m=1, num_fun(j,1)
                 do n=1, num_fun(j,1)
                    if(m.ne.1 .AND. n.ne.1)then	
                       mat_reg(j,i)%eigval_mat(m)=bmat_eigval(m-1)
                       mat_reg(j,i)%eigvec_mat(m,n)=bmat_eigvec(m-1,n-1)	
                    else
                       mat_reg(j,i)%eigval_mat(m)=0.d0
                       mat_reg(j,i)%eigvec_mat(m,n)=0.d0
                    endif
                 end do
              end do
              DEALLOCATE(bmat_eigvec,bmat_eigval)
           else
              mat_reg(j,i)%eigvec_mat(:,:) = mat_reg(j,i)%ke_mat(:,:)
              CALL LA_SYEV(mat_reg(j,i)%eigvec_mat,mat_reg(j,i)%eigval_mat,JOBZ='V')
           endif

           write(6,*)'MYID=',MYID,' j=',j,'  i=',i,'  Diagonalizing Kinetic-energy matrice....'

           !-----------------------------------------------------------------------	

           do k=1, num_fun(j,i)
              do m=1, num_fun(j,i)
                 Kin_eigvec(j,i,k,m)=mat_reg(j,i)%eigvec_mat(k,m)
                 Kin_eigvec_inv(j,i,k,m)=mat_reg(j,i)%eigvec_mat(m,k)
              end do
              Kin_eigval_dt(j,i,k)=exp(-dt*mat_reg(j,i)%eigval_mat(k))
              Kin_eigval_halfdt(j,i,k)=exp(-half*dt*mat_reg(j,i)%eigval_mat(k))
              Kin_eigval_quarterdt(j,i,k)=exp(-quarter*dt*mat_reg(j,i)%eigval_mat(k))
           end do


           !---------------------------------------------------------------------------------	
           !-------building the COMPACT Kinetic-energy operator------------------------------	
           !---------------------------------------------------------------------------------	


           do k=1, num_fun(j,i)
              do m=1, num_fun(j,i)
                 ctemp_mat1(k,m)=Kin_eigval_halfdt(j,i,k)*Kin_eigvec_inv(j,i,k,m)
                 ctemp_mat2(k,m)=Kin_eigval_quarterdt(j,i,k)*Kin_eigvec_inv(j,i,k,m)
                 ctemp_mat0(k,m)=Kin_eigval_dt(j,i,k)*Kin_eigvec_inv(j,i,k,m)
              end do
           end do


           do k=1, num_fun(j,i)
              do m=1, num_fun(j,i)

                 ctemp=zero
                 ctemp1=zero
                 ctemp2=zero
                 do n=1, num_fun(j,i)
                    ctemp=ctemp+Kin_eigvec(j,i,k,n)*ctemp_mat1(n,m)
                    ctemp1=ctemp1+Kin_eigvec(j,i,k,n)*ctemp_mat2(n,m)
                    ctemp2=ctemp2+Kin_eigvec(j,i,k,n)*ctemp_mat0(n,m)
                 end do
                 Kin_operator_dt(j,i,k,m)=ctemp2
                 Kin_operator_halfdt(j,i,k,m)=ctemp
                 Kin_operator_quarterdt(j,i,k,m)=ctemp1

              end do
           end do

           !-------------------------------------	


        END DO
     END DO

  ENDIF

  !---------------------------------------------------------------------------------
  !	Passing the KE-operator to each CPU
  !---------------------------------------------------------------------------------

  CALL MPI_BARRIER(comm2d, IERR)

  ntemp2=ndim*maxregion*maxfunc*maxfunc

  CALL MPI_BCAST(Kin_operator_dt(1,1,1,1),ntemp2,MPI_DOUBLE_PRECISION,0,comm2d, IERR)
  CALL MPI_BCAST(Kin_operator_halfdt(1,1,1,1),ntemp2,MPI_DOUBLE_PRECISION,0,comm2d, IERR)
  CALL MPI_BCAST(Kin_operator_quarterdt(1,1,1,1),ntemp2,MPI_DOUBLE_PRECISION,0,comm2d, IERR)

  CALL MPI_BARRIER(comm2d, IERR)


  !-------------------------------------	


  IF( MYID.eq.0 )THEN

     DEALLOCATE ( Kin_eigvec, Kin_eigvec_inv, Kin_eigval_halfdt, &
          Kin_eigval_quarterdt, Kin_eigval_dt )	

  ENDIF

  !--------------------------------------------------------------------------------
  !	Finishing the construction of Kinetic-energy operators
  !--------------------------------------------------------------------------------

  !--------------------------------------------------------------------------------
  ! Construct complete spatial grid for each dimension
  !--------------------------------------------------------------------------------

  maxpoint=1
  do j=1, ndim
     if(maxpoint.lt.ntot(j))maxpoint=ntot(j)
  end do

  ALLOCATE( grid(ndim,maxpoint) )

  write(6,*)
  write(6,*)'MYID=',MYID,'  Before setup_grid_mapping......'
  write(6,*)


  CALL setup_grid_mapping


  write(6,*)
  write(6,*)'MYID=',MYID,'  After setup_grid_mapping......'
  write(6,*)

  !-----------------Avoiding the singularity at the origin------------------

  do j=1, ndim	
     do i=1, ntot(j)
        if( dabs(grid(j,i)).le.1.d-15 )then
           grid(j,i)=grid(j,i)+1.d-8
        endif
     end do
  end do

  CALL MPI_BARRIER(comm2d, IERR)


  !-------------------------------------------------------------------------



  IF(MYID.eq.0)THEN
     open(55,file='space_points.dat')
     write(55,*)
     DO j=1, ndim
        write(55,*)'---------The ',j,'-th dimensional points are:'
        do i=1, ntot(j)
           write(55,801)i, grid(j,i)
        end do
        write(55,*)
     END DO
     close(55)
  ENDIF


  !=====================================================================
  !	Mapping the partial waves into one-dimensional array
  !=====================================================================


  ALLOCATE ( l1mat(500), l2mat(500), ltotmat(500) )


  l1mat(:)=0
  l2mat(:)=0
  ltotmat(:)=0

  !***************MASTER DOES THE FOLLOWING WORKS***************

  IF( MYID.EQ.0 )THEN


     mm=0
     DO LL=0,LLmax
        nparity=mod(LL,2)

        do l1=0,l1max
	   do l2=0,l2max

              IF( (LL.ge.abs(l1-l2)).AND.(L.le.abs(l1+l2)) )THEN

                 ntemp1=mod(l1+l2,2)

                 IF( nparity.eq.ntemp1 )THEN
                    mm=mm+1
                    ltotmat(mm)=LL
                    l1mat(mm)=l1
                    l2mat(mm)=l2
                 ENDIF

              ENDIF

	   end do
        end do

     END DO

     !--------------------------------
     nwave=mm
     !--------------------------------


     write(9,*)'=============================================='
     write(9,*)'       The partial waves used for expansion    '
     write(9,*)'=============================================='

     write(9,*)
     write(9,*)'The total number of L-l1-l2 combinations nwave:=',nwave

     write(9,*)
     write(9,*)'---------------------------------------------------------'
     write(9,*)'        I','         L','         l1','        l2'
     write(9,*)'---------------------------------------------------------'

     DO i=1,nwave

        write(9,*)
        write(9,810)i, ltotmat(i), l1mat(i), l2mat(i)

     END DO

     write(9,*)
     write(9,*)'================================================'

     call flush(5)

810  format(I10,I10,I10,I10)

     !***********************************************************************
     !	END of MASTER
     !***********************************************************************

  ENDIF

  !--------------------------------------------------

  CALL MPI_BARRIER(comm2d, IERR)

  CALL MPI_BCAST(nwave,1,MPI_INTEGER,0,comm2d, IERR)

  CALL MPI_BCAST(  l1mat(1:nwave),nwave,MPI_INTEGER,0,comm2d, IERR)
  CALL MPI_BCAST(  l2mat(1:nwave),nwave,MPI_INTEGER,0,comm2d, IERR)
  CALL MPI_BCAST(ltotmat(1:nwave),nwave,MPI_INTEGER,0,comm2d, IERR)

  CALL MPI_BARRIER(comm2d, IERR)

  write(6,*)
  write(6,*)'MYID=',MYID,'before diagonalizing potentials......'
  write(6,*)

  !--------------------------------------------------------------------------------
  !--------------------------------------------------------------------------------
  ! Construct the potential matrice
  !--------------------------------------------------------------------------------
  !--------------------------------------------------------------------------------

  ALLOCATE (    F0(0:nnp,nwave,nwave), F1(nwave,nwave), F2(nwave,nwave), &
       F111(nwave,nwave,ngridstart(1):ngridstop(1)), &
       F222(nwave,nwave,ngridstart(2):ngridstop(2)), &
       amat(nwave,nwave),BB(nwave,nwave),BT(nwave,nwave), &
       eiger(nwave), eigei(nwave),ivv1(nwave),fvv1(nwave), &
       VL(nwave,nwave,ngridstart(1):ngridstop(1),ngridstart(2):ngridstop(2)), &
       CGG(nwave,nwave,ngridstart(1):ngridstop(1),ngridstart(2):ngridstop(2)), &
       GG(nwave,nwave,ngridstart(1):ngridstop(1),ngridstart(2):ngridstop(2)), &
       GT(nwave,nwave,ngridstart(1):ngridstop(1),ngridstart(2):ngridstop(2)),&
       EG(nwave,ngridstart(1):ngridstop(1),ngridstart(2):ngridstop(2)) , &
       camat(nwave,nwave), &
       wavenorms(nwave), tmpwavenorms(nwave))

  !=========================================================================
  !=========================================================================

  !---------------------------------------------------------
  !    for F1 & F2
  !---------------------------------------------------------

  DO j=1,nwave

     l1=l1mat(j)
     l2=l2mat(j)
     LL=ltotmat(j)

     DO k=1,nwave

        l1p=l1mat(k)
        l2p=l2mat(k)
        Lp=ltotmat(k)

        IF(abs(LL-Lp).eq.1)THEN

	   temp1=dsqrt(9.d0*(2*l1p+1.d0)*(2*l2p+1.d0)*(2*Lp+1.d0)) &
	        /(two*twopi)

	   F1(j,k)=temp1*cleb(2,0,2*l1p,0,2*l1,0)* &
                cleb(0,0,2*l2p,0,2*l2,0)* &
                cleb(2,0,2*Lp,0,2*LL,0)* &
                ninej(2,2*l1p,2*l1,0,2*l2p,2*l2,2,2*Lp,2*LL)

	   F2(j,k)=temp1*cleb(0,0,2*l1p,0,2*l1,0)* &
                cleb(2,0,2*l2p,0,2*l2,0)* &
                cleb(2,0,2*Lp,0,2*LL,0)* &
                ninej(0,2*l1p,2*l1,2,2*l2p,2*l2,2,2*Lp,2*LL)
        ELSE
           F1(j,k)=0.d0
           F2(j,k)=0.d0
        ENDIF

     END DO
  END DO

  !---------------------------------------------------------
  !    set up angular momentum factors and
  !    radial parts for electron - electron interaction
  !---------------------------------------------------------

  !!! this only works if both dimensions actually have the same grid!!
  idim = 1
  NR = ntot(idim)
  NReff = NR-1
  ALLOCATE( V_r1r2lambda(NR,NR,0:nnp), Ttemp(2:NReff,2:NReff), Tinv(2:NReff,2:NReff), Tident(2:NReff,2:NReff), dvr_weights(NR) )
  
  dvr_weights(:) = 0.d0
  i1start = 1
  do ireg = 1, num_reg(idim)
     i1end = i1start + num_fun(idim,ireg) - 1
     dvr_weights(i1start:i1end) = dvr_weights(i1start:i1end) + mat_reg(idim,ireg)%wt(:)
     i1start = i1end
  end do
  if (i1end /= NR) then
     write(6,'(A,I5,A,I5)') 'i1end = ',i1end, ' /=  NR = ',NR
     stop 8
  end if

  DO lambda=0, nnp

     Ttemp(:,:) = 0.d0
     i1start = 2
     ireg = 1
     i1end = i1start + num_fun(idim,ireg) - 2 ! not -1 in first region because we don't take the starting point
     !----------- watch out! the factor here should be - 2.d0 because we actually want (d/dr^2 - lambda*(lambda+1)/r^2)
     !----------- but to have a positive definite (instead of a negative definite) matrix we take -T 
     !----------- (with T from Barry Schneider, Atomic.dvi)
     Ttemp(i1start:i1end,i1start:i1end) = Ttemp(i1start:i1end,i1start:i1end) + 2.d0 * mat_reg(idim,ireg)%ke_mat(2:,2:)
     i1start = i1end
     do ireg = 2, num_reg(idim)-1
        i1end = i1start + num_fun(idim,ireg) - 1
        !----------- the same thing here: +2 instead of -2
        Ttemp(i1start:i1end,i1start:i1end) = Ttemp(i1start:i1end,i1start:i1end) + 2.d0 * mat_reg(idim,ireg)%ke_mat(:,:)
        i1start = i1end
     end do
     ireg = num_reg(idim)
     i1end = i1start + num_fun(idim,ireg) - 2
     if (i1end /= NReff) stop 42
     !----------- the same thing here: +2 instead of -2
     ntmp = num_fun(idim,ireg) - 1
     Ttemp(i1start:i1end,i1start:i1end) = Ttemp(i1start:i1end,i1start:i1end) + 2.d0 * mat_reg(idim,ireg)%ke_mat(1:ntmp,1:ntmp)

     do ir1 = 2, NReff
        !----------- and here again - take +lambda*... instead of -lambda*...
        Ttemp(ir1,ir1) = Ttemp(ir1,ir1) + lambda*(lambda+1.d0)/grid(idim,ir1)**2
     end do
     
     ! invert Ttemp - use the expert driver with inverse iteration to get good results!
     ! matrix should be positive symmetric since we took -T
     
     ! set Tident to the identity matrix
     Tident(:,:) = 0.d0
     forall (ir1 = 2:NReff)
        Tident(ir1,ir1) = 1.d0
     end forall

     !SUBROUTINE LA_POSVX( A, B, X, UPLO=uplo, AF=af, FACT=fact, & 
     !     EQUED=equed,S=s, FERR=ferr, BERR=berr, & 
     !     RCOND=rcond, INFO=info ) 
     CALL LA_POSVX(Ttemp,Tident,Tinv,'U',FACT='E',RCOND=temp)
     ! rebuild lower part of matrix
     do ir1 = 2, NReff
        do ir2 = ir1+1, NReff
           Tinv(ir2,ir1) = Tinv(ir1,ir2)
        end do
     end do
     
     ! from barry schneider's calculation (express radial integral as solution to poisson's equation
     ! V_r1,r2^lambda = -(2*L+1) *  (T^(-1)_ir1,ir2) / (r_i * sqrt(w_i) * r_j * sqrt(w_j)) + (r_i*r_j)**lambda / r_N**(2*lambda+1)
     !----------- watch out - as the calculated Tinv is actually minus T^(-1) we have no minus in front of (2*lambda+1)

     !set first row and column to zero - not tested if this is necessary - in the time propagation they should be ignored anyway !?
     V_r1r2lambda(1,:,lambda) = 0.d0
     V_r1r2lambda(:,1,lambda) = 0.d0
     forall (ir1=2:NReff, ir2=2:NReff)
        V_r1r2lambda(ir1,ir2,lambda) = (2*lambda+1.d0) * Tinv(ir1,ir2) &
             &                       / (grid(idim,ir1) * grid(idim,ir2) * sqrt(dvr_weights(ir1)*dvr_weights(ir2)))
     end forall
     forall (ir1=2:NR, ir2=2:NR)
        V_r1r2lambda(ir1,ir2,lambda) = V_r1r2lambda(ir1,ir2,lambda) &
             &                       + (grid(idim,ir1) * grid(idim,ir2))**lambda / grid(idim,NR)**(2*lambda+1)
     end forall

     do tag=1,nwave

        l1=l1mat(tag)
        l2=l2mat(tag)
        LL=ltotmat(tag)

        do tagp=1,nwave

           l1p=l1mat(tagp)
           l2p=l2mat(tagp)
           Lp=ltotmat(tagp)

           if (LL == Lp) then
              ! taken from kenichi ishikawa's presentation
              F0(lambda,tag,tagp) = (-1)**(LL+l2+l2p) * &
                   &   sqrt( (2*l1+1.d0) * (2*l1p+1.d0) * (2*l2+1.d0) * (2*l2p+1.d0) ) * &
                   &   Wigner3J000(l1,lambda,l1p) * Wigner3J000(l2,lambda,l2p) * &
                   &   sixj(2*LL,2*l2p,2*l1p,2*lambda,2*l1,2*l2)
           else
              F0(lambda,tag,tagp) = 0.d0
           end if
        end do
     end do
  END DO
  
  DEALLOCATE(Ttemp, Tinv, Tident)

  !=============================================================================

  DO i=1,nwave

     DO j=1,nwave

        !========================
        IF(F1(i,j).ne.0.d0)THEN
           do ir1=ngridstart(1), ngridstop(1)
              F111(i,j,ir1)= grid(1,ir1)*F1(i,j)*7.255197457d0
           end do
        ELSE
           do ir1=ngridstart(1), ngridstop(1)
              F111(i,j,ir1)= 0.d0
           end do
        END IF

        !========================

        IF(F2(i,j).ne.0.d0)THEN
           do ir2=ngridstart(2), ngridstop(2)
              F222(i,j,ir2)= grid(2,ir2)*F2(i,j)*7.255197457d0
           end do
        ELSE
           do ir2=ngridstart(2), ngridstop(2)
              F222(i,j,ir2)= 0.d0
           end do
        END IF

        !========================

     END DO
  END DO

  !=========================================================================
  !=========================================================================
  !	Diagonalizing the potentials VL, F111 and F222
  !=========================================================================
  !=========================================================================

  DO ir1=ngridstart(1), ngridstop(1)

     write(6,*)'MYID=',MYID,'  ir1=',ir1,'  Diagonalizing VL...'

     DO ir2=ngridstart(2), ngridstop(2)


        do tag=1,nwave
           do tagp=1,nwave

              !------------------------------------------------------------------------

              LL=ltotmat(tag)
              Lp=ltotmat(tagp)

              IF(LL.eq.Lp)THEN

                 BB(tag,tagp)= sum(V_r1r2lambda(ir1,ir2,0:nnp) * F0(0:nnp,tag,tagp))

              ELSE

	         BB(tag,tagp)=0.d0

              ENDIF

              VL(tag,tagp,ir1,ir2) = BB(tag,tagp)

              !------------------------------------------------------------------------

           end do
        end do

        CALL LA_SYEV(BB,eiger,JOBZ='V')
        BT(:,:) = transpose(BB(:,:))

        do i=1,nwave
           EG(i,ir1,ir2)=exp(-dt*eiger(i))
           do j=1,nwave
              GG(i,j,ir1,ir2)=BB(i,j)
              GT(i,j,ir1,ir2)=BT(i,j)
	   end do
        end do

        if (MYID == 0) then
           BB(:,:) = matmul(BT(:,:),BB(:,:))
           do i=1,nwave
              BB(i,i) = BB(i,i) - 1.d0
           end do
           temp1 = sum(abs(BB(:,:)))
           if (temp1 > 1.d-10) then
              write(777,*) "BB not orthogonal, temp1 =", temp1
              do j=1,nwave
                 write(777,803) BB(:,j)
              end do
              write(777,*) ''; write(777,*) ''
           end if
        end if

        !---------------!---------------!---------------!---------------!---------------

     END DO
  END DO

  !---------------------------------------------------------------------------
  !	Combining the matrix GG*EG*GT ===> GG(nwave,nwave,ir1,ir2)
  !---------------------------------------------------------------------------

  DO ir1=ngridstart(1), ngridstop(1)	
     DO ir2=ngridstart(2), ngridstop(2)

        do i=1, nwave
           do j=1, nwave
              ctemp=zero
              do k=1, nwave
                 ctemp=ctemp+GG(i,k,ir1,ir2)*EG(k,ir1,ir2)*GT(k,j,ir1,ir2)
              end do
              camat(i,j)=ctemp	
           end do
        end do

        do i=1, nwave
           do j=1, nwave
              CGG(i,j,ir1,ir2)=camat(i,j)	
           end do
        end do

     END DO
  END DO


  DEALLOCATE ( EG, GT, GG, camat )


  !=========================================================================
  !=========================================================================

  ALLOCATE ( HH(nwave,nwave,ngridstart(1):ngridstop(1),ngridstart(2):ngridstop(2)), &
       HT(nwave,nwave,ngridstart(1):ngridstop(1),ngridstart(2):ngridstop(2)),&
       EH(nwave,ngridstart(1):ngridstop(1),ngridstart(2):ngridstop(2))     )


  DO ir1=ngridstart(1), ngridstop(1)

     write(6,*)'MYID=',MYID,'  ir1=',ir1,'  Diagonalizing F111+F222...'

     DO ir2=ngridstart(2), ngridstop(2)

        do i=1,nwave
           do j=1,nwave
              BB(i,j)=F111(i,j,ir1)+F222(i,j,ir2)
           end do
        end do

        !if (MYID==0) then
        !   temp1 = sum(abs(BB(:,:)-transpose(BB(:,:))))
        !   if (temp1 > 1.d-10) then
        !      write(777,*) "F111+F222 BB not symmetric, temp1 =", temp1
        !      write(777,*) "i =",i, " j =", j
        !      do j=1,nwave
        !         write(777,803) BB(:,j)
        !      end do
        !      write(777,*) ''; write(777,*) ''
        !   end if
        !end if

        CALL LA_SYEV(BB,eiger,JOBZ='V')
        BT(:,:) = transpose(BB(:,:))

        do i=1,nwave
           EH(i,ir1,ir2)=eiger(i)
           do j=1,nwave
              HH(i,j,ir1,ir2)=BB(i,j)
              HT(i,j,ir1,ir2)=BT(i,j)
	   end do
        end do

     END DO
  END DO

  !=========================================================================
  !=========================================================================

  DEALLOCATE ( F111, F222 )


  write(6,*)
  write(6,*)'MYID=',MYID,'After diagonalizing potentials......'
  write(6,*)




  !==============================================================================
  !==============================================================================
  !	Allocating spatial grids & coefficients (i.e., wavefunction) arrays
  !==============================================================================
  !==============================================================================


  write(6,*)
  write(6,*)'MYID=',MYID,'Before allocate psi......'
  write(6,*)

  ALLOCATE( psi3d( ngridstart(1):ngridstop(1), ngridstart(2):ngridstop(2), nwave ) )

  ALLOCATE( psi03d( ngridstart(1):ngridstop(1), ngridstart(2):ngridstop(2), nwave ) )

  ALLOCATE( coeff3d( ngridstart(1):ngridstop(1), ngridstart(2):ngridstop(2), nwave ))

  ALLOCATE( cpot3d( ngridstart(1):ngridstop(1), ngridstart(2):ngridstop(2), nwave ) )

  ALLOCATE( rpot3d( ngridstart(1):ngridstop(1), ngridstart(2):ngridstop(2), nwave ) )


  write(6,*)
  write(6,*)'MYID=',MYID,'After allocate psi and pot......'
  write(6,*)

  !=====================================================================
  !--------------Setting the abosorption edges--------------------------
  !=====================================================================

  ALLOCATE ( vv1i(ntot(1)), vv2i(ntot(2)) )

  do ix=1,ntot(1)

     if(ix.ge.(ntot(1)-nabx))then
        vv1i(ix)=1.d0 + VRi*dexp(-alphai*(grid(1,ntot(1))-grid(1,ix)))
     else
        vv1i(ix)=1.d0
     endif

  end do

  do iy=1,ntot(2)
     if(iy.ge.(ntot(2)-naby))then
        vv2i(iy)=1.d0 + VRi*dexp(-alphai*(grid(2,ntot(2))-grid(2,iy)))
     else
        vv2i(iy)=1.d0
     endif
  end do

  !---------------------------------------------------------
  IF(MYID.EQ.0)THEN

     open(887,file='absorr1.dat')
     open(888,file='absorr2.dat')

     do ix=1,ntot(1)
        write(887,*)grid(1,ix), vv1i(ix)
     end do
     do ix=1,ntot(2)
        write(888,*)grid(2,ix), vv2i(ix)
     end do

     close(887)
     close(888)

  ENDIF
  !---------------------------------------------------------


  !=========================================================================
  !	constructing the diagonal potential
  !=========================================================================

  write(6,*)
  write(6,*)'MYID=',MYID,'Before constructing diagonal potential cpot3d......'
  write(6,*)


  DO i=1,nwave
     l1=l1mat(i)
     l2=l2mat(i)
     DO ir1=ngridstart(1), ngridstop(1)
        DO ir2=ngridstart(2), ngridstop(2)

           rpot3d(ir1,ir2,i)=l1*(l1+one)/two/(grid(1,ir1)**2) + l2*(l2+one)/two/(grid(2,ir2)**2) &
                -Z1/grid(1,ir1) -Z1/grid(2,ir2)

           cpot3d(ir1,ir2,i)=exp(-dt*rpot3d(ir1,ir2,i))
           !*exp(-(vv1i(ir1)*vv2i(ir2)-1.d0))

	END DO
     END DO
  END DO


  write(6,*)
  write(6,*)'MYID=',MYID,'After constructing diagonal potential cpot3d......'
  write(6,*)



  !========================================================================================================
  !	Preparing temporary array for message passing
  !========================================================================================================

  npointx=(ngridstop(1)-ngridstart(1)) +1
  npointy=(ngridstop(2)-ngridstart(2)) +1
  npointz=nwave
  npointxy=npointx*npointy
  npointyz=npointz*npointy
  npointxz=npointz*npointx
  ALLOCATE ( prob2dtemp(npointx,npointy), probdens2d(ntot(1),ntot(2)) )
  ALLOCATE ( psitemp2dxy(npointx,npointy), psitemp2d(ntot(1),ntot(2)) )
  ALLOCATE ( psitemp2dyz(npointy,npointz), psitemp2dxz(npointx,npointz) )


  !========================================================================================================
  !========================================================================================================
  !	Loading initial wave packet & using imaginary dumping it to the GS (He atom)
  !========================================================================================================
  !========================================================================================================


  write(6,*)'MYID=',MYID,'  Before loading the initial wave packet....'

  ! TODO: check if GS.st exists and load it
  ! include checks for correct parameters?

  CALL MPI_BARRIER(comm2d, IERR)

  !write(0,*) 'jooodbg 0001000'
  
  do k=ngridstart(1), ngridstop(1)
     do m=ngridstart(2), ngridstop(2)

        !write(0,*) 'jooodbg 0001100, k=',k,'  m=',m

        ! use He^+ 1s^2 as guess for ground state
        coeff3d(k,m,1) = 2*Z1**1.5d0 * grid(1,k) * exp(-Z1*grid(1,k))/factor(1,k) * &
             &           2*Z1**1.5d0 * grid(2,m) * exp(-Z1*grid(2,m))/factor(2,m)

        if(k.eq.1 .OR. m.eq.1)then
           coeff3d(k,m,1) = zero
           coeff3d(k,m,2) = zero
        endif

        do n=2, nwave
           coeff3d(k,m,n)= zero
        end do
     end do
     !write(0,*) 'jooodbg 0002000'
  end do

  !write(0,*) 'jooodbg 0003000'

  CALL MPI_BARRIER(comm2d, IERR)

  !write(0,*) 'jooodbg 0004000'

  write(6,*)'MYID=',MYID,'  After loading the initial wave packet....'
  call flush(6)

  !========================================================================================================
  !-----------------------Normalization test-------------------------------------
  !========================================================================================================

  !------------------------------building back the wavefunction-------------------------	


  do k=ngridstart(1), ngridstop(1)
     do m=ngridstart(2), ngridstop(2)
        do n=1, nwave
           psi3d(k,m,n)=coeff3d(k,m,n)*factor(1,k)*factor(2,m)
        end do
     end do
  end do

  CALL MPI_BARRIER(comm2d, IERR)


  !------------------------------getting the wavefunction normalization--------------	


  write(6,*)'MYID=',MYID,'  Before initial normalization......'
  call flush(6)

  temp=0.d0
  tmpwavenorms(:) = 0.d0
  DO n=1, nwave
     DO j=regstart(1),regstop(1)
        DO i=regstart(2), regstop(2)

           do k=1, num_fun(1,j)
              do m=1, num_fun(2,i)
                 prob2d(k,m)=0.d0
              end do
           end do

           if(j.eq.1)then
	      startj=2
           else
	      startj=1
           endif

           if(i.eq.1)then
	      starti=2
           else
	      starti=1
           endif

           do k=startj, num_fun(1,j)
              do m=starti, num_fun(2,i)
                 prob2d(k,m) = psi3d(nindex(1,j,k),nindex(2,i,m),n)**2
              end do
           end do

           CALL integral2d(num_fun(1,j),num_fun(2,i),mat_reg(1,j)%wt, &
                mat_reg(2,i)%wt,prob2d,pnorm)

           tmpwavenorms(n) = tmpwavenorms(n) + pnorm
        END DO
     END DO
  END DO
  temp=sum(tmpwavenorms(:))

  write(6,*)'MYID=',MYID, ' Norm=',temp

  CALL MPI_REDUCE(tmpwavenorms, wavenorms, nwave, MPI_DOUBLE_PRECISION, MPI_SUM, 0, comm2d, IERR)
  CALL MPI_REDUCE(temp, temp2, 1, MPI_DOUBLE_PRECISION, MPI_SUM, 0, comm2d, IERR)

  IF(MYID.eq.0)THEN
     temp3=1.d0/sqrt(temp2)
  ENDIF

  CALL MPI_BARRIER(comm2d, IERR)
  CALL MPI_BCAST(temp3,1,MPI_DOUBLE_PRECISION,0,comm2d, IERR)
  CALL MPI_BARRIER(comm2d, IERR)

  if(MYID.eq.0)then
     write(9,*)'The initial norm of the trial wavefunction is: ',temp2
     write(9,*)
     call flush(5)
  endif

  !-----------------normalizing----------------------------
  do k=ngridstart(1), ngridstop(1)
     do m=ngridstart(2), ngridstop(2)
        do n=1, nwave
           coeff3d(k,m,n)=coeff3d(k,m,n)*temp3
        end do
     end do
  end do

  write(6,*)'MYID=',MYID,'  After initial normalization......'

  !==============================================================================================
  !	Storing the initial state
  !==============================================================================================

  do k=ngridstart(1), ngridstop(1)
     do m=ngridstart(2), ngridstop(2)
        do n=1, nwave
           psi03d(k,m,n)=coeff3d(k,m,n)*factor(1,k)*factor(2,m)
           psi3d(k,m,n)=psi03d(k,m,n)
        end do
     end do
  end do

  write(6,*)'MYID=',MYID,'  After storing the initial wave functions......'

  CALL MPI_BARRIER(comm2d, IERR)

  !====================================================================================
  !   Recording the initial wave packet
  !====================================================================================

  m=0
  do k=ngridstart(1), ngridstop(1)
     m=m+1
     l=0
     do n=ngridstart(2), ngridstop(2)
        l=l+1

        temp=0.d0
        do i=1, nwave
           temp=temp + psi3d(k,n,i)**2
        end do

        prob2dtemp(m,l)=temp
        if( prob2dtemp(m,l).le.1.d-35)prob2dtemp(m,l)=1.d-35

     end do
  end do


  if(MYID.ne.0)then

     CALL MPI_SSEND(ngridstart(1),1,MPI_INTEGER,0,MYID+2*ncount, comm2d, IERR)
     CALL MPI_SSEND(ngridstart(2),1,MPI_INTEGER,0,MYID+3*ncount, comm2d, IERR)
     CALL MPI_SSEND(prob2dtemp,npointxy,MPI_DOUBLE_PRECISION,0,MYID+4*ncount, comm2d, IERR)

  else

     do k=ngridstart(1), ngridstop(1)
        do n=ngridstart(2), ngridstop(2)
           probdens2d(k,n)=prob2dtemp(k,n)
        end do
     end do

     do j=1, nrcount

        CALL MPI_RECV(ntemp1,1,MPI_INTEGER,j,j+2*ncount, comm2d, status, IERR)
        CALL MPI_RECV(ntemp2,1,MPI_INTEGER,j,j+3*ncount, comm2d, status, IERR)
        CALL MPI_RECV(prob2dtemp,npointxy,MPI_DOUBLE_PRECISION,j,j+4*ncount, comm2d, status, IERR)

        do k=1, npointx
           do n=1, npointy
	      probdens2d(ntemp1+k-1,ntemp2+n-1)=prob2dtemp(k,n)
           end do
        end do

     end do


     !************Recording**********************************

     open(80,file='prob_0.dat')
     write(80,'(A,G10.5,A,I5,A,I5,A)')'#ZONE  T="',0.0,'", I=', ntot(2), ', J=', ntot(1),',  F=POINT'
     do i=1, ntot(1)
        do j=1, ntot(2)
           write(80,803)grid(1,i),grid(2,j),probdens2d(i,j)
        end do
     end do
     close(80)

  endif

  CALL MPI_BARRIER(comm2d, IERR)

  write(6,*)'MYID=',MYID,'  After recording the initial wave packets......'

  CALL MPI_BARRIER(comm2d, IERR)

  !==============================================================================
  !--------------At the beginning: Calculating the norm, <r1> <r2> -------------
  !==============================================================================


  time=0.d0


  temp=0.d0
  temp1=0.d0
  temp2=0.d0
  tmpwavenorms(:) = 0.d0
  DO n=1, nwave
     DO j=regstart(1),regstop(1)
        DO i=regstart(2), regstop(2)


           do k=1, num_fun(1,j)
              do m=1, num_fun(2,i)
                 prob2d(k,m)=0.d0
                 probx2d(k,m)=0.d0
                 proby2d(k,m)=0.d0
              end do
           end do

           if(j.eq.1)then
	      startj=2
           else
	      startj=1
           endif

           if(i.eq.1)then
	      starti=2
           else
	      starti=1
           endif

           do k=startj, num_fun(1,j)
              do m=starti, num_fun(2,i)
                 prob2d(k,m)=psi3d(nindex(1,j,k),nindex(2,i,m),n)**2
                 probx2d(k,m)=prob2d(k,m)*grid(1,nindex(1,j,k))
                 proby2d(k,m)=prob2d(k,m)*grid(2,nindex(2,i,m))
              end do
           end do


           CALL integral2d(num_fun(1,j),num_fun(2,i),mat_reg(1,j)%wt, &
                mat_reg(2,i)%wt,prob2d,pnorm)
           CALL integral2d(num_fun(1,j),num_fun(2,i),mat_reg(1,j)%wt, &
                mat_reg(2,i)%wt,probx2d,xexp)
           CALL integral2d(num_fun(1,j),num_fun(2,i),mat_reg(1,j)%wt, &
                mat_reg(2,i)%wt,proby2d,yexp)

           tmpwavenorms(n) = tmpwavenorms(n) + pnorm
           temp=temp+pnorm
           temp1=temp1+xexp
           temp2=temp2+yexp

        END DO
     END DO
  END DO

  CALL MPI_REDUCE(tmpwavenorms, wavenorms, nwave, MPI_DOUBLE_PRECISION, MPI_SUM, 0, comm2d, IERR)
  CALL MPI_REDUCE(temp,  pnorm, 1, MPI_DOUBLE_PRECISION, MPI_SUM, 0, comm2d, IERR)
  CALL MPI_REDUCE(temp1, xexp,  1, MPI_DOUBLE_PRECISION, MPI_SUM, 0, comm2d, IERR)
  CALL MPI_REDUCE(temp2, yexp,  1, MPI_DOUBLE_PRECISION, MPI_SUM, 0, comm2d, IERR)

  IF(MYID.eq.0)THEN
     write(10,802) time*convfs, pnorm, wavenorms
     write(11,802) time*convfs, xexp
     write(12,802) time*convfs, yexp
     write(14,803) xexp, yexp
     call flush(10)	
     call flush(11)	
     call flush(12)	
     call flush(14)	
  ENDIF

  CALL MPI_BARRIER(comm2d, IERR)

  !==============================================================================
  !--------------At the beginning: Calculating the <pr1> <pr2>-------------
  !==============================================================================

  !--------Radial momentum : <pr>=<-i hbar psi/r - i hbar d(psi)/dr >---------------------------------------

  time=0.d0

  ctemp4=zero
  ctemp5=zero

  DO l=1, nwave
     DO j=regstart(1),regstop(1)
        DO i=regstart(2), regstop(2)

           do k=1, num_fun(1,j)
              do m=1, num_fun(2,i)
                 cprobpx2d(k,m)=zero
                 cprobpy2d(k,m)=zero
              end do
           end do

           if(j.eq.1)then
	      startj=2
           else
	      startj=1
           endif

           if(i.eq.1)then
	      starti=2
           else
	      starti=1
           endif

           do k=startj, num_fun(1,j)
              do m=starti, num_fun(2,i)

                 ctemp=zero
                 do n=startj, num_fun(1,j)
                    ctemp=ctemp + mat_reg(1,j)%df(k,n)*psi3d(nindex(1,j,n),nindex(2,i,m),l)
                 end do
                 ctemp_mat4(k,m,1)=ctemp + psi3d(nindex(1,j,k),nindex(2,i,m),l)/grid(1,nindex(1,j,k))


                 ctemp=zero
                 do n=starti, num_fun(2,i)
                    ctemp=ctemp+ mat_reg(2,i)%df(m,n)*psi3d(nindex(1,j,k),nindex(2,i,n),l)
                 end do
                 ctemp_mat5(k,m,1)=ctemp + psi3d(nindex(1,j,k),nindex(2,i,m),l)/grid(2,nindex(2,i,m))

              end do
           end do

           do k=startj, num_fun(1,j)
              do m=starti, num_fun(2,i)	
                 cprobpx2d(k,m)= psi3d(nindex(1,j,k),nindex(2,i,m),l)*ctemp_mat4(k,m,1)
                 cprobpy2d(k,m)= psi3d(nindex(1,j,k),nindex(2,i,m),l)*ctemp_mat5(k,m,1)
              end do
           end do

           CALL integral2d(num_fun(1,j),num_fun(2,i),mat_reg(1,j)%wt,mat_reg(2,i)%wt,cprobpx2d,ctemp1)
           CALL integral2d(num_fun(1,j),num_fun(2,i),mat_reg(1,j)%wt,mat_reg(2,i)%wt,cprobpy2d,ctemp2)

           ctemp4=ctemp4+ctemp1
           ctemp5=ctemp5+ctemp2

        END DO
     END DO
  END DO

  IF(MYID.ne.0)THEN
     CALL MPI_SSEND(ctemp4,1,MPI_DOUBLE_PRECISION,0,MYID, comm2d, IERR)
     CALL MPI_SSEND(ctemp5,1,MPI_DOUBLE_PRECISION,0,MYID+ncount, comm2d, IERR)
  ELSE
     do j=1, nrcount
        CALL MPI_RECV(ctemp1,1,MPI_DOUBLE_PRECISION,j,j, comm2d, status, IERR)
        CALL MPI_RECV(ctemp2,1,MPI_DOUBLE_PRECISION,j,j+ncount, comm2d, status, IERR)
        ctemp4=ctemp4+ctemp1
        ctemp5=ctemp5+ctemp2
     end do

     write(111,802) time*convfs, ( ctemp4 )
     write(112,802) time*convfs, ( ctemp5 )
     write(114,803) ( ctemp4 ), ( ctemp5 )
     call flush(111)	
     call flush(112)	
     call flush(114)	

  ENDIF


  !==============================================================================

  CALL MPI_BARRIER(comm2d, IERR)

  !==============================================================================
  !==============================================================================
  !==============================================================================
  !==============================================================================
  !==============================================================================
  !==============================================================================
  !==============================================================================
  !	Starting propagation of the wave function (i.e., coefficients)
  !==============================================================================
  !==============================================================================
  !==============================================================================
  !==============================================================================
  !==============================================================================
  !==============================================================================
  !==============================================================================


  write(6,*)'  COME HERE STARTING TIME PROPAGATION.....'

  time=tstart
  npro=0
  starttime=MPI_Wtime()


  DO it=1, ntim

     ttx1=MPI_Wtime()

     time=time+dt


     IF(MYID.eq.0)THEN
        write(6,*)
        write(6,*)'it=',it,'  Time=',time*convfs
     ENDIF

     !------------getting the field amplitudes at time----------------------


     !	  CALL timfield( time, efldx)

     !	  if(MYID.eq.0 .AND. mod(it,1).eq.0)then
     !	   write(16,802)time*convfs, efldx
     !	   call flush(16)
     !	  endif


     !==============================================================================
     !==============================================================================
     !==============================================================================
     !---------------Kinetic-energy operation in multi-dimensions-----------
     !==============================================================================
     !==============================================================================
     !==============================================================================

     !==============================================================================
     !==============================================================================

     !------------------------------- 3D x-odd ------------------------------------------------------	

     CALL KE_quarterdt_3D_xodd

     CALL MPI_BARRIER(comm2d, IERR)

     ! 	write(6,*)'MYID=',MYID,' After first x-odd.....'

     !------------------------------------passing the boundary points---------------------------------------------------	

     m=0
     do i=ngridstart(2), ngridstop(2)
        m=m+1
        do j=1, nwave
           psitemp2dyz(m,j)=coeff3d(ngridstart(1),i,j)
        end do
     end do

     CALL MPI_SSEND(psitemp2dyz,npointyz,MPI_DOUBLE_PRECISION,nbrup,MYID, comm2d, IERR)

     CALL MPI_RECV(psitemp2dyz,npointyz,MPI_DOUBLE_PRECISION,nbrdown,nbrdown, comm2d,status,IERR)

     l=0
     do i=ngridstart(2), ngridstop(2)
        l=l+1
        do j=1, nwave
           coeff3d(ngridstop(1),i,j)=psitemp2dyz(l,j)
        end do
     end do


     CALL MPI_BARRIER(comm2d, IERR)

     ! 	write(6,*)'MYID=',MYID,' after first x-odd; passing/receiving.....'


     !------------------------------- 3D x-even ------------------------------------------------------	


     CALL KE_halfdt_3D_xeven

     CALL MPI_BARRIER(comm2d, IERR)

     ! 	write(6,*)'MYID=',MYID,' After x-even.....'

     !------------------------------------passing the boundary points---------------------------------------------------	

     m=0
     do i=ngridstart(2), ngridstop(2)
        m=m+1
        do j=1, nwave
           psitemp2dyz(m,j)=coeff3d(ngridstop(1),i,j)
        end do
     end do

     CALL MPI_SSEND(psitemp2dyz,npointyz,MPI_DOUBLE_PRECISION,nbrdown,MYID, comm2d, IERR)

     CALL MPI_RECV(psitemp2dyz,npointyz,MPI_DOUBLE_PRECISION,nbrup,nbrup, comm2d,status,IERR)

     l=0
     do i=ngridstart(2), ngridstop(2)
        l=l+1
        do j=1, nwave
           coeff3d(ngridstart(1),i,j)=psitemp2dyz(l,j)
        end do
     end do

     CALL MPI_BARRIER(comm2d, IERR)

     ! 	write(6,*)'MYID=',MYID,' After x-even; passing/receiving.....'

     !------------------------------- 3D x-odd AGAIN ------------------------------------------------	

     CALL KE_quarterdt_3D_xodd

     CALL MPI_BARRIER(comm2d, IERR)

     ! 	write(6,*)'MYID=',MYID,' After second x-odd.....'

     !------------------------------------passing the boundary points---------------------------------------------------	


     m=0
     do i=ngridstart(2), ngridstop(2)
        m=m+1
        do j=1, nwave
           psitemp2dyz(m,j)=coeff3d(ngridstart(1),i,j)
        end do
     end do

     CALL MPI_SSEND(psitemp2dyz,npointyz,MPI_DOUBLE_PRECISION,nbrup,MYID, comm2d, IERR)

     CALL MPI_RECV(psitemp2dyz,npointyz,MPI_DOUBLE_PRECISION,nbrdown,nbrdown, comm2d,status,IERR)

     l=0
     do i=ngridstart(2), ngridstop(2)
        l=l+1
        do j=1, nwave
           coeff3d(ngridstop(1),i,j)=psitemp2dyz(l,j)
        end do
     end do


     CALL MPI_BARRIER(comm2d, IERR)

     ! 	write(6,*)'MYID=',MYID,' After second x-odd; passing/receiving.....'


     !==============================================================================
     !	Carrying out y-dimension operations
     !==============================================================================

     !------------------------------- 3D y-odd ------------------------------------------------------	

     CALL KE_quarterdt_3D_yodd

     CALL MPI_BARRIER(comm2d, IERR)

     ! 	write(6,*)'MYID=',MYID,' after first y-odd.....'

     !------------------------------------passing the boundary points---------------------------------------------------	


     m=0
     do i=ngridstart(1), ngridstop(1)
        m=m+1
        do j=1, nwave
           psitemp2dxz(m,j)=coeff3d(i,ngridstart(2),j)
        end do
     end do

     CALL MPI_SSEND(psitemp2dxz,npointxz,MPI_DOUBLE_PRECISION,nbrleft,MYID, comm2d, IERR)

     CALL MPI_RECV(psitemp2dxz,npointxz,MPI_DOUBLE_PRECISION,nbrright,nbrright, comm2d,status,IERR)

     l=0
     do i=ngridstart(1), ngridstop(1)
        l=l+1
        do j=1, nwave
           coeff3d(i,ngridstop(2),j)=psitemp2dxz(l,j)
        end do
     end do


     CALL MPI_BARRIER(comm2d, IERR)

     ! 	write(6,*)'MYID=',MYID,' after first y-odd; passing/receiving.....'


     !------------------------------- 3D y-even ------------------------------------------------------	


     CALL KE_halfdt_3D_yeven

     CALL MPI_BARRIER(comm2d, IERR)

     ! 	write(6,*)'MYID=',MYID,' after y-even.....'

     !------------------------------------passing the boundary points---------------------------------------------------	


     m=0
     do i=ngridstart(1), ngridstop(1)
        m=m+1
        do j=1, nwave
           psitemp2dxz(m,j)=coeff3d(i,ngridstop(2),j)
        end do
     end do

     CALL MPI_SSEND(psitemp2dxz,npointxz,MPI_DOUBLE_PRECISION,nbrright,MYID, comm2d, IERR)

     CALL MPI_RECV(psitemp2dxz,npointxz,MPI_DOUBLE_PRECISION,nbrleft,nbrleft, comm2d,status,IERR)

     l=0
     do i=ngridstart(1), ngridstop(1)
        l=l+1
        do j=1, nwave
           coeff3d(i,ngridstart(2),j)=psitemp2dxz(l,j)
        end do
     end do


     CALL MPI_BARRIER(comm2d, IERR)

     !	write(6,*)'MYID=',MYID,' after y-even; passing/receiving.....'

     !------------------------------- 3D y-odd AGAIN ------------------------------------------------	

     CALL KE_quarterdt_3D_yodd

     CALL MPI_BARRIER(comm2d, IERR)

     ! 	write(6,*)'MYID=',MYID,' after second y-odd.....'

     !------------------------------------passing the boundary points---------------------------------------------------	


     m=0
     do i=ngridstart(1), ngridstop(1)
        m=m+1
        do j=1, nwave
           psitemp2dxz(m,j)=coeff3d(i,ngridstart(2),j)
        end do
     end do

     CALL MPI_SSEND(psitemp2dxz,npointxz,MPI_DOUBLE_PRECISION,nbrleft,MYID, comm2d, IERR)

     CALL MPI_RECV(psitemp2dxz,npointxz,MPI_DOUBLE_PRECISION,nbrright,nbrright, comm2d,status,IERR)

     l=0
     do i=ngridstart(1), ngridstop(1)
        l=l+1
        do j=1,  nwave
           coeff3d(i,ngridstop(2),j)=psitemp2dxz(l,j)
        end do
     end do


     CALL MPI_BARRIER(comm2d, IERR)

     ! 	write(6,*)'MYID=',MYID,' after second y-odd; passing/receiving.....'





     !==============================================================================
     !==============================================================================
     !==============================================================================
     !----------------------Propagating through the linear potential-------------------
     !==============================================================================
     !==============================================================================
     !==============================================================================


     !-----------------the diagonal part-----------------------------------

     do k=ngridstart(1), ngridstop(1)
        do m=ngridstart(2), ngridstop(2)
           do n=1, nwave

              IF(k.ne.1 .AND. m.ne.1)THEN
                 coeff3d(k,m,n)=coeff3d(k,m,n)*cpot3d(k,m,n)
              ENDIF

           end do
        end do
     end do

     !-----------------the 1/r12 interaction part-----------------------------------


     do k=ngridstart(1), ngridstop(1)
        do m=ngridstart(2), ngridstop(2)

           IF(k.ne.1 .AND. m.ne.1)THEN

              do n=1, nwave
                 ctemp=zero
                 do i=1, nwave
                    ctemp=ctemp+CGG(n,i,k,m)*coeff3d(k,m,i)
                 end do
                 ctemp_mat1(n,1)=ctemp
              end do

              do n=1, nwave
                 coeff3d(k,m,n)=ctemp_mat1(n,1)
              end do

           ENDIF


        end do
     end do

     !-----------------the field interaction part-----------------------------------

     GOTO 5678

     do k=ngridstart(1), ngridstop(1)
        do m=ngridstart(2), ngridstop(2)


	   do n=1, nwave
              ctemp=zero
              do i=1, nwave
                 ctemp=ctemp+HT(n,i,k,m)*coeff3d(k,m,i)
              end do
              ctemp_mat1(n,1)=ctemp*exp(-dt*EH(n,k,m)*efldx)
	   end do

	   do n=1, nwave
              ctemp=zero
              do i=1, nwave
                 ctemp=ctemp+HH(n,i,k,m)*ctemp_mat1(i,1)
              end do
              coeff3d(k,m,n)=ctemp
	   end do


        end do
     end do

5678 continue

     CALL MPI_BARRIER(comm2d, IERR)


     !==============================================================================
     !==============================================================================
     !==============================================================================
     !---------------AGAIN Kinetic-energy operation in multi-dimensions-----------
     !==============================================================================
     !==============================================================================
     !==============================================================================

     !------------------------------- 3D x-odd ------------------------------------------------------	

     CALL KE_quarterdt_3D_xodd

     CALL MPI_BARRIER(comm2d, IERR)

     ! 	write(6,*)'MYID=',MYID,' After first x-odd.....'

     !------------------------------------passing the boundary points---------------------------------------------------	

     m=0
     do i=ngridstart(2), ngridstop(2)
        m=m+1
        do j=1, nwave
           psitemp2dyz(m,j)=coeff3d(ngridstart(1),i,j)
        end do
     end do

     CALL MPI_SSEND(psitemp2dyz,npointyz,MPI_DOUBLE_PRECISION,nbrup,MYID, comm2d, IERR)

     CALL MPI_RECV(psitemp2dyz,npointyz,MPI_DOUBLE_PRECISION,nbrdown,nbrdown, comm2d,status,IERR)

     l=0
     do i=ngridstart(2), ngridstop(2)
        l=l+1
        do j=1, nwave
           coeff3d(ngridstop(1),i,j)=psitemp2dyz(l,j)
        end do
     end do


     CALL MPI_BARRIER(comm2d, IERR)

     ! 	write(6,*)'MYID=',MYID,' after first x-odd; passing/receiving.....'


     !------------------------------- 3D x-even ------------------------------------------------------	


     CALL KE_halfdt_3D_xeven

     CALL MPI_BARRIER(comm2d, IERR)

     ! 	write(6,*)'MYID=',MYID,' After x-even.....'

     !------------------------------------passing the boundary points---------------------------------------------------	

     m=0
     do i=ngridstart(2), ngridstop(2)
        m=m+1
        do j=1, nwave
           psitemp2dyz(m,j)=coeff3d(ngridstop(1),i,j)
        end do
     end do

     CALL MPI_SSEND(psitemp2dyz,npointyz,MPI_DOUBLE_PRECISION,nbrdown,MYID, comm2d, IERR)

     CALL MPI_RECV(psitemp2dyz,npointyz,MPI_DOUBLE_PRECISION,nbrup,nbrup, comm2d,status,IERR)

     l=0
     do i=ngridstart(2), ngridstop(2)
        l=l+1
        do j=1, nwave
           coeff3d(ngridstart(1),i,j)=psitemp2dyz(l,j)
        end do
     end do

     CALL MPI_BARRIER(comm2d, IERR)

     ! 	write(6,*)'MYID=',MYID,' After x-even; passing/receiving.....'

     !------------------------------- 3D x-odd AGAIN ------------------------------------------------	

     CALL KE_quarterdt_3D_xodd

     CALL MPI_BARRIER(comm2d, IERR)

     ! 	write(6,*)'MYID=',MYID,' After second x-odd.....'

     !------------------------------------passing the boundary points---------------------------------------------------	


     m=0
     do i=ngridstart(2), ngridstop(2)
        m=m+1
        do j=1, nwave
           psitemp2dyz(m,j)=coeff3d(ngridstart(1),i,j)
        end do
     end do

     CALL MPI_SSEND(psitemp2dyz,npointyz,MPI_DOUBLE_PRECISION,nbrup,MYID, comm2d, IERR)

     CALL MPI_RECV(psitemp2dyz,npointyz,MPI_DOUBLE_PRECISION,nbrdown,nbrdown, comm2d,status,IERR)

     l=0
     do i=ngridstart(2), ngridstop(2)
        l=l+1
        do j=1, nwave
           coeff3d(ngridstop(1),i,j)=psitemp2dyz(l,j)
        end do
     end do


     CALL MPI_BARRIER(comm2d, IERR)

     ! 	write(6,*)'MYID=',MYID,' After second x-odd; passing/receiving.....'


     !==============================================================================
     !	AGAIN Carrying out y-dimension operations
     !==============================================================================

     !------------------------------- 3D y-odd ------------------------------------------------------	

     CALL KE_quarterdt_3D_yodd

     CALL MPI_BARRIER(comm2d, IERR)

     ! 	write(6,*)'MYID=',MYID,' after first y-odd.....'

     !------------------------------------passing the boundary points---------------------------------------------------	


     m=0
     do i=ngridstart(1), ngridstop(1)
        m=m+1
        do j=1, nwave
           psitemp2dxz(m,j)=coeff3d(i,ngridstart(2),j)
        end do
     end do

     CALL MPI_SSEND(psitemp2dxz,npointxz,MPI_DOUBLE_PRECISION,nbrleft,MYID, comm2d, IERR)

     CALL MPI_RECV(psitemp2dxz,npointxz,MPI_DOUBLE_PRECISION,nbrright,nbrright, comm2d,status,IERR)

     l=0
     do i=ngridstart(1), ngridstop(1)
        l=l+1
        do j=1, nwave
           coeff3d(i,ngridstop(2),j)=psitemp2dxz(l,j)
        end do
     end do


     CALL MPI_BARRIER(comm2d, IERR)

     ! 	write(6,*)'MYID=',MYID,' after first y-odd; passing/receiving.....'


     !------------------------------- 3D y-even ------------------------------------------------------	


     CALL KE_halfdt_3D_yeven

     CALL MPI_BARRIER(comm2d, IERR)

     ! 	write(6,*)'MYID=',MYID,' after y-even.....'

     !------------------------------------passing the boundary points---------------------------------------------------	


     m=0
     do i=ngridstart(1), ngridstop(1)
        m=m+1
        do j=1, nwave
           psitemp2dxz(m,j)=coeff3d(i,ngridstop(2),j)
        end do
     end do

     CALL MPI_SSEND(psitemp2dxz,npointxz,MPI_DOUBLE_PRECISION,nbrright,MYID, comm2d, IERR)

     CALL MPI_RECV(psitemp2dxz,npointxz,MPI_DOUBLE_PRECISION,nbrleft,nbrleft, comm2d,status,IERR)

     l=0
     do i=ngridstart(1), ngridstop(1)
        l=l+1
        do j=1, nwave
           coeff3d(i,ngridstart(2),j)=psitemp2dxz(l,j)
        end do
     end do


     CALL MPI_BARRIER(comm2d, IERR)

     !	write(6,*)'MYID=',MYID,' after y-even; passing/receiving.....'

     !------------------------------- 3D y-odd AGAIN ------------------------------------------------	

     CALL KE_quarterdt_3D_yodd

     CALL MPI_BARRIER(comm2d, IERR)

     ! 	write(6,*)'MYID=',MYID,' after second y-odd.....'

     !------------------------------------passing the boundary points---------------------------------------------------	


     m=0
     do i=ngridstart(1), ngridstop(1)
        m=m+1
        do j=1, nwave
           psitemp2dxz(m,j)=coeff3d(i,ngridstart(2),j)
        end do
     end do

     CALL MPI_SSEND(psitemp2dxz,npointxz,MPI_DOUBLE_PRECISION,nbrleft,MYID, comm2d, IERR)

     CALL MPI_RECV(psitemp2dxz,npointxz,MPI_DOUBLE_PRECISION,nbrright,nbrright, comm2d,status,IERR)

     l=0
     do i=ngridstart(1), ngridstop(1)
        l=l+1
        do j=1,  nwave
           coeff3d(i,ngridstop(2),j)=psitemp2dxz(l,j)
        end do
     end do


     CALL MPI_BARRIER(comm2d, IERR)

     ! 	write(6,*)'MYID=',MYID,' after second y-odd; passing/receiving.....'


     !==============================================================================
     !==============================================================================
     !==============================================================================
     !---------Finishing one step (dt) evolution-----------------------------------
     !==============================================================================
     !==============================================================================
     !==============================================================================


     !*****************************************************************************
     !------------Building back the wave function----------------------------------
     !*****************************************************************************


     do k=ngridstart(1), ngridstop(1)
        do m=ngridstart(2), ngridstop(2)
           do n=1, nwave
              psi3d(k,m,n)=coeff3d(k,m,n)*factor(1,k)*factor(2,m)
           end do
        end do
     end do


     !*****************************************************************************
     !------------Normalizing the wave function for each time step-------------
     !*****************************************************************************


     temp=0.d0
     DO n=1, nwave
        DO j=regstart(1),regstop(1)
	   DO i=regstart(2), regstop(2)

              do k=1, num_fun(1,j)
                 do m=1, num_fun(2,i)
                    prob2d(k,m)=0.d0
                 end do
              end do

              if(j.eq.1)then
                 startj=2
              else
                 startj=1
              endif

              if(i.eq.1)then
                 starti=2
              else
                 starti=1
              endif

              do k=startj, num_fun(1,j)
                 do m=starti, num_fun(2,i)
                    prob2d(k,m)=psi3d(nindex(1,j,k),nindex(2,i,m),n)**2
                 end do
              end do

              CALL integral2d(num_fun(1,j),num_fun(2,i),mat_reg(1,j)%wt, &
                   mat_reg(2,i)%wt,prob2d,pnorm)
              temp=temp+pnorm

           END DO
        END DO
     END DO


     IF(MYID.ne.0)THEN
        CALL MPI_SSEND(temp,1,MPI_DOUBLE_PRECISION,0,MYID, comm2d, IERR)
     ELSE
        temp2=temp
        do j=1, nrcount
           CALL MPI_RECV(temp1,1,MPI_DOUBLE_PRECISION,j,j, comm2d, status, IERR)
           temp2=temp2+temp1
        end do
        temp3=1.d0/dsqrt(temp2)
     ENDIF

     CALL MPI_BARRIER(comm2d, IERR)
     CALL MPI_BCAST(temp3,1,MPI_DOUBLE_PRECISION,0,comm2d, IERR)
     CALL MPI_BARRIER(comm2d, IERR)


     !-----------------normalizing----------------------------
     do k=ngridstart(1), ngridstop(1)
        do m=ngridstart(2), ngridstop(2)
           do n=1, nwave
              coeff3d(k,m,n)=coeff3d(k,m,n)*temp3
              psi3d(k,m,n)=psi3d(k,m,n)*temp3
           end do
        end do
     end do


     CALL MPI_BARRIER(comm2d, IERR)


     !==============================================================================
     !--------------Calculating the <H>---------------------
     !==============================================================================


     IF ( mod(it,ntrec).eq.0 )THEN


        !--------------Calculating the <pr1^2/2> <pr2^2/2>-------------


        ctemp4=zero

        DO l=1, nwave
           DO i=ngstart(2),ngstop(2)
              DO j=regstart(1),regstop(1)

                 if(j.eq.1)then
                    startj=2
                 else
                    startj=1
                 endif

                 do k=startj, num_fun(1,j)
                    ctemp=zero
                    do n=startj, num_fun(1,j)
                       ctemp=ctemp+ mat_reg(1,j)%ke_mat(k,n)*coeff3d(nindex(1,j,n),i,l)
                    end do
                    ctemp_mat3(nindex(1,j,k))=ctemp
                 end do

                 ctemp=zero
                 do k=startj, num_fun(1,j)
                    ctemp= ctemp+ ( coeff3d(nindex(1,j,k),i,l)*ctemp_mat3(nindex(1,j,k)) )
                 end do

                 ctemp4=ctemp4+ctemp

              END DO
           END DO
        END DO

        !--------------------------------------------
        !--------------------------------------------
        ctemp5=zero

        DO l=1, nwave
           DO i=ngstart(1),ngstop(1)
              DO j=regstart(2),regstop(2)

                 if(j.eq.1)then
                    startj=2
                 else
                    startj=1
                 endif

                 do k=startj, num_fun(2,j)
                    ctemp=zero
                    do n=startj, num_fun(2,j)
                       ctemp=ctemp+ mat_reg(2,j)%ke_mat(k,n)*coeff3d(i,nindex(1,j,n),l)
                    end do
                    ctemp_mat3(nindex(2,j,k))=ctemp
                 end do

                 ctemp=zero
                 do k=startj, num_fun(2,j)
                    ctemp= ctemp+ ( coeff3d(i,nindex(1,j,k),l)*ctemp_mat3(nindex(2,j,k)) )
                 end do

                 ctemp5=ctemp5+ctemp

              END DO
           END DO
        END DO

        !--------------------------------------------

        IF(MYID.ne.0)THEN
	   CALL MPI_SSEND(ctemp4,1,MPI_DOUBLE_PRECISION,0,MYID, comm2d, IERR)
	   CALL MPI_SSEND(ctemp5,1,MPI_DOUBLE_PRECISION,0,MYID+ncount, comm2d, IERR)
        ELSE
	   do j=1, nrcount
              CALL MPI_RECV(ctemp1,1,MPI_DOUBLE_PRECISION,j,j, comm2d, status, IERR)
              CALL MPI_RECV(ctemp2,1,MPI_DOUBLE_PRECISION,j,j+ncount, comm2d, status, IERR)
              ctemp4=ctemp4+ctemp1
              ctemp5=ctemp5+ctemp2
           end do

           kinex= (ctemp4)
           kiney= (ctemp5)

        ENDIF


        !*************Calculating direct <pot>**********************************************

        temp=0.d0
        DO l=1, nwave
           DO j=regstart(1),regstop(1)
              DO i=regstart(2), regstop(2)

                 do k=1, num_fun(1,j)
                    do m=1, num_fun(2,i)
                       prob2d(k,m)=0.d0
                    end do
                 end do

                 if(j.eq.1)then
                    startj=2
                 else
                    startj=1
                 endif

                 if(i.eq.1)then
                    starti=2
                 else
                    starti=1
                 endif

                 do k=startj, num_fun(1,j)
                    do m=starti, num_fun(2,i)
                       prob2d(k,m) = psi3d(nindex(1,j,k),nindex(2,i,m),l)**2 &
                            &      * rpot3d(nindex(1,j,k),nindex(2,i,m),l)
                    end do
                 end do

                 CALL integral2d(num_fun(1,j),num_fun(2,i),mat_reg(1,j)%wt, &
                      mat_reg(2,i)%wt,prob2d,pnorm)
                 temp=temp+pnorm

              END DO
           END DO
        END DO


        IF(MYID.ne.0)THEN
	   CALL MPI_SSEND(temp,1,MPI_DOUBLE_PRECISION,0,MYID, comm2d, IERR)
        ELSE
	   pnorm=temp
	   do j=1, nrcount
              CALL MPI_RECV(temp1,1,MPI_DOUBLE_PRECISION,j,j, comm2d, status, IERR)
              pnorm=pnorm+temp1
           end do

	   pote = pnorm

        ENDIF


        !----------------------------------------------------------------------------------
        !*************Calculating 1/r12 <pot>**********************************************
        !----------------------------------------------------------------------------------

        temp=0.d0
        DO j=regstart(1),regstop(1)
	   DO i=regstart(2), regstop(2)

              do k=1, num_fun(1,j)
                 do m=1, num_fun(2,i)
                    prob2d(k,m)=0.d0
                 end do
              end do

              if(j.eq.1)then
                 startj=2
              else
                 startj=1
              endif

              if(i.eq.1)then
                 starti=2
              else
                 starti=1
              endif

              do k=startj, num_fun(1,j)
                 do m=starti, num_fun(2,i)

                    ir1=nindex(1,j,k)
                    ir2=nindex(2,i,m)

                    !------------------------------------------------------------------------

                    do n=1, nwave
                       ctemp1=zero
                       do l=1, nwave
                          ctemp1=ctemp1+VL(n,l,ir1,ir2)*psi3d(ir1,ir2,l)
                       end do
                       ctemp_mat3(n)=ctemp1
                    end do

                    ctemp1=zero
                    do l=1, nwave
                       ctemp1=ctemp1 + psi3d(ir1,ir2,l)*ctemp_mat3(l)
                    end do

                    prob2d(k,m)=( ctemp1 )

                    !------------------------------------------------------------------------

                 end do
              end do

              CALL integral2d(num_fun(1,j),num_fun(2,i),mat_reg(1,j)%wt, &
                   mat_reg(2,i)%wt,prob2d,pnorm)
              temp=temp+pnorm

	   END DO
        END DO

        !------------------------------------------------------------------------

        IF(MYID.ne.0)THEN
	   CALL MPI_SSEND(temp,1,MPI_DOUBLE_PRECISION,0,MYID, comm2d, IERR)
        ELSE
	   pnorm=temp
	   do j=1, nrcount
              CALL MPI_RECV(temp1,1,MPI_DOUBLE_PRECISION,j,j, comm2d, status, IERR)
              pnorm=pnorm+temp1
           end do

	   potee = pnorm

        ENDIF


9009    continue

        !------------------------------------------------------------------------


        CALL MPI_BARRIER(comm2d, IERR)


        !**********Recording the <H> and kinetic & potential energies***************

        IF( MYID.eq.0 ) THEN

           write(66, 802)time*convfs, kinex+kiney+pote+potee, kinex, kiney, pote, potee, kinex+kiney+pote
           call flush(66)

        ENDIF

        !--------------------------------------------------------------

     ENDIF

     !==============================================================================
     !--------------Calculating the norm, <r1> <r2> , and pop1---------------------
     !==============================================================================


     IF ( mod(it,ntrec).eq.0 )THEN
      
        
        temp=0.d0
        temp1=0.d0
        temp2=0.d0
        tmpwavenorms(:) = 0.d0
        DO n=1, nwave
           DO j=regstart(1),regstop(1)
              DO i=regstart(2), regstop(2)

                 do k=1, num_fun(1,j)
                    do m=1, num_fun(2,i)
                       prob2d(k,m)=0.d0
                       probx2d(k,m)=0.d0
                       proby2d(k,m)=0.d0
                    end do
                 end do

                 if(j.eq.1)then
                    startj=2
                 else
                    startj=1
                 endif

                 if(i.eq.1)then
                    starti=2
                 else
                    starti=1
                 endif

                 do k=startj, num_fun(1,j)
                    do m=starti, num_fun(2,i)
                       prob2d(k,m)= psi3d(nindex(1,j,k),nindex(2,i,m),n)**2
                       probx2d(k,m)=prob2d(k,m)*grid(1,nindex(1,j,k))
                       proby2d(k,m)=prob2d(k,m)*grid(2,nindex(2,i,m))
                    end do
                 end do

                 CALL integral2d(num_fun(1,j),num_fun(2,i),mat_reg(1,j)%wt, &
                      mat_reg(2,i)%wt,prob2d,pnorm)
                 CALL integral2d(num_fun(1,j),num_fun(2,i),mat_reg(1,j)%wt, &
                      mat_reg(2,i)%wt,probx2d,xexp)
                 CALL integral2d(num_fun(1,j),num_fun(2,i),mat_reg(1,j)%wt, &
                      mat_reg(2,i)%wt,proby2d,yexp)

                 tmpwavenorms(n) = tmpwavenorms(n) + pnorm
                 temp=temp+pnorm
                 temp1=temp1+xexp
                 temp2=temp2+yexp
                 
              END DO
           END DO
        END DO

        CALL MPI_REDUCE(tmpwavenorms, wavenorms, nwave, MPI_DOUBLE_PRECISION, MPI_SUM, 0, comm2d, IERR)
        CALL MPI_REDUCE(temp,  pnorm, 1,MPI_DOUBLE_PRECISION, MPI_SUM, 0, comm2d, IERR)
        CALL MPI_REDUCE(temp1,  xexp, 1,MPI_DOUBLE_PRECISION, MPI_SUM, 0, comm2d, IERR)
        CALL MPI_REDUCE(temp2,  yexp, 1,MPI_DOUBLE_PRECISION, MPI_SUM, 0, comm2d, IERR)
       
        IF(MYID.eq.0)THEN
           write(10,802) time*convfs, pnorm, wavenorms
           write(11,802) time*convfs, xexp
           write(12,802) time*convfs, yexp
           write(14,803) xexp, yexp
           call flush(10)	
           call flush(11)	
           call flush(12)	
           call flush(14)         
        ENDIF

        !--------calculating population left-------------------------------

        ctemp=zero
        DO l=1, nwave
           DO j=regstart(1),regstop(1)
              DO i=regstart(2), regstop(2)

                 do k=1, num_fun(1,j)
                    do m=1, num_fun(2,i)
                       cprob2d(k,m)=zero
                    end do
                 end do

                 if(j.eq.1)then
                    startj=2
                 else
                    startj=1
                 endif

                 if(i.eq.1)then
                    starti=2
                 else
                    starti=1
                 endif

                 do k=startj, num_fun(1,j)
                    do m=starti, num_fun(2,i)

                       cprob2d(k,m) = psi03d(nindex(1,j,k),nindex(2,i,m),l) &
                            &       * psi3d(nindex(1,j,k),nindex(2,i,m),l)
                    end do
                 end do

                 CALL integral2d(num_fun(1,j),num_fun(2,i),mat_reg(1,j)%wt, &
                      mat_reg(2,i)%wt,cprob2d,ctemp1)
                 ctemp=ctemp+ctemp1

              END DO
           END DO
        END DO


        IF(MYID.ne.0)THEN
	   CALL MPI_SSEND(ctemp,1,MPI_DOUBLE_PRECISION,0,MYID, comm2d, IERR)
        ELSE
	   do j=1, nrcount
              CALL MPI_RECV(ctemp1,1,MPI_DOUBLE_PRECISION,j,j, comm2d, status, IERR)
              ctemp=ctemp+ctemp1
           end do

	   temp = ctemp**2

           write(15,802) time*convfs, temp
           call flush(15)

        ENDIF


        !==============================================================================
        !--------------Calculating the <pr1> <pr2> -------------
        !==============================================================================


        ctemp4=zero
        ctemp5=zero

        DO l=1, nwave
           DO j=regstart(1),regstop(1)
              DO i=regstart(2), regstop(2)


                 do k=1, num_fun(1,j)
                    do m=1, num_fun(2,i)
                       cprobpx2d(k,m)=zero
                       cprobpy2d(k,m)=zero
                    end do
                 end do

                 if(j.eq.1)then
                    startj=2
                 else
                    startj=1
                 endif

                 if(i.eq.1)then
                    starti=2
                 else
                    starti=1
                 endif

                 do k=startj, num_fun(1,j)
                    do m=starti, num_fun(2,i)

                       ctemp=zero
                       do n=startj, num_fun(1,j)
                          ctemp=ctemp + mat_reg(1,j)%df(k,n)*psi3d(nindex(1,j,n),nindex(2,i,m),l)
                       end do
                       ctemp_mat4(k,m,1)=ctemp + psi3d(nindex(1,j,k),nindex(2,i,m),l)/grid(1,nindex(1,j,k))

                       ctemp=zero
                       do n=starti, num_fun(2,i)
                          ctemp=ctemp+ mat_reg(2,i)%df(m,n)*psi3d(nindex(1,j,k),nindex(2,i,n),l)
                       end do
                       ctemp_mat5(k,m,1)=ctemp + psi3d(nindex(1,j,k),nindex(2,i,m),l)/grid(2,nindex(2,i,m))

                    end do
                 end do


                 do k=startj, num_fun(1,j)
                    do m=starti, num_fun(2,i)
                       cprobpx2d(k,m)= psi3d(nindex(1,j,k),nindex(2,i,m),l)*ctemp_mat4(k,m,1)
                       cprobpy2d(k,m)= psi3d(nindex(1,j,k),nindex(2,i,m),l)*ctemp_mat5(k,m,1)
                    end do
                 end do

                 CALL integral2d(num_fun(1,j),num_fun(2,i),mat_reg(1,j)%wt,mat_reg(2,i)%wt,cprobpx2d,ctemp1)
                 CALL integral2d(num_fun(1,j),num_fun(2,i),mat_reg(1,j)%wt,mat_reg(2,i)%wt,cprobpy2d,ctemp2)

                 ctemp4=ctemp4+ctemp1
                 ctemp5=ctemp5+ctemp2

              END DO
           END DO
        END DO


        IF(MYID.ne.0)THEN
	   CALL MPI_SSEND(ctemp4,1,MPI_DOUBLE_PRECISION,0,MYID, comm2d, IERR)
	   CALL MPI_SSEND(ctemp5,1,MPI_DOUBLE_PRECISION,0,MYID+ncount, comm2d, IERR)
        ELSE
	   do j=1, nrcount
              CALL MPI_RECV(ctemp1,1,MPI_DOUBLE_PRECISION,j,j, comm2d, status, IERR)
              CALL MPI_RECV(ctemp2,1,MPI_DOUBLE_PRECISION,j,j+ncount, comm2d, status, IERR)
              ctemp4=ctemp4+ctemp1
              ctemp5=ctemp5+ctemp2
           end do

           write(111,802) time*convfs, ( ctemp4 )
           write(112,802) time*convfs, ( ctemp5 )
           write(114,803) ( ctemp4 ), ( ctemp5 )
           call flush(111)	
           call flush(112)	
           call flush(114)	

        ENDIF


        !-------------------------ntrec end---------------------------------

     ENDIF


     !==============================================================================
     !--------------Recording the temporal probability in two-dimension-------------
     !==============================================================================


     IF ( mod(it,ntprob).eq.0 )THEN

        npro=npro+1

        m=0
        do k=ngridstart(1), ngridstop(1)
           m=m+1
           l=0
           do n=ngridstart(2), ngridstop(2)
              l=l+1

              temp=0.d0
              do i=1, nwave
                 temp=temp+ psi3d(k,n,i)**2
              end do

              prob2dtemp(m,l)=temp
              if( prob2dtemp(m,l).le.1.d-45)prob2dtemp(m,l)=1.d-45

           end do
        end do


	IF(MYID.ne.0)then

	   CALL MPI_SSEND(ngridstart(1),1,MPI_INTEGER,0,MYID+2*ncount, comm2d, IERR)
	   CALL MPI_SSEND(ngridstart(2),1,MPI_INTEGER,0,MYID+3*ncount, comm2d, IERR)
	   CALL MPI_SSEND(prob2dtemp,npointxy,MPI_DOUBLE_PRECISION,0,MYID+4*ncount, comm2d, IERR)

	ELSE

           do k=ngridstart(1), ngridstop(1)
              do n=ngridstart(2), ngridstop(2)
                 probdens2d(k,n)=prob2dtemp(k,n)
              end do
           end do

	   do j=1, nrcount

              CALL MPI_RECV(ntemp1,1,MPI_INTEGER,j,j+2*ncount, comm2d, status, IERR)
              CALL MPI_RECV(ntemp2,1,MPI_INTEGER,j,j+3*ncount, comm2d, status, IERR)
              CALL MPI_RECV(prob2dtemp,npointxy,MPI_DOUBLE_PRECISION,j,j+4*ncount, comm2d, status, IERR)

              do k=1, npointx
                 do n=1, npointy
                    probdens2d(ntemp1+k-1,ntemp2+n-1)=prob2dtemp(k,n)
                 end do
              end do

	   end do


           !************Recording**********************************

           if(npro.eq.1)then
              open(80,file='prob_1.dat')
              write(80,*)'ZONE  T="',time*convfs,'", I=', ntot(2), ', J=', ntot(1),',  F=POINT'
              do i=1, ntot(1)
                 do j=1, ntot(2)
                    write(80,803)grid(1,i),grid(2,j),probdens2d(i,j)
                 end do
              end do
              close(80)
           endif

           if(npro.eq.2)then
              open(80,file='prob_2.dat')
              write(80,*)'ZONE  T="',time*convfs,'", I=', ntot(2), ', J=', ntot(1),',  F=POINT'
              do i=1, ntot(1)
                 do j=1, ntot(2)
                    write(80,803)grid(1,i),grid(2,j),probdens2d(i,j)
                 end do
              end do
              close(80)
           endif

           if(npro.eq.3)then
              open(80,file='prob_3.dat')
              write(80,*)'ZONE  T="',time*convfs,'", I=', ntot(2), ', J=', ntot(1),',  F=POINT'
              do i=1, ntot(1)
                 do j=1, ntot(2)
                    write(80,803)grid(1,i),grid(2,j),probdens2d(i,j)
                 end do
              end do
              close(80)
           endif

           if(npro.eq.4)then
              open(80,file='prob_4.dat')
              write(80,*)'ZONE  T="',time*convfs,'", I=', ntot(2), ', J=', ntot(1),',  F=POINT'
              do i=1, ntot(1)
                 do j=1, ntot(2)
                    write(80,803)grid(1,i),grid(2,j),probdens2d(i,j)
                 end do
              end do
              close(80)
           endif

           if(npro.eq.5)then
              open(80,file='prob_5.dat')
              write(80,*)'ZONE  T="',time*convfs,'", I=', ntot(2), ', J=', ntot(1),',  F=POINT'
              do i=1, ntot(1)
                 do j=1, ntot(2)
                    write(80,803)grid(1,i),grid(2,j),probdens2d(i,j)
                 end do
              end do
              close(80)
           endif

           if(npro.eq.6)then
              open(80,file='prob_6.dat')
              write(80,*)'ZONE  T="',time*convfs,'", I=', ntot(2), ', J=', ntot(1),',  F=POINT'
              do i=1, ntot(1)
                 do j=1, ntot(2)
                    write(80,803)grid(1,i),grid(2,j),probdens2d(i,j)
                 end do
              end do
              close(80)
           endif

           if(npro.eq.7)then
              open(80,file='prob_7.dat')
              write(80,*)'ZONE  T="',time*convfs,'", I=', ntot(2), ', J=', ntot(1),',  F=POINT'
              do i=1, ntot(1)
                 do j=1, ntot(2)
                    write(80,803)grid(1,i),grid(2,j),probdens2d(i,j)
                 end do
              end do
              close(80)
           endif

           if(npro.eq.8)then
              open(80,file='prob_8.dat')
              write(80,*)'ZONE  T="',time*convfs,'", I=', ntot(2), ', J=', ntot(1),',  F=POINT'
              do i=1, ntot(1)
                 do j=1, ntot(2)
                    write(80,803)grid(1,i),grid(2,j),probdens2d(i,j)
                 end do
              end do
              close(80)
           endif

           if(npro.eq.9)then
              open(80,file='prob_9.dat')
              write(80,*)'ZONE  T="',time*convfs,'", I=', ntot(2), ', J=', ntot(1),',  F=POINT'
              do i=1, ntot(1)
                 do j=1, ntot(2)
                    write(80,803)grid(1,i),grid(2,j),probdens2d(i,j)
                 end do
              end do
              close(80)
           endif

           if(npro.eq.10)then
              open(80,file='prob_10.dat')
              write(80,*)'ZONE  T="',time*convfs,'", I=', ntot(2), ', J=', ntot(1),',  F=POINT'
              do i=1, ntot(1)
                 do j=1, ntot(2)
                    write(80,803)grid(1,i),grid(2,j),probdens2d(i,j)
                 end do
              end do
              close(80)
           endif

           !************END of Recording**********************************

        ENDIF


     ENDIF
     

     !-------------Timing------------------------------------------

     CALL MPI_BARRIER(comm2d, IERR)

     ttx2=MPI_Wtime()

     IF(MYID.eq.0 .AND. mod(it,20).eq.0)THEN
        write(6,*)'it=',it,' The Total Time cost for this time-step =',ttx2-ttx1
        write(6,*)	
     ENDIF

     !--------------Ending of the time-loop-------------------------------------------

  END DO
  
  !-------------------------------------------------------

  CALL MPI_BARRIER(comm2d, IERR)

  endtime=MPI_Wtime()
  temp=(endtime-starttime)/dfloat(ntim)

  IF(MYID.eq.0 )THEN
     write(9,*)	
     write(9,*)'The total time cost for ', ntim,' steps Tot=',(endtime-starttime)/3600.0,' (hours)'
     write(9,*)'The average time cost per step is :=', temp,' (seconds)'
     write(9,*)	
     call flush(5)
  ENDIF


  CALL MPI_BARRIER(comm2d, IERR)


  !-------------------------------------------------------------------------------
  !-------------------------------------------------------------------------------
  !------Recording the final wave function----------------------------------------
  !-------------------------------------------------------------------------------
  !-------------------------------------------------------------------------------

  IF ( MYID.ne.0 ) THEN
     CALL MPI_SSEND(ngridstart(1),1,MPI_INTEGER,0,MYID       , comm2d, IERR)
     CALL MPI_SSEND(ngridstart(2),1,MPI_INTEGER,0,MYID+ncount, comm2d, IERR)
  ELSE
     do j=1, nrcount
        ntarget(j)=j*npointz
        CALL MPI_RECV(nst1(j),1,MPI_INTEGER,j,j       , comm2d,status,IERR)
        CALL MPI_RECV(nst2(j),1,MPI_INTEGER,j,j+ncount, comm2d,status,IERR)
     end do
  ENDIF

  !-------------------------------------------------------------------------------

  CALL MPI_BARRIER(comm2d, IERR)

  !-------------------------------------------------------------------------------
  IF ( MYID == 0 ) ALLOCATE ( wholepsi3d(ntot(1),ntot(2),nwave) )

  DO i=1, nwave

     IF ( MYID.ne.0 ) THEN

        j=0
        do m=ngridstart(1), ngridstop(1)
           j=j+1
           k=0
           do n=ngridstart(2), ngridstop(2)
              k=k+1
              psitemp2dxy(j,k)=psi3d(m,n,i)
           end do
        end do


        CALL MPI_SSEND(psitemp2dxy,npointxy,MPI_DOUBLE_PRECISION,0,MYID, comm2d, IERR)

     ELSE

        !*******************Receiving & Recording********************************

        do m=ngridstart(1), ngridstop(1)
           do n=ngridstart(2), ngridstop(2)
              psitemp2d(m,n)=psi3d(m,n,i)
           end do
        end do

        DO j=1, nrcount

           CALL MPI_RECV(psitemp2dxy,npointxy,MPI_DOUBLE_PRECISION,j, j, comm2d,status,IERR)

           do m=2, npointx
	      do n=2, npointy
                 psitemp2d(m+nst1(j)-1,n+nst2(j)-1)=psitemp2dxy(m,n)
	      end do
           end do

        END DO

        wholepsi3d(:,:,i) = psitemp2d(:,:)

        !write(46) psitemp2d(:,:) ! write binary data to GS.st

        !do n=1, ntot(2)
        !   do m=1, ntot(1)
        !      write(47,'(9999(G20.12E3,2x))') grid(1,m), grid(1,n), psitemp2d(m,n)
        !   end do
        !   write(47,*) ''
        !end do
        !write(47,*) ''; write(47,*) ''

        write(6,*)'3rd-dimension index i=',i,' Receiving the final wave function....'

     ENDIF

     CALL MPI_BARRIER(comm2d, IERR)


  END DO

  IF ( MYID == 0 ) THEN
     open(46,file='GS.st',status='unknown',form='unformatted')
     write(46) wholepsi3d(:,:,:) ! write binary data to GS.st
     close(46)

     open(47,file='groundstate.dat',status='unknown')
     write(47,'(2(A15,7x),9999(3x,I3.3,2x,I3.3,2x,I3.3,6x))') "#grid(x)", "grid(y)", (ltotmat(i),l1mat(i),l2mat(i),i=1,nwave)
     !                        3 + 3 + 2 + 3 + 2 + 3 + 6 = 22
     do n=1, ntot(2)
        do m=1, ntot(1)
           write(47,'(9999(G20.12E3,2x))') grid(1,m), grid(1,n), wholepsi3d(m,n,:)
        end do
        write(47,*) ''
     end do
     write(47,*) ''; write(47,*) ''
     close(47)

     DEALLOCATE ( wholepsi3d )
  END IF

  !----------------------------------------

  CALL MPI_BARRIER(comm2d, IERR)

  !----------------------------------------

  IF(MYID.eq.0)write(9,*)'Myid=',MYID,' ALL CALCULATIONS HAVE BEEN FINISHED !!'

  CALL MPI_FINALIZE(IERR)

END program Parallel_FEDVR_All
