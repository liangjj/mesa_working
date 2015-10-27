!***********************************************************************
! Read_FEDVR_Data_Module
!**begin prologue     Read_FEDVR_Data_Module
!**date written       090119   (yymmdd)
!**revision date               (yymmdd)
!**keywords           time, dvr, Iterative, Lanczos, propagation
!**
!**author             schneider, b. i.(nsf)
!**source             
!**purpose            Read in FEDVR Data
!***references
!***modules needed    See USE statements below
!***comments          
!***                  
!***                  
!***                  
!***                  
!***end prologue      Read_FEDVR_Data_Module
!***********************************************************************
!***********************************************************************
                           MODULE Read_FEDVR_Module
                           USE FEDVR_Global
  IMPLICIT NONE
!***********************************************************************
!***********************************************************************
                              CONTAINS
!***********************************************************************
!***********************************************************************
!deck Read_Grid_Parameters
!***begin prologue     Read_Grid_Parameters
!***date written       021231   (yymmdd)
!***revision date               (yymmdd)
!***keywords           dvr, basis, orthogonal polynomial, input
!***author             schneider, b. i.(nsf)
!***source             dvrlib
!***purpose            input subroutine for dvr basis sets
!***description        user interface for dvr library
!***references
!***routines called    
!***end prologue       Read_Grid_Parameters
  SUBROUTINE Read_Grid_Parameters(grid)
  IMPLICIT NONE
  TYPE(coordinates)                :: grid  
  CHARACTER(LEN=3)                 :: itoc
  INTEGER                          :: intkey
  INTEGER                          :: i
  INTEGER                          :: j
  INTEGER                          :: nply
  INTEGER                          :: nsubr
  INTEGER                          :: nblock
  INTEGER                          :: begin
  INTEGER                          :: ntest
  INTEGER                          :: n_sub
  LOGICAL                          :: logkey
  LOGICAL                          :: automte
  LOGICAL                          :: dollar
  REAL(idp)                        :: fpkey
  REAL(idp)                        :: step
  REAL(idp)                        :: boundl
  REAL(idp)                        :: boundr
  automte=logkey(card,'automate',.false.,' ')
  reuse_sector_information=logkey(card,'reuse_sector_information',.false.,' ')
  nfix=intkey(card,'number_of_fixed_points',2,' ')
  fix(1)=.false.
  fix(2)=.false.
  drop(1)=.false.
  drop(2)=.false.
  IF(nfix /= 0) THEN
     fix(1)=logkey(card,'left_fixed_point',.false.,' ')
     fix(2)=logkey(card,'right_fixed_point',.false.,' ')
     drop(1)=logkey(card,'drop_left_function',.false.,' ')
     drop(2)=logkey(card,'drop_right_function',.false.,' ')
  END IF
  bcl=1
  bcr=1
  IF(drop(1)) THEN
     bcl=0
  END IF
  IF(drop(2)) THEN
     bcr=0
  END IF
  IF(.not.automte) THEN
     nreg=intkey(card,'number_of_regions',1,' ')
     CALL fparr(card,'region_boundaries',edge,nreg+1,' ')
     CALL intarr(card,'polynomial_order_per_region',n,nreg,' ')
     npt=n+1
     nrq=npt 
     CALL intarr(card,'number_of_reference_quadrature_'//  &
                      'points_per_region',nrq,nreg,' ')
  ELSE
     write(iout,*)
     write(iout,*)
     write(iout,*) '                   Automated Selection of', &
                   ' Steps   '
     nreg=0
     nblock=intkey(card,'number_of_major_blocks',1,' ')       
     DO i=1,nblock 
        IF ( dollar('$'//grid%label//'_block_'//itoc(i),card, cpass,inp) ) THEN
             skip=logkey(card,'skip',.false.,' ')
             IF(skip) THEN
                write(iout,2) i
             ELSE 
                nply=intkey(card,'default_order',10,' ')
                nsubr=intkey(card,'number_of_subregions',1,' ')
                boundl=fpkey(card,'left_boundary',0.d0,' ')
                boundr=fpkey(card,'right_boundary',0.d0,' ')
                n_sub=intkey(card,'number_of_reference_quadrature_'//  &
                                  'points_per_sub_region',2*nply,' ')
                write(iout,1) i, nply, nsubr, boundl, boundr
                step = ( boundr - boundl ) / nsubr
                nreg = nreg + 1
                begin = nreg
                edge(nreg)=boundl
                DO j=1,nsubr
                   nrq(nreg) = n_sub
                   nreg = nreg + 1
                   edge(nreg) = edge(nreg-1) + step
                END DO
                nreg = nreg - 1
                n(begin:nreg)=nply
             END IF
        END IF
     END DO     
     npt = n + 1
  END IF
1   FORMAT(/,1x, 'block = ',i3, &
            /,15x,'polynomial order     = ',i4, &
            /,15x,'number of subregions = ',i4, & 
            /,15x,'left hand boundary   = ',e15.8, &
            /,15x,'right hand boundary  = ',e15.8)
2   FORMAT(/,1x,'skipping input block = ',i4)
  END SUBROUTINE Read_Grid_Parameters
!***********************************************************************
!***********************************************************************
!deck Read_Cartesian
!***begin prologue     Read_Cartesian
!***date written       021231   (yymmdd)
!***revision date               (yymmdd)
!***keywords           dvr, basis, orthogonal polynomial, input
!***author             schneider, b. i.(nsf)
!***source             dvrlib
!***purpose            
!***description        
!***                   
!***references
!***routines called    
!***end prologue       Read_Cartesian
  SUBROUTINE Read_Cartesian(grid)
  IMPLICIT NONE
  TYPE(coordinates)           :: grid  
  LOGICAL                     :: dollar
  LOGICAL                     :: logkey
  INTEGER                     :: intkey
  CHARACTER(LEN=80)           :: chrkey
  IF ( dollar('$cartesian_'//grid%label,card,cpass,inp) ) THEN
       nodiag=logkey(card,'do_not_diagonalize',.false.,' ')
       diag=logkey(card,'diagonalize',.true.,' ')
       typwt=chrkey(card,'weight_type','one',' ')
       unit_weight=logkey(card,'use_unit_weight',.false.,' ')
       atomic=logkey(card,'atomic_calculation',.false.,' ')
       typke='dvr'
       drctv=chrkey(card,'type_calculation','all_integrals',' ')
       if(drctv == 'poisson') then
          dentyp=chrkey(card,'type_density','exponential',' ')
       end if
       proj=.false.
       angmom=intkey(card,'angular_momentum',0,' ')
       reuse=logkey(card,'reuse_space_data',.false.,' ')
       Call Read_Grid_Parameters(grid)
       write(iout,1) nfix, fix, drop
!
!      Determine nphy
!
       call ptcal(physical_points,global_points,typwt)
       grid%num_fixed = nfix
       grid%fix_pt(1:2) = fix(1:2)
       grid%drop_pt(1:2) = drop(1:2)
       grid%num_reg = nreg
       ALLOCATE( grid%num_pts_reg(1:nreg) )
       grid%num_pts_reg(1:nreg) = npt(1:nreg)
!
  ELSE
       call lnkerr('no coordinate keyword found for label = '//grid%label)
  END If
1 Format(/,1x,'number of fixed points = ',i1,/15x,                   &
              'left point fixed       = ',l1,/,15x,                  &
              'right point fixed      = ',l1,/,15x,                  &
              'left point dropped     = ',l1,/15x,                   &
              'right point dropped    = ',l1)
  END SUBROUTINE Read_Cartesian
!***********************************************************************
!***********************************************************************
!deck Read_Spherical_Data
!***begin prologue     Read_Spherical_Data
!***date written       021231   (yymmdd)
!***revision date               (yymmdd)
!***keywords           dvr, basis, orthogonal polynomial, input
!***author             schneider, b. i.(nsf)
!***source             dvrlib
!***purpose            
!***description        
!***                   
!***references
!***routines called    
!***end prologue       Read_Spherical_Data
  SUBROUTINE Read_Spherical_Data
  IMPLICIT NONE
  LOGICAL                     :: dollar
  LOGICAL                     :: logkey
  INTEGER                     :: intkey
  CHARACTER(LEN=80)           :: chrkey
  IF ( dollar('$spherical',card,cpass,inp) ) THEN
       nodiag=logkey(card,'do_not_diagonalize',.false.,' ')
       diag=logkey(card,'diagonalize',.true.,' ')
       atomic=logkey(card,'atomic_calculation',.false.,' ')
       typke='dvr'
       drctv=chrkey(card,'type_calculation','all_integrals',' ')
       if(drctv == 'poisson') then
          dentyp=chrkey(card,'type_density','exponential',' ')
       end if
       proj=.false.
       l_max=intkey(card,'maximum_l_value',0,' ')
       m_max=intkey(card,'maximum_m_value',l_max,' ')
  ELSE
       call lnkerr('no keyword for spherical')
  END If
1 FORMAT(/,20x,'Basic Spherical Data',                              &
          /,5x,'maximum l value                  = ',i3,            &
          /,5x,'maximum m value                  = ',i4 )
  END SUBROUTINE Read_Spherical_Data
!***********************************************************************
!***********************************************************************
!deck Read_Spherical
!***begin prologue     Read_Spherical
!***date written       021231   (yymmdd)
!***revision date               (yymmdd)
!***keywords           dvr, basis, orthogonal polynomial, input
!***author             schneider, b. i.(nsf)
!***source             dvrlib
!***purpose            
!***description        
!***                   
!***references
!***routines called    
!***end prologue       Read_Spherical
  SUBROUTINE Read_Spherical(grid)
  IMPLICIT NONE
  TYPE(coordinates)           :: grid  
  LOGICAL                     :: dollar
  LOGICAL                     :: logkey
  INTEGER                     :: intkey
  CHARACTER(LEN=80)           :: chrkey
  IF ( dollar('$spherical_'//grid%label,card,cpass,inp) ) THEN
       typwt=chrkey(card,'weight_type','one',' ')
       unit_weight=logkey(card,'use_unit_weight',.false.,' ')
       typke='dvr'
       reuse=logkey(card,'reuse_space_data',.false.,' ')
       drctv=chrkey(card,'type_calculation','all_integrals',' ')
       if(drctv == 'poisson') then
          dentyp=chrkey(card,'type_density','exponential',' ')
       end if
       proj=.false.
       Call Read_Grid_Parameters(grid)
       write(iout,1) nfix, fix, drop
!
!      Determine nphy
!
       call ptcal(physical_points,global_points,typwt)
       grid%num_fixed = nfix
       grid%fix_pt(1:2) = fix(1:2)
       grid%drop_pt(1:2) = drop(1:2)
       grid%num_reg = nreg
       ALLOCATE( grid%num_pts_reg(1:nreg) )
       grid%num_pts_reg(1:nreg) = npt(1:nreg)
!
  ELSE
       call lnkerr('no coordinate keyword found for label = '//grid%label)
  END If
1 Format(/,1x,'number of fixed points = ',i1,/15x,                   &
              'left point fixed       = ',l1,/,15x,                  &
              'right point fixed      = ',l1,/,15x,                  &
              'left point dropped     = ',l1,/15x,                   &
              'right point dropped    = ',l1)
  END SUBROUTINE Read_Spherical
!***********************************************************************
!***********************************************************************
!deck Read_Fourier
!***begin prologue     Read_Fourier
!***date written       021231   (yymmdd)
!***revision date               (yymmdd)
!***keywords           dvr, basis, orthogonal polynomial, input
!***author             schneider, b. i.(nsf)
!***source             dvrlib
!***purpose            
!***description        
!***                   
!***references
!***routines called    
!***end prologue       Read_Fourier
  SUBROUTINE Read_Fourier(grid)
  IMPLICIT NONE
  TYPE(coordinates)                :: grid  
  INTEGER                          :: ntest
  LOGICAL                          :: dollar
  nreg=1
  nfix=0
  fix(1)=.false.
  fix(2)=.false.
  drop(1)=.false.
  drop(2)=.false.
  bcl=1
  bcr=1
  edge(1)=0.d0
  edge(2)=2.d0*pi
  IF ( dollar('$fourier',card,cpass,inp) ) THEN
       CALL fparr(card,'region_boundaries',edge,nreg+1,' ')
       CALL intarr(card,'number_of_points',n,nreg,' ')
  ELSE
       call lnkerr('no coordinate keyword found for label = fourier')
  END IF
!
!             if n is not odd fix it.
!
  ntest = n(1) - 2 * (n(1)/2)
  IF(ntest == 0 ) THEN
     n(1) = n(1) + 1
  END IF 
  npt(1)=n(1)
  nrq(1)=npt(1) 
  typwt='fourier'
  call ptcal(physical_points,global_points,'fourier')
  WRITE(iout,1) edge(1), edge(2)
  WRITE(iout,2) npt(1:nreg)
  grid%num_fixed = nfix
  grid%fix_pt(1:2) = fix(1:2)
  grid%drop_pt(1:2) = drop(1:2)
  grid%num_reg = nreg
  ALLOCATE( grid%num_pts_reg(1:nreg) )
  grid%num_pts_reg(1:nreg) = npt(1:nreg)
1 FORMAT(/,1x,'     Fourier : Region is (',e15.8,',',e15.8,')')
2 FORMAT(/,15x,'    Polynomial order          = ',                   &
                    (/,15x,5(i4,1x)))
  END SUBROUTINE Read_Fourier
!***********************************************************************
!***********************************************************************
!deck Read_Hermite
!***begin prologue     Read_Hermite
!***date written       021231   (yymmdd)
!***revision date               (yymmdd)
!***keywords           dvr, basis, orthogonal polynomial, input
!***author             schneider, b. i.(nsf)
!***source             dvrlib
!***purpose            input subroutine for dvr basis sets
!***description        user interface for dvr library
!***references
!***routines called    
!***end prologue       Read_Hermite
  SUBROUTINE Read_Hermite(grid)
  IMPLICIT NONE
  TYPE(coordinates)                     :: grid  
  INTEGER                               :: intkey
  LOGICAL                               :: logkey
  LOGICAL                               :: dollar
  write(iout,1)
  IF ( dollar('$hermite',card,cpass,inp) ) THEN
       nreg=1
       nfix=0
       fix(1)=.false.
       fix(2)=.false.
       drop(1)=.false.
       drop(2)=.false.
       CALL intarr(card,'number_of_points',n,nreg,' ')
       npt(1)=n(1)
  ELSE
       call lnkerr('no coordinate keyword found for label = hermite')
       Call lnkerr('no keyword found')
  END IF
  call ptcal(physical_points,global_points,'hermite')
  grid%num_fixed = nfix
  grid%fix_pt(1:2) = fix(1:2)
  grid%drop_pt(1:2) = drop(1:2)
  grid%num_reg = nreg
  ALLOCATE( grid%num_pts_reg(1:nreg) )
  grid%num_pts_reg(1:nreg) = npt(1:nreg)
1 FORMAT(/,1x,'     Hermite : Region is - infinity to + infinity')
!***********************************************************************
!***********************************************************************
  END SUBROUTINE Read_Hermite
!***********************************************************************
!***********************************************************************
!deck read_laguerre
!***begin prologue     read_laguerre
!***date written       021231   (yymmdd)
!***revision date               (yymmdd)
!***keywords           dvr, basis, orthogonal polynomial, input
!***author             schneider, b. i.(nsf)
!***source             dvrlib
!***purpose            input subroutine for dvr basis sets
!***description        user interface for dvr library
!***references
!***routines called    
!***end prologue       read_laguerre
  SUBROUTINE Read_Laguerre(grid)
  IMPLICIT NONE
  TYPE(coordinates)           :: grid  
  INTEGER                     :: intkey
  LOGICAL                     :: logkey
  nreg=1
  nfix=intkey(card,'number_of_fixed_points',0,' ')
  CALL intarr(card,'number_of_points',n,nreg,' ')
  IF(nfix == 1) THEN
     fix(1)=logkey(card,'left_fixed_point',.true.,' ')
     drop(1)=logkey(card,'drop_left_function',.false.,' ')
  END IF
  bcl=1
  IF(drop(1)) THEN
     bcl=0
  END IF
  npt(1)=n(1)
  call ptcal(physical_points,global_points,'laguerre')
  grid%num_fixed = nfix
  grid%fix_pt(1:2) = fix(1:2)
  grid%drop_pt(1:2) = drop(1:2)
  grid%num_reg = nreg
  ALLOCATE( grid%num_pts_reg(1:nreg) )
  grid%num_pts_reg(1:nreg) = npt(1:nreg)
!*************************************************************************************
  END SUBROUTINE Read_Laguerre
!*************************************************************************************
!*************************************************************************************
!deck Read_Legendre
!***begin prologue     Read_Legendre
!***date written       021231   (yymmdd)
!***revision date               (yymmdd)
!***keywords           dvr, basis, orthogonal polynomial, input
!***author             schneider, b. i.(nsf)
!***source             dvrlib
!***purpose            input subroutine for dvr basis sets
!***description        user interface for dvr library
!***references
!***routines called    
!***end prologue       Read_Legendre
  SUBROUTINE Read_Legendre(grid)
  IMPLICIT NONE
  TYPE(coordinates)                     :: grid  
  INTEGER                               :: intkey
  LOGICAL                               :: logkey
  LOGICAL                               :: dollar
  IF ( dollar('$legendre',card,cpass,inp) ) THEN
       m_val=intkey(card,'legendre_m',0,' ')
       skip=logkey(card,'read_grid_parameters',.false.,' ')
       write(iout,1) m_val
       nreg=1
!
!      This will do what needs to be done for one region
!      If there is more than one region the user needs to call
!      the grid_parameters subroutine and explicitly set things up.
!
       IF (skip) THEN
           grid%label='legendre'
           CALL Read_Grid_Parameters(grid)
       ELSE
           CALL intarr(card,'number_of_points',n,nreg,' ')
           npt(1)=n(1)
           nrq(1)=npt(1) 
           edge(1)=-1.d0
           edge(2)=1.d0
           IF(m_val == 0 ) THEN
              nfix=0
              fix(1)=.false.
              fix(2)=.false.
              drop(1)=.false.
              drop(2)=.false.
              bcl=1
              bcr=1
           ELSE 
              nfix=2
              fix(1)=.true.
              fix(2)=.true.
              drop(1)=.true.
              drop(2)=.true.
              bcl=0
              bcr=0
           END IF 
       END IF
       WRITE(iout,2) npt(1:nreg)
  ELSE
       call lnkerr('no coordinate keyword found for label = legendre')
  END IF
  call ptcal(physical_points,global_points,'legendre')
  grid%num_fixed = nfix
  grid%fix_pt(1:2) = fix(1:2)
  grid%drop_pt(1:2) = drop(1:2)
  grid%num_reg = nreg
  ALLOCATE( grid%num_pts_reg(1:nreg) )
  grid%num_pts_reg(1:nreg) = npt(1:nreg)
1 Format(/,1x,'legendre m value = ',i2)
2 FORMAT(/,15x,'    Polynomial order          = ',(/,15x,5(i4,1x)))
!***********************************************************************
!***********************************************************************
  END SUBROUTINE Read_Legendre
!***********************************************************************
!***********************************************************************
!deck Read_Theta
!***begin prologue     Read_Theta
!***date written       021231   (yymmdd)
!***revision date               (yymmdd)
!***keywords           dvr, basis, orthogonal polynomial, input
!***author             schneider, b. i.(nsf)
!***source             dvrlib
!***purpose            input subroutine for dvr basis sets
!***description        user interface for dvr library
!***references
!***routines called    
!***end prologue       Read_Theta
  SUBROUTINE Read_Theta(grid)
  IMPLICIT NONE
  TYPE(coordinates)                     :: grid  
  INTEGER                               :: intkey
  LOGICAL                               :: logkey
  LOGICAL                               :: dollar
  IF ( dollar('$theta',card,cpass,inp) ) THEN
       nreg=1
       m_val=intkey(card,'legendre_m',0,' ')
       write(iout,1) m_val
!
!     This will do what needs to be done for one region
!     If there is more than one region the user needs to call
!     the grid_parameters subroutine and explicitly set things up.
!
       skip=logkey(card,'read_grid_parameters',.false.,' ')
       IF (skip) THEN
           grid%label='theta'
           CALL Read_Grid_Parameters(grid)
       ELSE
           CALL intarr(card,'number_of_points',n,nreg,' ')
           npt(1)=n(1)
           nrq(1)=npt(1) 
           edge(1)=0.d0
           edge(2)=pi
           IF(m_val == 0 ) THEN
              nfix=0
              fix(1)=.false.
              fix(2)=.false.
              drop(1)=.false.
              drop(2)=.false.
              bcl=1
              bcr=1
           ELSE 
              nfix=2
              fix(1)=.true.
              fix(2)=.true.
              drop(1)=.true.
              drop(2)=.true.
              bcl=0
              bcr=0
           END IF 
       END IF
       WRITE(iout,2) npt(1:nreg)
  ELSE
       call lnkerr('no coordinate keyword found for label = theta')
  END IF
  call ptcal(physical_points,global_points,'theta')
  grid%num_fixed = nfix
  grid%fix_pt(1:2) = fix(1:2)
  grid%drop_pt(1:2) = drop(1:2)
  grid%num_reg = nreg
  ALLOCATE( grid%num_pts_reg(1:nreg) )
  grid%num_pts_reg(1:nreg) = npt(1:nreg)
1 Format(/,1x,'m value = ',i2)
2 FORMAT(/,15x,'    Polynomial order          = ',(/,15x,5(i4,1x)))
!***********************************************************************
!***********************************************************************
  END SUBROUTINE Read_Theta
!***********************************************************************
!***********************************************************************
!deck Read_Spheroidal
!***begin prologue     Read_Spheroidal
!***date written       021231   (yymmdd)
!***revision date               (yymmdd)
!***keywords           dvr, basis, orthogonal polynomial, input
!***author             schneider, b. i.(nsf)
!***source             dvrlib
!***purpose            input subroutine for spheroidal dvr basis sets
!***description        sets up the unperturbed dv hamiltonian for a two center
!***                   problem where charge Z_a is at -R/2 and Z_b at +R/2
!***references
!***routines called    
!***end prologue       Read_Spheroidal
  SUBROUTINE Read_Spheroidal(grid)
  IMPLICIT NONE
  TYPE(coordinates)                :: grid  
  LOGICAL                          :: dollar
  IF ( dollar('$spheroidal_'//grid%label,card,cpass,inp) ) THEN
!
       Call Read_Grid_Parameters(grid)
       write(iout,1) nfix, fix, drop
!
!      Determine nphy
!
       call ptcal(physical_points,global_points,'spheroidal')

!      Compute the Lobatto regional polynomials for even and odd types.
!
       grid%num_fixed = nfix
       grid%fix_pt(1:2) = fix(1:2)
       grid%drop_pt(1:2) = drop(1:2)
       grid%num_reg = nreg
       ALLOCATE( grid%num_pts_reg(1:nreg) )
       grid%num_pts_reg(1:nreg) = npt(1:nreg)
!
  ELSE
       call lnkerr('no coordinate keyword found for label = '//grid%label)
  END IF
1 Format(/,1x,'number of fixed points = ',i1,/15x,                   &
              'left point fixed       = ',l1,/,15x,                  &
              'right point fixed      = ',l1,/,15x,                  &
              'left point dropped     = ',l1,/15x,                   &
              'right point dropped    = ',l1)
  END SUBROUTINE Read_Spheroidal
!***********************************************************************
!***********************************************************************
!deck Read_Spheroidal_Data
!***begin prologue     Read_Spheroidal_Data
!***date written       021231   (yymmdd)
!***revision date               (yymmdd)
!***keywords           dvr, basis, orthogonal polynomial, input
!***author             schneider, b. i.(nsf)
!***source             dvrlib
!***purpose            input subroutine for spheroidal dvr basis sets
!***description        sets up the unperturbed dv hamiltonian for a two center
!***                   problem where charge Z_a is at -R/2 and Z_b at +R/2
!***references
!***routines called    
!***end prologue       Read_Spheroidal_Data
  SUBROUTINE Read_Spheroidal_Data
  IMPLICIT NONE
  INTEGER                          :: intkey
  LOGICAL                          :: dollar
  REAL(idp)                        :: fpkey
!
  IF ( dollar('$spheroidal',card,cpass,inp) ) THEN
!
!      z_a is located on the negative z axis, z_b on the positive.
!
       z_a=fpkey(card,'charge_nucleus_a',1.d0,' ')
       z_b=fpkey(card,'charge_nucleus_b',1.d0,' ')
       R_ab=fpkey(card,'internuclear_distance',1.d0,' ')
       l_max=intkey(card,'maximum_l_value',0,' ')
       m_max=intkey(card,'maximum_m_value',l_max,' ')
       WRITE(iout,1) z_a, z_b, R_ab, l_max, m_max
       typwt='one'
  ELSE
       Call lnkerr('no keyword found for spheroidal')
  END IF
1 FORMAT(/,20x,'Basic Spheroidal Data',                              &
          /,5x,'charge on nucleus a              = ',e15.8,          & 
          /,5x,'charge on nucleus b              = ',e15.8,          & 
          /,5x,'internuclear distance            = ',e15.8,          & 
          /,5x,'maximum l value                  = ',i3,             &
          /,5x,'maximum m value                  = ',i4 )
  END SUBROUTINE Read_Spheroidal_Data
!***********************************************************************
!***********************************************************************
!deck ptcal
!***begin prologue     ptcal
!***date written       021231   (yymmdd)
!***revision date               (yymmdd)
!***keywords           dvr, basis, orthogonal polynomial, input
!***author             schneider, b. i.(nsf)
!***source             dvrlib
!***purpose            determine the size of the gid with all constraints
!***description
!***references
!***routines called    
!***end prologue       ptcal
  SUBROUTINE ptcal(nphy,nglobal,typwt)
  IMPLICIT NONE
  INTEGER                    :: nphy
  INTEGER                    :: nglobal
  INTEGER                    :: i
  CHARACTER(LEN=*)           :: typwt
  IF(typwt == 'hermite') THEN
     nphy = npt(1)
     nglobal=nphy
  ELSE IF(typwt == 'laguerre') THEN
     nphy = npt(1)
     nglobal=nphy
     IF(bcl == 0) THEN
        nphy=nphy-1
     END IF
  ELSE IF(typwt == 'fourier') THEN
     nphy=npt(1)
     nglobal=nphy
  ELSE
     nphy=0
     DO  i=1,nreg
!
!        number of internal functions
! 
         nphy = nphy + npt(i) - 2
     END DO
! 
!     add one bridge function between each interval.
!
     nphy = nphy + nreg - 1
! 
!     add the extreme left and extreme right points
!     we have not yet dropped any functions at the endpoints.
!
     nphy=nphy + 2
     nglobal=nphy
     IF(bcl == 0) THEN
        nphy=nphy-1
     END IF
     IF(bcr == 0) THEN
        nphy=nphy-1
     END IF
  END IF
  Write(iout,1) nglobal, nphy
1 Format(/,20x,'Number of Global Grid Points   = ',i5,/,20x,   &
               'Number of Physical Grid Points = ',i5)
END SUBROUTINE ptcal
!***********************************************************************
!***********************************************************************
           END MODULE Read_FEDVR_Module
!***********************************************************************
!***********************************************************************
