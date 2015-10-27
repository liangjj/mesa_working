! \documentclass{article}
! \usepackage{graphicx}
! \setkeys{Gin}{width=\linewidth}
! \title{Main Program for One Dimensional R-matrix Code}
! \author{Barry I. Schneider}
! \date{}
! \def \<{\langle}
! \def \>{\rangle}
! \begin{document}
! \maketitle
!deck r_matrix_main
!**begin prologue     r_matrix_main
!**date written       030525   (yymmdd)
!**revision date               (yymmdd)
!**keywords           r-matrix, dvr
!**
!**author             schneider, b. i.(nsf)
!**source             R_Matrix_f90
!**purpose            r-matrix calculation of phase shifts for a one
!**                   dimensional schroedinger equation using drv basis sets.
!**description
!**references
!**routines called    iosys, util and mdutil
!**end prologue       r_matrix_main
  PROGRAM r_matrix_main
  USE io
  USE r_matrix_global
  USE dvr_shared
  USE dvr_global
  IMPLICIT NONE
  CHARACTER (LEN=80)                     :: chrkey
  LOGICAL                                :: dollar, logkey, set
  REAL*8                                 :: fpkey, elow, ehigh, dele, ratio, tmp
  REAL*4                                 :: secnds
  INTEGER                                :: intkey, i, acc=30
  CHARACTER*4096                         :: ops 
  CALL Drum
  CALL IOsys('read character options from rwf',-1,0,0,ops)
  pi=acos(-1.d0)
  WRITE(iout,*)
  open(unit=99,file='phase_shift',status='unknown')
  write(99,4) '$res'
  IF ( dollar('$r_matrix_basis',card,cpass,inp) ) then
!
!         Set program options
!
        system=chrkey(card,'coordinate-system','cartesian',' ')
        typke=chrkey(card,'kinetic-energy-type','dvr',' ')
        diag=logkey(card,'get-eigenpairs',.false.,' ')
        pr_main='print=r-matrix='//pr_main
        pr_main(6)=chrkey(card,'print=r-matrix=',pr_main(6),' ')
        IF(pr_main(6) == 'all') THEN
           CALL setprn(log_main,5)
        ELSE
           CALL setlog(log_main,pr_main,card,5)
        END IF
        coord(1)=chrkey(card,'dimension','x',' ')
        set=logkey(card,'explicit-read',.false.,' ')
        eunit=chrkey(card,'energy-units','k**2/2',' ')
        if(set) then
           nen=intkey(card,'number-of-energies',1,' ')
           ALLOCATE(energy(nen))
           call fparr(card,'energies',energy,nen,' ')   
        else
           elow=fpkey(card,'lowest-energy',.001,' ')
           ehigh=fpkey(card,'highest-energy',5.0,' ')
           dele=fpkey(card,'energy-step',.001,' ')
           nen=(ehigh-elow)/dele + 1
           ALLOCATE(energy(nen))
           energy(1)=elow
           do i=2,nen
              energy(i) = energy(i-1) + dele
           end do
        endif
        IF(eunit == 'k') then
           energy = energy*energy*.5d0
        END IF 
        angmom=intkey(card,'angular-momentum',0,' ')
        WRITE(iout,1) system, coord(1), eunit, typke, diag, &
                      nen, angmom 
        WRITE(iout,2) energy(1:nen)
!
!       Get all of the one-dimensional matrices needed to construct
!       the spatial part of the hamiltonian and associated quantities.
!
!
        ALLOCATE(grid(1))
        CALL dvr_input(nphy(1),nglobal(1),coord(1))
!
!       Allocate the needed DVR arrays.
!
         ALLOCATE(grid(1)%pt(nphy(1)),            &
                  grid(1)%wt(nphy(1)),            &
                  grid(1)%f(nphy(1),nphy(1)),     &
                  grid(1)%df(nphy(1),nphy(1)),    &
                  grid(1)%ddf(nphy(1),nphy(1)),   &
                  grid(1)%ke(nphy(1),nphy(1)),    &
                  grid(1)%p_mom(nphy(1),nphy(1)), &
                  grid(1)%h(nphy(1),nphy(1)),     &
                  grid(1)%v(nphy(1)),grid(1)%srf_prm(2))
         IF(diag) then
            ALLOCATE(grid(1)%eigv_0(nphy(1)),           &
                     grid(1)%eigvec_0(nphy(1),nphy(1)), &
                     grid(1)%eigv(nphy(1)),             &
                     grid(1)%eigvec(nphy(1),nphy(1)),   &
                     grid(1)%srf_0(nphy(1),2),          &
                     grid(1)%srf(nphy(1),2))
         END IF
!
!
!             Compute the DVR points, weights, functions, first and second
!             derivatives, kinetic energy matrix, eigenvalues and eigenvectors
!             of the kinetic energy matrix, full one-particle Hamiltonian
!             matrix,eigenvalues and eigenvectors of the Hamiltonian matrix, 
!             one-body potential, value of DVR functions at the endpoints, 
!             and value of eigenvectors of the kinetic energy and one-body 
!             Hamiltonian at the endpoints.
!

         call dvr_basis(pt_0(1),grid(1)%pt,grid(1)%wt,grid(1)%f,grid(1)%df,  &
                        grid(1)%ddf,grid(1)%ke,grid(1)%p_mom,grid(1)%eigv_0, &
                        grid(1)%eigvec_0,grid(1)%h,grid(1)%eigv,             &
                        grid(1)%eigvec,grid(1)%v,grid(1)%srf_prm,            &
                        grid(1)%srf_0,grid(1)%srf,coord(1),nphy(1),nglobal(1)) 
        row(1)=2*nphy(1) - 2
!
  ELSE
       write(iout,3)
  END IF
!
!             Calculate phase shifts and cross section information
!
  tmp=angmom*acc
  ltop=angmom+sqrt(tmp)
  ltop=max(ltop,angmom)
  ALLOCATE(phase(nen),kmat(nen),jbes(0:ltop),djbes(0:ltop), &
           ybes(0:ltop),dybes(0:ltop))
  DO i=1,nen
     call r_mat(i)
  END DO
  write(99,5) energy
  write(99,6) phase
  write(99,4) "$"
  phase=phase+pi
  write(99,6) phase
!
!                   Deallocate all of the memory and quit.
!
  DEALLOCATE(grid(1)%pt,                     &
             grid(1)%wt,                     &
             grid(1)%f,                      &
             grid(1)%df,                     &
             grid(1)%ddf,                    &
             grid(1)%ke,                     &
             grid(1)%p_mom,                  &
             grid(1)%h,                      &
             grid(1)%v,                      &
             grid(1)%srf_prm)
  IF(diag) then
     DEALLOCATE(grid(1)%eigv_0,              &
                grid(1)%eigvec_0,            &
                grid(1)%eigv,                &
                grid(1)%eigvec,              &
                grid(1)%srf_0,               &
                grid(1)%srf)
  END IF
  call chainx(0)
  stop
1    FORMAT(/,20X,'one dimensional r-matrix code', &
            /,10x,'coordinate system    = ',a32, &
            /,10x,'coordinate           = ',a32, &
            /,10x,'energy units         = ',a10, &
            /,10x,'kinetic energy basis = ',a16, &
            /,10x,'diagonalize          = ',l1,  &
            /,10x,'number of energies   = ',i4,  &
            /,10x,'angular momentum     = ',i2)
2    FORMAT(/,'energies = ',/,(4(1x,e15.8,1x)))
3    FORMAT(/,5x,'no basis card section')
4    FORMAT(a4)
5    FORMAT('elist=',5(e15.8,','))
6    FORMAT('del=',5(e15.8,','))
END PROGRAM r_matrix_main
