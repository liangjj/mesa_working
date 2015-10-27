! \documentclass{article}
! \usepackage{graphicx}
! \setkeys{Gin}{width=\linewidth}
! \title{Main Program for Exact Propagator Code}
! \author{Barry I. Schneider}
! \date{}
! \def \<{\langle}
! \def \>{\rangle}
! \begin{document}
! \maketitle
!deck eigen_main
!**begin prologue     eigen_main
!**date written       960718   (yymmdd)
!**revision date               (yymmdd)
!**keywords           time, propagate, dvr, fd
!**
!**author             schneider, b. i.(nsf)
!**source             arnoldi
!**purpose            driver for solution of time dependent
!**                   schroedinger equation using exact propagation
!**description
!**references
!**routines called    iosys, util and mdutil
!**end prologue       eigen_main

  PROGRAM eigen_main
  USE arnoldi_global
  IMPLICIT NONE
  CHARACTER (LEN=1)                      :: itoc
  CHARACTER (LEN=80)                     :: chrkey
  LOGICAL                                :: dollar, logkey
  REAL*8                                 :: fpkey
  REAL*4                                 :: secnds
  INTEGER                                :: intkey, i, iostat, tim
  CHARACTER*4096                         :: ops 
  CALL Drum
  CALL IOsys('read character options from rwf',-1,0,0,ops)
  WRITE(iout,*)
  spdim=1
  IF ( dollar('$eigen_basis',card,cpass,inp) ) then
!
!         Set program options
!
        ntreg=intkey(card,'number-of-time-regions',1,' ')
        genpts=logkey(card,'automate-points',.false.,' ')
        system=chrkey(card,'coordinate-system','cartesian',' ')
        typke=chrkey(card,'kinetic-energy-type','dvr',' ')
        plot=logkey(card,'plot',.false.,' ')
        space=logkey(card,'no-spatial-hamiltonian',.false.,' ')
        diag=.true.
        imtime=logkey(card,'imaginary-time',.false.,' ')
        pr_main(1:2)='print=diag='//pr_main(1:2)
        pr_main(9)=chrkey(card,'print=diag=',pr_main(9),' ')
        IF(pr_main(9) == 'all') THEN
           CALL setprn(log_main,2)
        ELSE
           CALL setlog(log_main,pr_main,card,2)
        END IF
        IF(plot) THEN
           OPEN (UNIT=iplot(1),FILE='spatial_points', &
                 ACCESS='sequential',FORM='formatted', &
                 IOSTAT=IOSTAT,STATUS='unknown')
           OPEN (UNIT=iplot(2),FILE='time_points', &
                 ACCESS='sequential',FORM='formatted', &
                 IOSTAT=IOSTAT,STATUS='unknown')
           OPEN (UNIT=iplot(3),FILE='real_'//plot_name, &
                 ACCESS='sequential',FORM='formatted', &
                 IOSTAT=IOSTAT,STATUS='unknown')
           OPEN (UNIT=iplot(4),FILE='imaginary_'//plot_name, &
                 ACCESS='sequential',FORM='formatted', &
                 IOSTAT=IOSTAT,STATUS='unknown')
           OPEN (UNIT=iplot(5),FILE='absolute_'//plot_name, &
                 ACCESS='sequential',FORM='formatted', &
                 IOSTAT=IOSTAT,STATUS='unknown')
           OPEN (UNIT=iplot(6),FILE='auto_corr', &
                 ACCESS='sequential',FORM='formatted', &
                 IOSTAT=IOSTAT,STATUS='unknown')
           IF(IOSTAT /= 0) THEN
              CALL lnkerr('error in file handling')
           END IF
        END IF
        proj=logkey(card,'projections',.false.,' ')
        WRITE(iout,1) spdim
        WRITE(iout,2)
        WRITE(iout,3) space, ntreg
        DO i=1, spdim
           coord(i)=chrkey(card,'dimension-'//itoc(i),'x',' ')
        END DO   
        ALLOCATE(grid(spdim))
!
!       Get all of the one-dimensional matrices needed to construct
!       the spatial part of the hamiltonian and associated quantities.
!
        call iosys ('read character "bec filename" from rwf', &
                    -1,0,0,filbec)
        call iosys ('open bec as new',0,0,0,filbec)
        call iosys('rewind all on bec read-and-write',0,0,0,' ')
        time(1)=secnds(0.0)
        DO  i=1,spdim
            IF(typke == 'dvr'.OR.typke == 'packed') THEN
!
!          We are going to use a DVR basis.
!          Read the input file for the DVR data
!
              CALL dvr_input(nphy(i),nglobal(i),coord(i))
!
!             Allocate the needed DVR arrays.
!
              ALLOCATE(grid(i)%pt(nphy(i)), &
                       grid(i)%wt(nphy(i)), &
                       grid(i)%f(nphy(i),nphy(i)), &
                       grid(i)%df(nphy(i),nphy(i)), &    
                       grid(i)%ddf(nphy(i),nphy(i)), &    
                       grid(i)%ke(nphy(i),nphy(i)), &   
                       grid(i)%h(nphy(i),nphy(i)), &     
                       grid(i)%v(nphy(i)),grid(i)%srf_prm(2), &
                       grid(i)%eigv_0(nphy(i)), &    
                       grid(i)%eigvec_0(nphy(i),nphy(i)), &    
                       grid(i)%eigv(nphy(i)), &    
                       grid(i)%eigvec(nphy(i),nphy(i)), &    
                       grid(i)%srf_0(nphy(i),2), &
                       grid(i)%srf(nphy(i),2))    
!
!             Compute the DVR points, weights, functions, first and second
!             derivatives, kinetic energy matrix, eigenvalues and eigenvectors
!             of the kinetic energy matrix, full one-particle Hamiltonian
!             matrix,eigenvalues and eigenvectors of the Hamiltonian matrix, 
!             one-body potential, value of DVR functions at the endpoints, 
!             and value of eigenvectors of the kinetic energy and one-body 
!             Hamiltonian at the endpoints.
!
              CALL dvr_basis(pt_0(i),grid(i)%pt,grid(i)%wt, &
                             grid(i)%f,grid(i)%df, &
                             grid(i)%ddf,grid(i)%ke, &
                             grid(i)%eigv_0,grid(i)%eigvec_0, &
                             grid(i)%h,grid(i)%eigv,grid(i)%eigvec, &
                             grid(i)%v,grid(i)%srf_prm,grid(i)%srf_0, &
                             grid(i)%srf,coord(i),nphy(i),nglobal(i))
              row(i)=2*nphy(i) - 2
!
!             We are going to use either a 3, 5, or 7 point finite 
!             difference approximation to the second derivative.
!
            ELSE IF(typke == 'fd') THEN
!
!             Read the finite difference input
!
              CALL fd_input(nphy(i),nglobal(i),row(i),coord(i))
!
!             Allocate the needed memory for the points, weights, 
!             kinetic energy matrix, and potential.
!
              ALLOCATE(grid(i)%pt(nphy(i)), &
                       grid(i)%wt(nphy(i)), &
                       grid(i)%ke(row(i),nphy(i)), &    
                       grid(i)%h(row(i),nphy(i)), &    
                       grid(i)%v(nphy(i)), & 
                       grid(i)%eigv_0(nphy(i)), &    
                       grid(i)%eigvec_0(nphy(i),nphy(i)), &    
                       grid(i)%eigv(nphy(i)), &    
                       grid(i)%eigvec(nphy(i),nphy(i)))
!
!             Calculate the finite difference quantities.
!
              CALL fd_basis(pt_0(i),grid(i)%pt,grid(i)%wt, &
                            grid(i)%ke,grid(i)%eigv_0,grid(i)%eigvec_0, &
                            grid(i)%h,grid(i)%eigv,grid(i)%eigvec, &
                            grid(i)%v,nphy(i),nglobal(i),row(i),coord(i))
            ELSE
              WRITE(iout,4)
              STOP
            END IF
            time(2)=secnds(0.0)
            delta(1)=time(2)-time(1)
            WRITE(iout,5) delta(1)
        END DO
  ELSE
       write(iout,6)
       stop
  END IF
!
!               We now assume we have available the
!               eigenvalues and vectors of the hamiltonian.
!
!
!               Now do the actual Propagation
!
  IF( dollar('$time',card,cpass,inp) ) then
!  
!     Read in the number of time regions and their edges
!  
      DO  i=3,8
          pr_main(i)='print=prop='//pr_main(i)
      END DO
      pr_main(9)=chrkey(card,'print=prop=',pr_main(9),' ')
      IF(pr_main(9) == 'all') THEN
         CALL setprn(log_main(3),6)
      ELSE
         CALL setlog(log_main(3),pr_main(3),card,6)
      END IF
      ntreg=intkey(card,'number-of-time-regions',1,' ')
      ALLOCATE(edge(ntreg))
      IF(logkey(card,'automate',.false.,' ')) THEN
         edge(1)=fpkey(card,'first-time',0.d0,' ')
         delt=fpkey(card,'time-interval',.01D0,' ')
         DO  tim=2,ntreg+1
             edge(tim) = edge(tim-1) + delt
         END DO
      ELSE
         CALL fparr(card,'time-points',edge,ntreg+1,' ')
      END IF
      WRITE(iout,7) ntreg, (edge(i),i=1,ntreg+1)
!
!     This is the driving routine for the propagation
!
      CALL eigen_prop
  ELSE
      write (iout,8)
      stop
  END IF
  IF(plot) THEN
      CLOSE (UNIT=iplot(1),IOSTAT=IOSTAT)
      CLOSE (UNIT=iplot(2),IOSTAT=IOSTAT)
      CLOSE (UNIT=iplot(3),IOSTAT=IOSTAT)
      CLOSE (UNIT=iplot(4),IOSTAT=IOSTAT)
      CLOSE (UNIT=iplot(5),IOSTAT=IOSTAT)
      CLOSE (UNIT=iplot(6),IOSTAT=IOSTAT)
      IF(IOSTAT /= 0) THEN
         write(iout,9)
         stop 
      END IF
  END IF
!
!                   Deallocate all of the memory and quit.
!
  do i=1,spdim
     IF(typke == 'dvr'.OR.typke == 'packed') THEN
        DEALLOCATE(grid(i)%pt,grid(i)%wt,grid(i)%f,           &
                   grid(i)%df,grid(i)%ddf,grid(i)%ke,         &
                   grid(i)%h,grid(i)%v,grid(i)%srf_prm, &
                   grid(i)%eigv_0,grid(i)%eigvec_0, &
                   grid(i)%eigv,grid(i)%eigvec,     &
                   grid(i)%srf_0,grid(i)%srf  )    
     ELSE IF(typke == 'fd') THEN
        DEALLOCATE(grid(i)%pt,grid(i)%wt,grid(i)%ke,grid(i)%v, &    
                   grid(i)%eigv_0,grid(i)%eigvec_0, &    
                   grid(i)%eigv,grid(i)%eigvec)
     END IF
  END DO
  DEALLOCATE(grid)
  call chainx(0)
  stop
1    FORMAT(/,20X,'time-dependent basis function code',//,20X,  &
                  'number of spatial dimensions = ',i1)
2    FORMAT(/,15X,'calculation = solve time-dependent schrodinger'  &
                  ' equation')
3    FORMAT(/,5X,'time-dependent data',/,5X,  &
                 'no spatial hamiltonian   = ',l1, &
            /,5X, 'number of time intervals = ',i3)
4    FORMAT(/,5x,'basis type error')
5    FORMAT('***********************************************' &
            '*************************'                       &
            /,10X,'time to compute the spatial Hamiltonian '  &
            /,10x,'and associated quantities = ',f15.8,/,     &
            '***********************************************' &
            '*************************')
6    FORMAT(/,5x,'no basis card section')
7    FORMAT(/,5X,'number of time intervals = ',i4,/,5X,  &
                 'time intervals = ',(/,5X,5F15.8) )
8    FORMAT(/,5x,'no propagation card section')
9    FORMAT(/,1x,'error in file handling')
END PROGRAM eigen_main
