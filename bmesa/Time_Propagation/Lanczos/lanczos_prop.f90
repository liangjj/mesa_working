!deck lanczos_prop
!**begin prologue     lanczos_prop
!**date written       960718   (yymmdd)
!**revision date               (yymmdd)
!**keywords           time, dvr, Arnoldi, propagation
!**
!**author             schneider, b. i.(nsf)
!**source             Lanczos
!**purpose            time dependent schroedinger equation
!**                   using Lanczos method with finite difference or
!**                   dvr space representation.
!**references
!**routines called    iosys, util and mdutil
!**end prologue       lanczos_prop
  SUBROUTINE lanczos_prop
  USE lanczos_global
  USE dvr_shared
  USE dvr_global
  IMPLICIT NONE
  CHARACTER (LEN=16)                      :: fptoc
  CHARACTER (LEN=2)                       :: itoc
  CHARACTER(LEN=8)                        :: cntrl 
  INTEGER                                 :: i, bigs, bigv, t
  INTEGER                                 :: iostat, length, len
  INTEGER                                 :: it
  COMPLEX*16                              :: conj
  REAL*8, DIMENSION(:), ALLOCATABLE       :: rtemp
  REAL*8                                  :: error, t_int, total_t
  REAL*4                                  :: secnds
!
  time(1)=secnds(0.0)
!
!     Pack the Hamiltonian into hbuf and ihbuf for each
!     dimension if requested.
!
!     Allocate memory for the variable type buf and the associated arrays.
!
  ALLOCATE(buf(spdim))
  IF(typke == 'packed') THEN
     ALLOCATE( buf(1)%d(nphy(1)), &
               buf(1)%hbuf(nphy(1)*nphy(1)), &
               buf(1)%hibuf(2,nphy(1)*nphy(1)))
     CALL pack_h(1)
     IF(spdim > 1) THEN
        ALLOCATE(buf(2)%d(nphy(2)), &
                 buf(2)%hbuf(nphy(2)*nphy(2)), &
                 buf(2)%hibuf(2,nphy(2)*nphy(2)))
        CALL pack_h(2)
        IF(spdim > 2) THEN
           ALLOCATE(buf(3)%d(nphy(3)), &
                    buf(3)%hbuf(nphy(3)*nphy(3)), &
                    buf(3)%hibuf(3,nphy(2)*nphy(2)))
           CALL pack_h(3)
        END IF
     END IF
  END IF
  time(2)=secnds(0.0)
  delta(1)=time(2)-time(1)
  WRITE(iout,1) delta(1) 
!
!      Get total dimensionality in space.  nphy contains the 
!      information for each dimension.
!
  n3d=1
  DO  i=1,spdim
      n3d=n3d*nphy(i)
  END DO
!
! Read in any data for the arnoldi iterations.
!
  CALL lanczos_dat
!
!              Allocate memory for all of the major variables needed in the
!              propagation.  The variable n3d is the length of the vector and
!              maxvec is the maximum size of the Arnoldi space.
!
  ALLOCATE(psi0(n3d),soln_0(n3d),chi(n3d),v_tot(n3d),vec(n3d,0:maxit), &
           hvec(n3d),a_lanczos(maxit),b_lanczos(maxit),bwrk(maxit), &
           eig(maxit),b(maxit,maxit),u(n3d,maxit),vscr(n3d),work(n3d*maxit), &
           tim_pts(ntreg+1),auto_corr(ntreg),rtemp(max(n3d,ntreg)))
!
!     Set up a file to hold the wavefunction for plotting
!
  call iosys('open plot as scratch',0,0,0,' ')
  call iosys('create real wavefunction on plot',ntreg*2*n3d,0,0,' ')
  DO t=1,ntreg+1
     tim_pts(t)=edge(t)
  END DO
  IF(plot) then
     DO i=1,spdim
        write(iplot(1),*) grid(i)%pt
     END DO
     write(iplot(2),*) (tim_pts(i),i=1,ntreg+1)
  END IF
  DO  t=1,ntreg
      keywrd=itoc(t)
      LEN=length(keywrd)
      keywrd='t'//keywrd(1:LEN)
      LEN=length(keywrd)
      keywrd='$v0('//keywrd(1:LEN)//')'
      t0=edge(t)
      t1=edge(t+1)
      deltat=t1-t0
      time(1)=secnds(0.0)
      v_tot = 0.d0
      CALL v_tim
      CALL pert
      time(2)=secnds(0.0)
      delta(1)=time(2)-time(1)
      CALL cp_psi0(t)
      time(3)=secnds(0.0)
      delta(2)=time(3)-time(2)
      IF(vtyp(2) /= 'none') THEN
         CALL v_nl
      END IF
      time(4)=secnds(0.0)
      delta(3)=time(4)-time(3)
! 
!     initialize the first lanczos vector
!
      error = 1.d+10
      write(iout,2) error
      cntrl='continue'
      vec(:,0) = psi0(:)     
      it=0
      del_t(1:4)=0.e0
      DO WHILE ( cntrl == 'continue'.AND.it < maxit )
!        
         it = it + 1
!
!
!        perform a lanczos iteration
!
         eltim(1)=secnds(0.0)
         CALL lanczos(it)
         eltim(2)=secnds(0.0)
         del_t(1) = del_t(1) + eltim(2) - eltim(1)
         IF(it > maxit) THEN
             WRITE(iout,3) (del_t(i),i=1,4)
             WRITE(iout,4)
             RETURN
         END IF
         WRITE(iout,5) it
         cntrl='continue'
!
!        calculate the eigenvalues of the tridiagonal matrix
!       
         call  tri_diag(it)
         eltim(3) = secnds(0.0)
         del_t(2) = del_t(2) + eltim(3) - eltim(2)
         CALL ecbcx(u,n3d,vec,n3d,b,maxit,n3d,it,it) 
         eltim(4) = secnds(0.0)
         del_t(3) = del_t(3) + eltim(4) - eltim(3)
         CALL exp_prop('subtract',it)
         eltim(5) = secnds(0.0)
         del_t(4) = del_t(4) + eltim(5) - eltim(4)
         IF(it == 1) THEN
            vscr=chi
         ELSE
            CALL chk_con(cntrl)
         END IF
         IF(cntrl == 'finished') THEN
            WRITE(iout,3) (del_t(i),i=1,4)
            WRITE(iout,6)
         END IF
         IF(it >= maxit) THEN
            WRITE(iout,3) (del_t(i),i=1,5)
            WRITE(iout,7)
            RETURN
         END IF
      END DO
!
!        compare the approximate and exact solutions where possible
!
     IF(i0stat == 'gaussian-pulse') THEN
        CALL moment
     END IF
! Form the total solution and compute the autocorrelation
! function.
     CALL soln(t,rtemp)
     time(5)=secnds(0.0)
     delta(4)=time(5)-time(4)
     total_t = delta(1) + delta(2) + delta(3) + delta(4)
     WRITE(iout,8) t, (delta(i), i=1,4), total_t
     call chk_nrm(chi,n3d)
  END DO
!
!     End the propagation
!
   do i=1,ntreg
      rtemp(i) = 1.d0 - auto_corr(i) * conjg(auto_corr(i))
   END DO
   write(iout,9)
   do t=1,ntreg,plot_step
      write(iout,11) edge(t), rtemp(t)
   end do
   CALL iosys('rewind all on plot read-and-write',0,0,0,' ')
   IF(plot) then
      write(iplot(6),*) (rtemp(i),i=1,ntreg-1)
   END IF
!     Release the memory
!
   DEALLOCATE(psi0,soln_0,chi,v_tot,vec,hvec,a_lanczos,b_lanczos,bwrk, &
              eig,b,u,vscr,work,tim_pts,auto_corr,rtemp)
  IF(typke == 'packed') THEN
     DEALLOCATE( buf(1)%d,buf(1)%hbuf,buf(1)%hibuf)
     IF(spdim > 1) THEN
        DEALLOCATE(buf(2)%d,buf(2)%hbuf,buf(2)%hibuf)
        IF(spdim > 2) THEN
           DEALLOCATE(buf(3)%d,buf(3)%hbuf,buf(3)%hibuf)
        END IF
     END IF
  END IF
  DEALLOCATE(buf)
  do i=1,spdim
     IF(typke == 'dvr'.OR.typke == 'packed') THEN
        DEALLOCATE(grid(i)%pt,grid(i)%wt,grid(i)%f,           &
                   grid(i)%df,grid(i)%ddf,grid(i)%ke,         &
                   grid(i)%p_mom,grid(i)%h,grid(i)%v, &
                   grid(i)%srf_prm)
        if(diag) then          
           DEALLOCATE(grid(i)%eigv_0,grid(i)%eigvec_0, &
                      grid(i)%eigv,grid(i)%eigvec,     &
                      grid(i)%srf_0,grid(i)%srf  )    
        endif
     ELSE IF(typke == 'fd') THEN
        DEALLOCATE(grid(i)%pt,grid(i)%wt,grid(i)%ke,grid(i)%v)    
        if(diag) then
           DEALLOCATE(grid(i)%eigv_0,grid(i)%eigvec_0, &    
                      grid(i)%eigv,grid(i)%eigvec)
        endif
     END IF
  END DO
  DEALLOCATE(grid)
1 FORMAT('***********************************************'   &
         '*************************'                         &
         /,10X,'Time to pack the spatial Hamiltonian = '     &
               ,f15.8,/                                      &
         '***********************************************'   &
         '*************************')
2    FORMAT(/,5X,'beginning arnoldi iterations:',/,5X,  &
                 'initial error = 'e15.8)
3    FORMAT(/,10X,'total time for lanczos iteration                 = '    &
                                                                    f15.8, &
            /,10X,'total time for eigenvalues of tridiagonal matrix = '    &
                                                                    f15.8, &
            /,10X,'total time for transformation of vectors         = '    &
                                                                    f15.8, &
            /,10X,'total time for exponentiation                    = '    &
                                                                    f15.8)
4    FORMAT(/,5X,'exiting arnoldi routine. maximal small matrix') 
5    FORMAT (/,5X,'cycle = ',i4) 
6    FORMAT(/,5X,'finished this time step')
7    FORMAT(/,5X,'iteration limit exceeded.',/,      &
              5X,'quit and return to main')
8    FORMAT('***********************************************'   &
            '*************************'                         &
            /,10X,'Time Summary for Interval         = ',i4,    &
            /,10X,'time for linear perturbation      = ',f15.8,   &
            /,10X,'time for right hand side          = ',f15.8,   &
            /,10X,'time for non-linear perturbation  = ',f15.8,   &
            /,10X,'time to propagate                 = ',f15.8,/, &
            /,10X,'Total                             = ',f15.8,/, &
            '***********************************************'     &
            '*************************')
9 FORMAT(/5x,'Survival Probablity',/,/15x,'Time', &
          15x,'Probability')
11 FORMAT(/10x,e15.8,5x,e15.8)
END SUBROUTINE lanczos_prop
