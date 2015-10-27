!
  program m275
  USE Time_Independent_Module
  USE Crank_Nicholson_Module
  USE SO_Module
  USE data
  IMPLICIT NONE
  Real(idp)               :: fpkey
  INTEGER                 :: intkey
  INTEGER                 :: lenth
  INTEGER                 :: len
  LOGICAL                 :: logkey
  LOGICAL                 :: dollar
  CHARACTER(LEN=80)       :: chrkey
!                                                                               
!                                                                               
!     initialize /io/
  CALL drum !    Open the input and output files
  CALL iosys ('read character options from rwf',-1,0,0,ops)  ! Read in the options string and the Variables
  CALL iosys ('read character "tdse filename" from rwf',-1,0,0,td_filename) ! Open the time dependent file
  CALL iosys ('open tdse as new',0,0,0,td_filename) ! Open the time dependent file
  Call pakstr(td_filename,len)
  Write(iout,5) td_filename(1:len)
  IF ( dollar('$time_dependent',card,cpass,inp )) THEN
       Step_Size = fpkey(card,'step_size',.01d0,' ')
       No_time_Steps = intkey(card,'number_of_time_steps',5,' ')
       IL = intkey(card,'lower_index',1,' ')
       IU = intkey(card,'upper_index',20,' ')
       Method = chrkey(card,'method_of_solution','cn_length',' ')     
       call pakstr(method,len)
       pulse = chrkey(card,'pulse','square',' ')     
       phase = fpkey(card,'phase',0.d0,' ')     
       omega = fpkey(card,'omega',.148d0,' ')     
       E_0   = fpkey(card,'electric_field',0.d0,' ')     
       E_con   = fpkey(card,'energy_convergence_parameter',1.d-08,' ')     
       delta_t   = fpkey(card,'time_step',.08d0,' ')     
       Prnt(1)=logkey(card,'print_eigenvectors',.false.,' ')
       Prnt(2)=logkey(card,'print_solution',.false.,' ')
       Prnt(3)=logkey(card,'print_propagators',.false.,' ')
       Prnt(4)=logkey(card,'print_scaled_vectors',.false.,' ')
       left_end = fpkey(card,'left_end',0.d0,' ')
       right_end = fpkey(card,'right_end',1.d0,' ')
       M_Size = (right_end - left_end)/step_size - 1 
       M_Eig = IU - IL + 1
       M_Eig = min(M_Size,M_Eig)
       Get_Eigenvectors = logkey(card,'get_eigenvectors',.false.,' ')
       IF ( Get_Eigenvectors == .true.) THEN
            Number_of_Eigenvectors=intkey(card,'number_of_eigenvectors',IU-IL+1,' ')
            Number_of_Eigenvectors=min(Number_of_Eigenvectors,M_Size)
            IU=min(IU,Number_of_Eigenvectors)
       END IF
  END IF
  IF ( method(1:len) == 'cn_length'                   .or.                     &
       method(1:len) == 'cn_velocity'                 .or.                     &
       method(1:len) == 'cn_length_and_cn_velocity'   .or.                     &
       method(1:len) == 'real_time_so_length'         .or.                     &
       method(1:len) == 'real_time_so_velocity'       .or.                     &
       method(1:len) == 'imaginary_time_so'           .or.                     &
       method(1:len) == 'real_time_so'                .or.                     &
       method(1:len) == 'diagonalize' )                        THEN
       write(iout,1) 
       Write(iout,2) method, left_end, right_end, m_size, get_eigenvectors
       write(iout,1) 
       IF (Method == 'cn_length' .or.  Method == 'cn_velocity' .or. Method == 'cn_length_and_cn_velocity') THEN
           Call H_0
           write(iout,1) 
           write(iout,3) No_time_Steps, delta_t, E_0, pulse, phase, omega
           write(iout,1) 
           write(iout,*) '          Print_Solution = ',prnt(2)
           Call CN_Driver
       ELSE IF (Method == 'real_time_so_length' .or. Method == 'real_time_so_velocity' .or.              &
                Method == 'imaginary_time_so') THEN
           CALL H_0
           write(iout,1) 
           write(iout,4) delta_t
           write(iout,*) '         Print_Solution       = ',prnt(2)
           write(iout,*) '         Print_Propagators    = ',prnt(3)
           write(iout,*) '         Print_Scaled_Vectors = ',prnt(4)
           CALL SO_Driver
      ELSE IF (Method == 'diagonalize') THEN
           Call H_0
      END IF    
  ELSE
      write(iout,*) '                        Error in input method'
      write(iout,*) 'Allowable methods are:  cn_length, cn_velocity, cn_length_and_cn_velocity'    
      write(iout,*) '                        real_time_so_length, real_time_so_velocity,'          
      write(iout,*) '                        imaginary_time_so and diagonalize'
      write(iout,*)  method(1:len)
      Call lnkerr('Error')
  END IF
  Write(iout,*) '               Calculation Finished'
  Call chainx(0)
1 Format(//,'********************************************************************************')
2 Format(/,30x,'Method = ',a32,1x,                                                                  &
         //,10x,'left boundary = ',f10.5,1x,'right boundary = 'f10.5,1x,                            &
         //,10x,'number of spatial points = ',i6,1x,'get_eigenvectors = ',l1)
3 Format(/,10x,'number of time steps = ',i8,1x, 'time step = ',f10.5,1x,                            &
         /,10x,'electric field = ',f10.5,1x,1x,'pulse = ',a16,1x,'phase = ',f10.5,1x,               &
                'omega = ',f10.5)
4 Format(/,10x,'time step = ',f10.5)
5 Format(/5x,'Opening Time Dependent file = ',a)
!********************************************************************************
  END PROGRAM m275
!********************************************************************************
