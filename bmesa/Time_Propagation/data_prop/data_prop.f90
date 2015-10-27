  PROGRAM data_input_propagator
  IMPLICIT NONE
  CHARACTER (LEN=8)                        :: parity
  CHARACTER (LEN=8), DIMENSION(10)         :: coord
  CHARACTER (LEN=8)                        :: itoc
  CHARACTER (LEN=80)                       :: chrkey, units, typke, tchar
  CHARACTER(LEN=80)                        :: input_file_name
  CHARACTER(LEN=80)                        :: pr_dvr
  CHARACTER(LEN=80)                        :: line, title
  CHARACTER(LEN=16)                        :: typwt
  CHARACTER(LEN=16)                        :: fptoc
  LOGICAL                                  :: logkey
  LOGICAL                                  :: genpts, nlse, plot
  LOGICAL                                  :: space, keep_diag, diag
  LOGICAL                                  :: absorb, proj
  LOGICAL                                  :: reuse, nodiag, automte
  LOGICAL, DIMENSION(2)                    :: fix, drop
  CHARACTER(LEN=8)                         :: prop_order
  INTEGER                                  :: plot_step, intkey, iostat
  INTEGER                                  :: spdim, i, j, k, ii
  INTEGER                                  :: length, ilen, jlen, nfix
  INTEGER                                  :: nreg, nblock, ntreg
  INTEGER                                  :: angmom, order, steps
  CHARACTER (LEN=8), DIMENSION(500)        :: n
  CHARACTER(LEN=8), DIMENSION(500)         :: edge
  CHARACTER(LEN=8)                         :: mass, t_init, deltat
  REAL*8, DIMENSION(10)                    :: omega_d
  REAL*8                                   :: natoms, mass_os
  REAL*8                                   :: scattering_length
  REAL*8                                   :: tmp
!
  WRITE(5,*) '          Terminal Input of Variables'
  WRITE(5,*) '          Variables with ? require a .t. or .f. answer'
  WRITE(5,*) '          Strings with embedded blanks and special characters'
  WRITE(5,*) '          must be enclosed in double quotes'
  WRITE(5,*) '          Input Filename'
  READ(5,*) input_file_name
  OPEN (UNIT=10,FILE=input_file_name,ACCESS='sequential',FORM='formatted',       &
        IOSTAT=IOSTAT,STATUS='new')
  tchar='$route'
  WRITE(10,1) tchar
  line=''
  WRITE(10,1) line
  tchar='$end'
  WRITE(10,2) tchar
  tchar='$nonstd'
  WRITE(10,3) tchar
  WRITE(5,*) 'Enter Route'
  tchar=''
  DO WHILE(tchar /= '$end')
     READ(5,*) tchar
     IF(tchar /= '$end') THEN
        WRITE(10,3) tchar
     END IF
  END DO
  tchar='$end'
  WRITE(10,2) tchar
  tchar='$title'
  write(10,1) tchar
  WRITE(5,*) 'Enter Title Card'
  READ(5,*) title
  WRITE(10,5) title
  tchar='$end'
  WRITE(10,2) tchar
  tchar='$dvrprop_basis'
  WRITE(10,4) tchar
  WRITE(5,*) '$dvrprop_basis keyword'
  WRITE(5,*) 'number-of-space-variables, propagation-order'
  READ(5,*) spdim, prop_order
  line='number-of-space-variables='//itoc(spdim)
  i=length(line) 
  line=line(1:i)//' propagation-order='//prop_order
  i=length(line)
  WRITE(10,*) line(1:i)
  WRITE(5,*) 'coordinate system'
  READ(5,*) line
  i=length(line)  
  line='coordinate-system='//line(1:i)
  i=length(line)  
  WRITE(10,*) line(1:i)
  DO i=1, spdim
     tchar=itoc(i)
     ilen=length(tchar)
     WRITE(5,*) 'label-for-space-variable-'//tchar(1:ilen)
     READ(5,*) line
     j=length(line)
     coord(i)=line(1:j)
     jlen=length(coord(i))
     line='space-variable-'//tchar(1:ilen)//'='//coord(i)(1:jlen)
     j=length(line)
     WRITE(10,*) line(1:j)
  END DO
  WRITE(5,*) 'units, kinetic-energy-type, automate-points ?,'
  WRITE(5,*) 'non-linear-equations ?, plot ?'     
  READ(5,*) units, typke, genpts, nlse, plot
  i=length(units)
  line='units='//units(1:i)
  i=length(line)
  j=length(typke)
  line= line(1:i)//' kinetic-energy-type='//typke(1:j)
  i=length(line)
  WRITE(10,*) line(1:i)
  line=''
  IF(genpts) THEN
     line='automate-points'
  END IF
  IF(nlse) THEN
     i=length(line)
     line = line(1:i)//' non-linear-equations'
  END IF
  IF(plot) THEN
     i=length(line)
     line = line(1:i)//' plot'
     WRITE(5,*) 'plot step'
     READ(5,*)  plot_step     
     tchar=itoc(plot_step)
     ilen=length(tchar)
     i=length(line)
     line=line(1:i)//' plot_step='//tchar(1:ilen)
  END IF
  i=length(line)
  IF(i /= 0) then
     WRITE(10,*) line(1:i)
  END IF
  line=''
  WRITE(5,*) 'no-spatial-hamiltonian ?, keep-diagonals ?, get-eigenpairs ?'
  READ(5,*) space, keep_diag, diag
  IF(space) THEN
     line='no-spatial-hamiltonian'
  END IF
  IF(keep_diag) THEN
     i=length(line)
     line=line(1:i)//' keep-diagonals'
  END IF
  IF(diag) THEN
     i=length(line)
     line=line(1:i)//' get-eigenpairs'
  END IF
  i=length(line)
  WRITE(10,*) line(1:i)
  line=''
  WRITE(5,*) 'print parameter'
  READ(5,*) pr_dvr
  line='print='//pr_dvr
  i=length(line)
  WRITE(10,*) line(1:i)
  line=''
  WRITE(5,*) 'add-absorbing-potential ? projections ?'
  READ(5,*) absorb, proj
  IF(absorb) THEN
     i=length(line)
     line=line(1:i)//' add-absorbing-potential'
  END IF
  IF(proj) THEN
     i=length(line)
     line=line(1:i)//' projections'
  END IF
  i=length(line)
  IF( i/=0) THEN
     WRITE(10,*) line(1:i)
  END IF
  tchar='$end'
  WRITE(10,2) tchar
  IF(nlse) THEN
     tchar='$nlse'
     WRITE(10,10) tchar
     WRITE(5,*) 'oscillator-frequencies-in-hertz, ' 
     READ(5,*) (omega_d(i), i=1,spdim)
     WRITE(5,*) 'oscillator-mass, number-of-atoms, '
     WRITE(5,*) 'scattering-length'
     READ(5,*)  mass_os, natoms, scattering_length
     WRITE(5,*) 'The frequencies will now be reordered from '
     WRITE(5,*) 'smallest to largest.  The entered grid needs '
     WRITE(5,*) 'to reflect this reordering'
!
!    Order them smallest to largest
!
     IF ( spdim > 1 ) THEN
          DO  ii=2,spdim
              i=ii-1
              k=i
              tmp=omega_d(i)
              DO  j=ii,spdim
                  IF(omega_d(j) < tmp) THEN
                     k=j
                     tmp=omega_d(j)
                  END IF
              END DO
              IF(k /= i) THEN
                 omega_d(k) = omega_d(i)
                 omega_d(i) = tmp
              END IF
          END DO
     END IF
     IF(spdim == 1) THEN
        WRITE(10,11) omega_d(1)
     ELSE IF(spdim == 2) THEN
        WRITE(10,12) (omega_d(i), i=1,2)
     ELSE IF(spdim == 3) THEN
        WRITE(10,13) (omega_d(i), i=1,3)
     END IF
     WRITE(10,14) natoms, mass_os
     WRITE(10,15) scattering_length
     tchar='$end'
     WRITE(10,4) tchar
  END IF
  DO i=1, spdim
     j=length(coord(i))
     tchar='$h0('//coord(i)(1:j)//')'
     WRITE(10,5) tchar
     WRITE(5,*) '$h0('//coord(i)(1:j)//') keyword'
     WRITE(5,*) 'main print variable control'
     READ(5,*) pr_dvr
     line='print='//pr_dvr
     WRITE(5,*) 'sector-print'
     READ(5,*) pr_dvr
     j=length(line)
     line=line(1:j)//' sector-print='//pr_dvr
     j=length(line)
     WRITE(10,*) line(1:j)
     IF(units /= 'atomic-units') THEN
        WRITE(5,*) 'mass'
        READ(5,*) mass
        line='mass='//mass
     END IF
     line=''
     WRITE(5,*) 'angular-momentum, parity'
     READ(5,*) angmom, parity
     j=length(line)
     tchar=itoc(angmom)
     ilen=length(tchar)    
     line=line(1:j)//' angular-momentum='//tchar(1:ilen)
     j=length(line)
     line=line(1:j)//' parity='//parity
     j=length(line)
     WRITE(10,*) line(1:j)
     line=''
     WRITE(5,*) 'reuse-space-data ?, do-not-diagonalize ?'
     WRITE(5,*) 'automate ?, weight-type'
     READ(5,*) reuse, nodiag, automte, typwt
     IF(reuse) THEN
        line='reuse-space-data'
     END IF
     IF(nodiag) THEN
        j=length(line)
        line=line(1:j)//' do-not-diagonalize'
     END IF
     IF(automte) THEN
        j=length(line)
        line=line(1:j)//' automate'
     END IF
     j=length(line)
     line=line(1:j)//' weight-type='//typwt
     j=length(line)
     IF( j /= 0 ) THEN
         WRITE(10,*) line(1:j)    
     END IF
     line=''
     WRITE(5,*) 'number-of-fixed-points'
     READ(5,*) nfix
     tchar=itoc(nfix)
     ilen=length(tchar)
     line='number-of-fixed-points='//tchar(1:ilen)
     j=length(line)
     write(10,*) line(1:j)
     line=''
     IF(nfix /= 0) THEN
        WRITE(5,*) 'left-fixed-point ?, right-fixed-point ?'
        READ(5,*) fix(1), fix(2)
        IF(fix(1)) THEN
           line='left-fixed-point'
           WRITE(5,*) 'drop-left-function ?'
           READ(5,*) drop(1)
           IF(drop(1)) THEN
              j=length(line)
              line=line(1:j)//' drop-left-function'
           END IF
        END IF
        IF(fix(2)) THEN
           j=length(line)
           line=line(1:j)//' right-fixed-point'
           WRITE(5,*) 'drop-right-function ?'
           READ(5,*) drop(2)
           IF(drop(2)) THEN
              j=length(line)
              line=line(1:j)//' drop-right-function'
           END IF
        END IF
        j=length(line)
        IF( j/=0 ) THEN
           WRITE(10,*) line(1:j)
        END IF
     END IF
     IF(typke == 'dvr'.OR.typke == 'packed') THEN
        IF(.not.automte) then
           WRITE(5,*) 'number-of-regions, edges'
           READ(5,*) nreg,(edge(j),j=1,nreg+1)
           tchar=itoc(nreg)
           ilen=length(tchar)
           line='number-of-regions='//tchar(1:ilen)
           j=length(line)
           line=line(1:j)//' region-boundaries=('//edge(1)
           DO j=2,nreg
              k=length(line)
              line=line(1:k)//','//edge(j)
           END DO
           k=length(line)
           line=line(1:k)//','//edge(nreg+1)
           k=length(line)
           line=line(1:k)//')'
           k=length(line)
           WRITE(10,*) line(1:k)
           WRITE(5,*) 'polynomial-order-per-region'
           READ(5,*) (n(j),j=1,nreg)
           ilen=length(n(1))
           line='polynomial-order-per-region=('//n(1)(1:ilen)
           k=length(line)
           line=line(1:k)
           DO j=2,nreg-1
              ilen=length(n(j))
              k=length(line)
              line=line(1:k)//','//n(j)(1:ilen)
           END DO
           k=length(line)
           ilen=length(n(nreg))
           line=line(1:k)//','//n(nreg)(1:ilen)//')'
           k=length(line)
           WRITE(10,*) line(1:k)
           WRITE(5,*) 'number-of-reference-quadrature-points-per-region'
           READ(5,*) (n(j),j=1,nreg)
           ilen=length(n(1))
           line='number-of-reference-quadrature-points-per-region='// &
                 '('//n(1)(1:ilen)
           k=length(line)
           line=line(1:k)
           DO j=2,nreg-1
              ilen=length(n(j))
              k=length(line)
              line=line(1:k)//','//n(j)(1:ilen)
           END DO
           k=length(line)
           ilen=length(n(nreg))
           line=line(1:k)//','//n(nreg)(1:ilen)//')'
           k=length(line)
           WRITE(10,*) line(1:k)
           tchar='$end'
           WRITE(10,2) tchar
           j=length(coord(i))
           line='$v_reg_1('//coord(i)(1:j)//')'
           WRITE(10,5) line
           WRITE(5,*) 'Enter regional potential information'
           line=''
           DO WHILE(line /= '$end')
              READ(5,*) line
              j=length(line)
              IF(line(1:j) /= '$end') THEN
                 WRITE(10,*) line(1:j)
              END IF
           END DO
           tchar='$end'
           WRITE(10,2) tchar
        ELSE
           nreg=0
           WRITE(5,*) 'Automated Selection of Steps'
           WRITE(5,*)'number-of-major-blocks'
           READ(5,*) nblock
           tchar=itoc(nblock)
           ilen=length(tchar)
           line='number-of-major-blocks='//tchar(1:ilen)
           k=length(line)
           WRITE(10,*) line(1:k)          
           tchar='$end'
           WRITE(10,2) tchar
           DO j=1,nblock
              tchar=itoc(j)
              ilen=length(tchar)
              tchar='$block-'//tchar(1:ilen)
              ilen=length(tchar)
              WRITE(5,*) tchar(1:ilen)//' keyword'
              WRITE(10,9) tchar
              WRITE(5,*) 'number-of-subregions, default-order,'
              WRITE(5,*) 'left-boundary, right-boundary'
              READ(5,*) n(1), n(2), edge(1), edge(2)
              ilen=length(n(1))
              line='number-of-subregions='//n(1)(1:ilen)
              k=length(line)
              ilen=length(n(2))
              line=line(1:k)//' default-order='//n(2)(1:ilen)
              k=length(line)
              WRITE(10,*) line(1:k)
              ilen=length(edge(1))
              line=' left-boundary='//edge(1)(1:ilen)
              k=length(line)
              ilen=length(edge(2))
              line=line(1:k)//' right-boundary='//edge(2)(1:ilen)
              k=length(line)
              WRITE(10,*) line(1:k)
              tchar='$end'
              WRITE(10,2) tchar
              tchar=itoc(j)
              ilen=length(tchar)
              jlen=length(coord(i))
              tchar='$v_reg_'//tchar(1:ilen)
              ilen=length(tchar)
              tchar=tchar(1:ilen)//'('//coord(i)(1:jlen)//')'
              WRITE(10,5) tchar
              WRITE(5,*) 'Enter regional potential information'
              line=''
              DO WHILE(line /= '$end')
                 READ(5,*) line
                 k=length(line)
                 IF(line(1:k) /= '$end') THEN
                    WRITE(10,*) line(1:k)
                 END IF
              END DO
              tchar='$end'
              WRITE(10,2) tchar
           END DO
        END IF
     ELSE
        WRITE(5,*) 'order-of-finite-difference-formula, number-of-steps'
        READ(5,*) order, steps
        tchar=itoc(order)
        ilen=length(tchar)
        line='order-of-finite-difference-formula='//tchar(1:ilen)
        j=length(line)
        tchar=itoc(steps)
        ilen=length(tchar)
        line=line(1:j)//' number-of-steps='//tchar(1:ilen)
        j=length(line)        
        WRITE(10,*) line(1:j)
        tchar='$end'
        WRITE(10,2) tchar
        line='$v_reg_1'//'('//coord(i)(1:jlen)//')'
        j=length(line)
        WRITE(10,5) line
        WRITE(5,*) 'Enter regional potential information'
        line=''
        DO WHILE(line /= '$end')
           READ(5,*) line
           j=length(line)
           IF(line(1:j) /= '$end') THEN
              WRITE(10,*) line(1:j)
           END IF
        END DO
        tchar='$end'
        WRITE(10,2) tchar
     END IF
  END DO
  WRITE(5,*) '$time keyword'
  tchar='$time'
  WRITE(10,1) tchar
  WRITE(5,*) 'initial-time, number-of-time-regions, time-interval'
  READ(5,*) t_init, ntreg, deltat
  line='initial-time='//t_init
  j=length(line)
  tchar=itoc(ntreg)
  ilen=length(tchar)
  line=line(1:j)//' number-of-time-regions='//tchar(1:ilen)
  j=length(line)
  line=line(1:j)//' time-interval='//deltat
  j=length(line)
  WRITE(10,*) line(1:j)
  WRITE(5,*) ' Main print control parameter'
  READ(5,*) pr_dvr
  line=' print=main='//pr_dvr
  j=length(line)
  WRITE(10,*) line(1:j)
  tchar='$end'
  WRITE(10,2) tchar
  tchar='$v0(t1)'
  WRITE(10,3) tchar
  WRITE(5,*) 'Enter pure time potential information'
  line=''
  DO WHILE(line /= '$end')
     READ(5,*) line
     j=length(line)
     IF(line(1:j) /= '$end') THEN
        WRITE(10,*) line(1:j)
     END IF
  END DO
  tchar='$end'
  WRITE(10,2) tchar
  tchar='$v_couple'
  WRITE(10,8) tchar
  WRITE(5,*) 'Enter Information on the nature of the coupling potential' 
  WRITE(5,*) 'and any parameters needed to describe its mathematical form' 
 line=''
  DO WHILE(line /= '$end')
     READ(5,*) line
     j=length(line)
     IF(line(1:j) /= '$end') THEN
        WRITE(10,*) line(1:j)
     END IF
  END DO 
  tchar='$end'
  WRITE(10,2) tchar
  tchar='$initial-state'
  WRITE(10,4) tchar
  WRITE(5,*) 'Enter Information on initial state' 
  line=''
  DO WHILE(line /= '$end')
     READ(5,*) line
     j=length(line)
     IF(line(1:j) /= '$end') THEN
        WRITE(10,*) line(1:j)
     END IF
  END DO 
  tchar='$end'
  WRITE(10,2) tchar

1  FORMAT(a6)
2  FORMAT(a4)
3  FORMAT(a7)
4  FORMAT(a14)
5  FORMAT(a24)
7  FORMAT(a80)
8  FORMAT(a9)
9  FORMAT(a10)
10 FORMAT(a5)
11 FORMAT('oscillator-frequencies-in-hertz=',d15.8)
12 FORMAT('oscillator-frequencies-in-hertz=('d15.8,',',  &
                                            d15.8,')')
13 FORMAT('oscillator-frequencies-in-hertz=('d15.8,',',  &
                                             d15.8,',',  &
                                             d15.8,')')
14 FORMAT('number-of-atoms=',d15.8,' mass=',d15.8)
15 FORMAT('scattering-length=',d15.8)
  stop
END PROGRAM data_input_propagator






