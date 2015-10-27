!
MODULE input_output
!***begin prologue     input_output
!***date written       021231   (yymmdd)
!***revision date               (yymmdd)
!***keywords           
!***author             schneider, b. i.(nsf)
!***source             
!***purpose            To set up input and output files and to allow file substitution
!                      on the execute line.  Very convenient
!***references

!***routines called    
!***end prologue       input_output

  IMPLICIT NONE
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!                   Some Constants and Units
!
   INTEGER                                         :: inp=8
   INTEGER                                         :: iout=9
   INTEGER                                         :: input
   INTEGER                                         :: output
   INTEGER                                         :: rows_to_print
   INTEGER                                         :: columns_to_print
   INTEGER                                         :: eigenvectors_to_print
   LOGICAL                                         :: print_paramete
   INTEGER                                         :: nunits
   INTEGER                                         :: n_dir
   CHARACTER (LEN=80), DIMENSION(:), ALLOCATABLE   :: names
   INTEGER, DIMENSION(:), ALLOCATABLE              :: len_name
   CHARACTER (LEN=80), DIMENSION(:), ALLOCATABLE   :: File_Directory
   INTEGER, DIMENSION(:), ALLOCATABLE              :: len_dir
   CHARACTER(LEN=8),  DIMENSION(:), ALLOCATABLE    :: rowlab, collab
   common / io / input, output   
!
!
!***********************************************************************
!***********************************************************************
                              CONTAINS
!***********************************************************************
!***********************************************************************
!***begin prologue     Command_Line
!***date written       021231   (yymmdd)
!***revision date               (yymmdd)
!***keywords           
!***author             schneider, b. i.(nsf)
!***source             Modules
!***purpose            
!***description        
!***                   
!
!***references

!***routines called    
!***end prologue       Command_Line
!  Subroutine Command_Line(home_directory,nunits,n_direct)
  Subroutine Command_Line(home_directory)
  IMPLICIT NONE
  INTEGER                                        :: i
  INTEGER                                        :: st
  INTEGER                                        :: input
  INTEGER                                        :: output
  CHARACTER(*), OPTIONAL                         :: home_directory
  common / io / input, output   
  input = inp
  output = iout
  names(1) = 'home'
  names(2) = 'inp'
  names(3) = 'out'
  call comand (nunits,names)
  IF( PRESENT(home_directory) ) THEN  
       File_Directory(1) = home_directory      ! Default home directory
       names(1) =home_directory
  ELSE
       Call lnkerr('no home directory specified')
  END IF
  DO i = 1,nunits
     call pakstr(names(i),len_name(i))
  END DO
  IF (names(2) == 'inp') THEN
      names(2)='input'
  END IF
  IF (names(3) == 'out') THEN
      names(3)='output'
  END IF
  call pakstr(File_Directory(1),len_dir(1))
  DO i=1,nunits
     Call pakstr(names(i),len_name(i))
  END DO
!
!                The other directories
!
  File_Directory(2) = File_Directory(1)(1:len_dir(1))//'/Inp'
  File_Directory(3) = File_Directory(1)(1:len_dir(1))//'/Out'
  File_Directory(4) = File_Directory(1)(1:len_dir(1))//'/Data_Files'
  File_Directory(5) = File_Directory(1)(1:len_dir(1))//'/Run'
  File_Directory(6) = File_Directory(1)(1:len_dir(1))//'/Library'
  DO i=1,6
     Call pakstr(File_Directory(i),len_dir(i))
  END DO
!
!
!     initialize /io/
!
  open (input,file=File_Directory(2)(1:len_dir(2))//'/'//names(2)(1:len_name(2)),access='sequential', &
                   form='formatted',iostat=st,status='old')
  IF (st > 0 ) THEN
      call abort('could not open the input file')
  END IF
  open (output,file=File_Directory(3)(1:len_dir(3))//'/'//names(3)(1:len_name(3)),access='sequential', &    
                    form='formatted',iostat=st,status='unknown')
  IF (st > 0 ) THEN
      call abort('could not open the output file')
  END IF
!
!
  Call versn(iout)  !  Some nice printout
  Write(iout,1) 
  Write(iout,2) input, output
  Write(iout,3) File_Directory(1:6)
  call wind(output)
1 Format(/,10x,'********************************************************************************'  &
         /,35x,'Setting up Directories and Files'                                                  &
         /,10x,'********************************************************************************'  )
2 Format(/10x,'input file number = ',i1,10x,'output file number = 'i1)
3 Format(/,10x,'home directory             = ',a64,                                                &
         /,10x,'input file directory       = ',a64,                                                & 
         /,10x,'output file directory      = ',a64                                                 &
         /,10x,'matrix directory           = ',a64                                                 &
         /,10x,'run files directory        = ',a64,                                                & 
         /,10x,'library files directory    = ',a64)
END SUBROUTINE Command_Line
!***********************************************************************
!***********************************************************************
!deck @(#)versn.f	1.2  5/30/91
      subroutine versn(iout)
!
!***begin prologue     versn
!***date written       850601  (yymmdd)
!***revision date      870207  (yymmdd)
!
!***keywords 
!***author             martin, richard (lanl)
!***source             @(#)versn.f	1.2   5/30/91
!***purpose            prints the date and current version of the mesa system.
!***description
!                      call versn(iout)
!                         iout   output file.
!
!***references
!***routines called    (none)
!***end prologue       versn
!
  IMPLICIT NONE
!
  INTEGER, PARAMETER            :: nmx=2
  INTEGER                       :: iout
  INTEGER                       :: cskipb      
  INTEGER                       :: pos
  INTEGER                       :: lenth
  INTEGER                       :: l_1
  INTEGER                       :: l_2
  CHARACTER(LEN=24)             :: today
  CHARACTER(LEN=120)            :: copyr,authrz,line
  CHARACTER(LEN=120)            :: machin, mch(nmx)
  CHARACTER(LEN=120)            :: mchtyp(nmx), mx
  CHARACTER(LEN=120)            :: site='National Science Foundation;'
  CHARACTER(LEN=16)             :: user
  INTEGER                       :: i
!
!
 1000 format(35x,a80)
 1020 format(25x,'User: 'a80,/)
 1010 format(25x,'Date: '80a1)
 1030 format(10x,a120)
 authrz ='B. I. Schneider.'
 mch(1) = 'rocks-fe'
 mch(2) = 'fermi'
 mchtyp(1) = 'Quad core APPRO infinband cluster(80 cores): AMD processor'
 mchtyp(2) = 'Quad core MICROWAY workstation(8 cores): Intel Westmere processors'
!
!     write the version information.
!
 write(iout,1020) authrz
!
!     accumulate and write operating conditions.
!
 call dattim(today)
 write(iout,1010) (today(i:i),i=1,lenth(today))
 line=site
!     if (hostnm(mx).ne.0) mx='unknown'
!
 mx='fermi'
 DO i=1,nmx
    IF (mx == mch(i)) THEN
        machin=mchtyp(i)
    END IF
 END DO 
 Call pakstr(mx,l_1)
 Call pakstr(machin,pos)
 call getlog(user)
 Call pakstr(user,l_2)
 line=machin(1:pos)//'('//mx(1:l_1)//'):'//user(1:l_2)
 write(iout,1030) line
!
!
  END Subroutine versn
!***********************************************************************
!***********************************************************************
  END MODULE input_output
!***********************************************************************
!***********************************************************************
