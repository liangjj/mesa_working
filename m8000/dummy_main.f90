!*deck mesa_data_transfer
!***begin prologue     m8000
!***date written       041206
!***revision date      (yymmdd)
!                      take molecular data and reformat for other codes.                   
!***keywords           
!***author             schneider, barry (nsf)
!***source             m8000
!***purpose            massage molecular data
!***
!
!***references
!
!***routines called    iosys, util and mdutil
!***end prologue       m8000
      PROGRAM dummy_main
      USE mesa_global
      IMPLICIT NONE
      CHARACTER(LEN=8)                      :: presnt, chrkey, drctv
      CHARACTER(LEN=16)                     :: fptoc
      CHARACTER(LEN=24)                     :: ftitl
      CHARACTER(LEN=32)                     :: xform
      CHARACTER(LEN=4)                      :: itoc
      CHARACTER(LEN=128)                    :: r_matrix_file_name
      LOGICAL                               :: dollar, logkey
      LOGICAL                               :: rearr, rearrc
      REAL*8                                :: fpkey, rmin, rmax, ksqmax
      REAL*8                                :: ene, reschg
      REAL*8, DIMENSION(100)                :: echan, energy
      INTEGER                               :: i, intkey
      CALL DRUM
      call iosys ('read character options from rwf',-1,0,0,ops)
      write (iout,500)
      call iosys('read integer "number of basis functions" from rwf', &
                  1,nbf,0,' ')
      drop=logkey(ops,'drop',.false.,' ') 
      IF(drop) THEN
         call iosys('read integer "old nbf" from rwf',1,oldnbf,0,' ')         
         idrop=1
      ELSE
         oldnbf=nbf
         idrop=0
      END IF
      IF ( dollar('$mesa_data_transfer',card,cpass,inp) ) then
          call intarr(card,'orbitals-kept',list,nmotot,' ')
          pntbuf=intkey(card,'point-buffer',100000,' ')
          DO i=1,15           
             prnloc(i)=logkey(card,prnkey(i),.false.,' ')
          END DO
      END IF
      CALL iosys('read character "r-matrix filename" from rwf',       &
                  -1,0,0,r_matrix_file_name)
      CALL iosys('open rmtrx as unknown',0,0,0,r_matrix_file_name)
!----------------------------------------------------------------------
!                 get the basis function information                        
!----------------------------------------------------------------------
      ALLOCATE(bflabl(maxnbf))
      call iosys('read character "basis function labels" from rwf',   &
                 -1,0,0,bflabl)
      call iosys('read integer "number of atoms" from rwf',1,         &
                  nat,0,' ')
      call iosys('length of exponents on rwf',nprim,0,0,' ')
      call iosys('length of "contraction coefficients" on rwf',       &
                  ncont,0,0,' ')
      call iosys('length of "number of pure functions" on rwf',       &
                  ntypes,0,0,' ')
      call iosys('read integer "number basis types" from rwf',        &
                  1,nbtype,0,' ')
      call iosys('length of "power of x" on rwf',lenxyz,0,0,' ')
!----------------------------------------------------------------------
!              divide core for basis set information                   
!----------------------------------------------------------------------
      nnp=nbf*(nbf+1)/2
      ALLOCATE(coords(3,nat),ex(nprim),cont(ncont),ptprim(nat,ntypes),  &
               noprim(nat,ntypes),nocont(nat,ntypes),                   &
               ptcont(nat,ntypes),start(nat,ntypes),zan(nat),           &
               nocart(ntypes),nobf(ntypes),maxmom(ntypes),              &
               mintyp(ntypes),nx(lenxyz),ny(lenxyz),nz(lenxyz),         &
               minmom(ntypes),index(oldnbf)) 
      CALL basout      
      call chainx(0)
      stop
  500 format (/,15x,'m8000: process and reformat mesa data')
END PROGRAM dummy_main
