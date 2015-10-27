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
      PROGRAM mesa_data_transfer
      USE mesa_global
      IMPLICIT NONE
      CHARACTER(LEN=128)                    :: r_matrix_file_name
      LOGICAL                               :: dollar, logkey
      LOGICAL                               :: rearr
      INTEGER                               :: i, intkey
      CALL DRUM
      call iosys ('read character options from rwf',-1,0,0,ops)
      write (iout,10)
      call iosys('read integer "number of basis functions" from rwf', &
                  1,nbf,0,' ')
      drop=logkey(ops,'drop',.false.,' ') 
      IF(drop) THEN
         call iosys('read integer "old nbf" from rwf',1,oldnbf,0,' ')    
         idrop=1
         WRITE(iout,20)
      ELSE
         oldnbf=nbf
         idrop=0
         WRITE(iout,30)
      END IF
      IF ( dollar('$mesa_data_transfer',card,cpass,inp) ) then
          nmotot=intkey(card,'number-of-molecular-orbitals',1,' ')
          rearr=logkey(card,'rearrange-or-drop-mos',.false.,' ')
          IF(rearr) THEN
             call intarr(card,'orbitals-kept',list,nmotot,' ')
          END IF
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
!----------------------------------------------------------------------
!              put out transformation matrix                            
!                        and other                                    
!                   required matrices                                
!----------------------------------------------------------------------
      ALLOCATE(c(nbf,nbf),s(nbf,nbf),scr(nbf,nbf))
!----------------------------------------------------------------------
!             all matrices are written in ao and mo form               
!                          to rmtrx                                    
!----------------------------------------------------------------------
      call iosys('read character "transformation vector" from rwf',    &
                  -1,0,0,xform)
      call iosys('write character "transformation vector" to rmtrx',   &
                  0,0,0,xform)
      call iosys('read real '//xform//' from rwf',nbf*nbf,c,0,' ')
      if (rearr) then
          call newmat(c,scr,scr,'ao-mo')
      end if
      if (prnloc(2)) then
          title='transformation vectors'
          call prntrm(title,c,nbf,nmotot,nbf,nmotot,iout)
      endif
      call iosys ('write real '//xform//' to rmtrx',nbf*nmotot,    &
                   c,0,' ')
10    format (/,15x,'m8000: process and reformat mesa data')
20    format (/,20x,'***** dropping basis functions *****')
30    format (/,20x,'***** no basis functions dropped *****')
      call chainx(0)
      stop
END PROGRAM MESA_DATA_TRANSFER
