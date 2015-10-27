!deck dvr_dvd
!***begin prologue     dvr_dvd
!***date written       010828   (yymmdd)
!***revision date      yymmdd   (yymmdd)
!***keywords           iterative eigenvalue
!***author             schneider, barry (nsf)
!***source             
!***purpose            find the eigenvalues of a real symmetric matrix
!
!***                   using the davidson algorithm.  this version is
!                      specialized to dvr basis sets.  
!***
!***description        the davidson algorithm as programmed here is
!                      capable of getting the lowest m roots or of
!                      climbing the eigenvalue tree in a piecewise 
!                      fashion by getting the first m1 roots, 
!                      then the next m2 roots
!                      then the next.... roots until all are obtained.
!                      it can often happen that roots will be missed, 
!                      especially if one tries to compute too many 
!                      at each pass.  so, the routine must be used 
!                      intelligently.  it is typically simple to 
!                      see if a root has been missed by inspection.  
!                      the method avoids diagonalizing a large 
!                      matrix for the higher roots by
!                      explicitly orthogonalizing the guesses to 
!                      the converged eigenvectors of the lower roots.
!
!***references         
!
!***routines called    
!***end prologue       dvr_dvd
  SUBROUTINE dvr_dvd
  USE io
  USE dvd_global
  USE dvd_prnt
  IMPLICIT NONE
  INTEGER                 :: i
  character*5 itoc
  character*8 cntrl
!
  write(iout,1) nroot, nattim, ntrial, maxit, maxvec, cnverg
!
!         find out how many passes are needed to get all the roots
!      
  npass=nroot/nattim
  nleft=nroot-npass*nattim
  if(nleft /= 0) then
     npass=npass+1
  else
     nleft=nattim
  END IF            
!
!     read in the trials and write then out in a more convenient form.
!
  call trilin

!
  rtdone=0
  root0=0
  rootn=0
  num2do=nattim
  nin = nattim
  used = nin      
  call iosys('read real trials from rwf without rewinding',  &
              nin*n_dvd,vec,0,' ')
  DO ipass=1,npass
     IF(ipass == npass) then
        num2do=nleft
     END IF
     root0=rootn+1
     rootn=rootn+num2do    
     write(iout,2) ipass, root0, rootn      
!--------------------------------------------------------------------
!                 Initialization Section                             
!--------------------------------------------------------------------
     orth=.true.
     begin = 1
     drctv = 'initialize'
     call init(code(1))

!--------------------------------------------------------------------
!                 Iteration Sequence                               
!     iteration is continued until all of the roots are converged    
!     or if convergence is not achieved some preset maximum number of  
!--------------------------------------------------------------------
     iter=0
     error=1.d+10
     write(iout,3) error
     cntrl='continue'
     DO WHILE ( cntrl == 'continue'.and.iter < maxit )
        iter = iter + 1
!        Step 1:
!           get eigenvalues and eigenvectors of the small matrix.
!
!                bwrk holds the initial matrix which is destroyed.
!                svec has the transformation matrix.
!                note that resid is used as temporary storage in rdiag.
!
!            title='iteration = '//itoc(iter)//' diagonalizing '// &
!                  'matrix of dimension N = '//itoc(size)
!
        if(log_dvd(4)) then
           title='initial small matrix'
           call prntfm(title,bwrk,size,size,maxvec,maxvec,iout)
        END IF
!        
        rep = 0.d0
        call rdiag
!        Step 2:
!           transform vec and hvec to the new basis defined by the
!           diagonalization in rdiag.  the small matrix becomes diagonal under
!           this transformation and we fill it with the eigenvalues.
!           then form the residuals and check for convergence.  the converged
!           eigenvalues are stored and the eigenpairs written to disk.  the
!           unconverged residuals are moved from their current positions in 
!           resid to the leading positions.
        write(iout,4) iter, size
        call frmres(code(2))
!           check to see if all num2do roots are converged
!           if so, we are done and can print results and then proceed to 
!           the next set of roots.
        if(con == num2do) then
           cntrl='finished'
           write(iout,5) num2do
           rtdone = rtdone + num2do
           if(rtdone == nroot) then
!                 we are not only DOne with these roots but are finished
!                 with all the roots.
              write(iout,6)
              return
           else
!                 we need to prepare for the next set of roots.
!                 begin next set of roots with the best available
!                 vectors.  these are the the remaining vectors coming 
!                 from the diagonalization of the small matrix.
!
              left=size-num2do
              left=min(left,maxvec)
              if(left /= 0) then
                 write(iout,7) left
                 vec(:,1:left) = vec(:,num2do+1:left+num2do)
              elseif(used < ntrial) then
                 left = ntrial - used
                 n2rd=min(num2do,left)
                 write(iout,8) n2rd
                 used = used + n2rd
                 call iosys('read real trials from rwf '//               &
                            'without rewinding',n_dvd*n2rd,              &
                             vec,0,' ')
              else
                 write(iout,9)
                 call lnkerr('quit. no available trials')
              END IF
           END IF
        else
!              all roots are not converged.  set the error to the largest
!              current error and either restart or continue the 
!              iteration sequence.
!
           error=min(error,maxerr)
!              how many new vectors could be added in principle
           numnew = maxvec - size
!              how many will we add
           addvec = min(numnew,uncon)
!              check if the number of old plus new vectors will exceed
!              the maxvec limit to determine if a restart is needed.
           chkno = size + addvec    
           if(chkno > maxvec) then
!              the maximum number of vectors is exceeded.
              min2kp = con + uncon + uncon
              size = min(min2kp,size)
              write(iout,11) size
              numnew = maxvec - size
              addvec = min(numnew,uncon)
           END IF
           write(iout,12) addvec
!              maximum number of vectors is still within the allowed
!              limits.  add vectors to the set from the unconverged
!              residuals and put them after the current vectors.
           if(addvec /= 0) then
              begin = size + 1
              call genvec(vec(1,begin))
           else
!                 the only way to continue is if there are more
!                 available trial vectors
              if(used > ntrial) then
                 call lnkerr('cannot contine. no more vectors')
              else
                 left = ntrial - used
                 n2rd=min(nattim,left)
                 write(iout,8) n2rd
                 used = used + n2rd
                 call iosys('read real trials from rwf '//               &
                            'without rewinding',n_dvd*n2rd,              &
                             vec(1,begin),0,' ')
              END IF            
           END IF
           orth=.false.
           drctv = 'fill'
           call init(code(3))
           if(log_dvd(11)) then
              call tstovl(vec,n_dvd,size)
           END IF
        END IF
        nin = size
     END DO
     if(iter >= maxit) then
        write(iout,13)   
        return
     END IF        
  END DO
  write(iout,14)
  DO i=1,nroot
     write(iout,15) i, eig(i)
  END DO           
1 format(/,1x,'davidson eigenvalue solver',/,10x,                        & 
              'number of roots               = ',i4,/,10x,               &
              'number of roots at a time     = ',i4,/,10x,               &
              'number of trials              = ',i4,/,10x,               &
              'maximum number of iterations  = ',i4,/,10x,               &
              'maximum number of vectors     = ',i4,/,10x,               &
              'convergence criterion for rms = ',e15.8)      
2 format(/,5x,'pass = ',i3,/,5x,                                         &
              'processing root = ',i4,' to root =',i4)
3 format(/,5x,'beginning davidson iterations:',/,5x,                     &
              'initial error = 'e15.8) 
4 format (/,5x,'cycle = ',i4,2x,'no. vectors = ',i4)
5 format(/,5x,'all',1x,i5,1x,'roots are converged')
6 format(/,5x,'we are DOne')
7 format(/,5x,'preparing for the next set of roots',                     &
         /,5x,'number of added vectors from unconverged '                &
              'roots = ',i5)
8 format(/,5x,'there are no vectors left from the '                      &
              'diagonalization',                                         &
         /,5x,'read ',i4,' trial vectors for next guesses')
9 format(/,5x,'no new trial vectors available for next pass')
11 format(/,5x,'maximum number of vectors exceeded',                     &
          /5x,'contract back to ',1x,i4,' vectors')
12 format(/,5x,'number of vectors added = ',i4)                   
13 format(/,5x,'iteration limit exceeded.',/,5x,                         &
               'quit and return to main')
14 format(///,30x,'final summary',/,17x,'root',22x,'energy')
15 format(15x,i4,20x,f15.8)
END SUBROUTINE dvr_dvd       














