c \documentclass{article}
c \usepackage{graphicx}
c \setkeys{Gin}{width=\linewidth}
c \title{DVRDVD: A Davidson Code Specialized to the DVR Representation}
c \author{Barry I. Schneider}
c \date{}
c \def \<{\langle}
c \def \>{\rangle}
c \begin{document}
c \maketitle

*deck dvrdvd.f
c***begin prologue     dvrdvd
c***date written       010828   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c***keywords           iterative eigenvalue
c***author             schneider, barry (nsf)
c***source             
c***purpose            find the eigenvalues of a real symmetric matrix
c***                   using the davidson algorithm.  this version is
c***                   specialized to dvr basis sets.  
c***
c***description        the davidson algorithm as programmed here is
c***                   capable of getting the lowest m roots or of
c***                   climbing the eigenvalue tree in a piecewise fashion
c***                   by getting the first m1 roots, then the next m2 roots
c***                   then the next.... roots until all are obtained.
c***                   it can often happen that roots will be missed, 
c***                   especially if one tries to compute too many 
c***                   at each pass.  so, the routine must be used 
c***                   intelligently.  it is typically simple to see if a
c***                   root has been missed by inspection.  the method avoids
c***                   diagonalizing a large matrix for the higher roots by
c***                   explicitly orthogonalizing the guesses to the converged
c***                   eigenvectors of the lower roots.
c***
c***     variable      type       size             description
c***     _______       ____       ____             ___________
c***      hx            real       nx         one body hamiltonian in
c***                                          x coordinate with diagonal
c***                                          elements set to zero
c***      hy            real       ny         same in y coordinate
c***      hz            real       nz         same in z coordinate
c***      diag          real       n          diagonal part of matrix
c***      vec           real    n*maxvec      davidson vectors
c***      hvec          real    n*maxvec      action of the hamiltonian on pvec
c***      b, bwrk       real    maxvec*maxvec small hamiltonian matrix and a 
c***                                          copy
c***      eigwrk        real    maxvec        eigenvalues of small matrix
c***      work          real    max(5*maxvec, scratch array
c***                                n*maxvec)
c***      svec          real    maxvec*maxvec small scratch matrix
c***      resid         real    n*maxvec      residual vectors
c***      eig           real    nroot         converged eigenvalues
c***      cnverg        real                  convergence criterion for a root
c***      thresh        real                  overlap tolerance for accepting 
c***                                          a new davidson vector
c***      n             integer               matrix size
c***      nx            integer               dimension of x matrix
c***      ny            integer               same for y matrix
c***      nz            integer               same for z matrix
c***      nroot         integer               number of roots to get
c***      ntrial        integer               number of trial vectors
c***      nattim        integer               number of roots to converge at a 
c***                                          pass
c***      maxit         integer               maximum size of davidson space
c***      maxvec        integer               maximum number of vectors 
c***      nblck         integer               preconditioner maximum block size
c***                                          available for storage
c***      prnt          logical               print flags   
c***
c***      precon        character             type of preconditioner
c***      npass         integer               number of passes needed to get
c***                                          all roots 
c***      num2do        integer               number of roots to converge
c***                                          on this pass
c***      rtdone        integer               number of roots converged
c***      size          integer               current size of davidson space
c***references         
c
c***routines called    
c***end prologue       dvrdvd
      subroutine dvrdvd(hx,hy,hz,diag,eig,vec,hvec,resid,b,bwrk,
     1                  eigwrk,work,svec,cnverg,thresh,dim,
     2                  nx,ny,nz,n,nroot,ntrial,nattim,maxit,maxvec,
     3                  nblck,prnt,precon,pre)
      implicit integer (a-z)
      real*8 hx, hy, hz, diag, eig, vec, hvec, resid
      real*8 b, bwrk, eigwrk, work, svec
      real*8 error, cnverg, thresh
      real*8 maxerr
      integer*8 pre
      logical orth, prnt
      character*5 itoc
      character*8 cntrl, code
      character*80 title
      character*(*) precon
      dimension hx(nx,nx), hy(ny,ny), hz(nz,nz)
      dimension diag(n), eig(nroot)
      dimension vec(n,maxvec), hvec(n,maxvec), resid(n,maxvec)
      dimension b(maxvec,maxvec), bwrk(maxvec,maxvec)
      dimension eigwrk(maxvec), work(*), svec(maxvec,*)
      dimension code(3)
      dimension prnt(11)
      common/io/inp, iout
      data code / 'trial:', 'dvr-ene:', 'dvr-vec:' /
c
      write(iout,1) nroot, nattim, ntrial, maxit, maxvec, cnverg
c
c     find out how many passes are needed to get all the roots
c      
      npass=nroot/nattim
      nleft=nroot-npass*nattim
      if(nleft.ne.0) then
         npass=npass+1
      else
         nleft=nattim
      endif            
c
c     read in the trials and write then out in a more convenient form.
c
      call trilin(vec,code(1)(1:6),n,ntrial,maxvec,prnt(1))

c
      rtdone=0
      root0=0
      rootn=0
      num2do=nattim
      nin = nattim
      used = nin      
      call iosys('read real trials from rwf without rewinding',nin*n,
     1            vec,0,' ')
      do 1000 ipass=1,npass
         if(ipass.eq.npass) then
            num2do=nleft
         endif
         root0=rootn+1
         rootn=rootn+num2do    
         write(iout,2) ipass, root0, rootn      
c-----------------------------------------------------------------------c
c                                                                       c
c                    Initialization Section                             c
c                                                                       c
c-----------------------------------------------------------------------c
         orth=.true.
         call init(hx,hy,hz,diag,vec,hvec,b,bwrk,resid,thresh,
     1             'initialize',dim,n,nin,1,size,nx,ny,nz,
     2             rtdone,code(3),maxvec,orth,prnt)
c----------------------------------------------------------------------c
c                                                                      c
c                    Iteration Sequence                                c
c                                                                      c
c     iteration is continued until all of the roots are converged      c 
c     or if convergence is not achieved some preset maximum number of  c
c     iterations are performed.                                        c
c                                                                      c
c----------------------------------------------------------------------c
         iter=0
         error=1.d+10
         write(iout,3) error
         cntrl='continue'
         do while ( cntrl.eq.'continue'.and.iter.lt.maxit )
            iter = iter + 1
c        Step 1:
c
c           get eigenvalues and eigenvectors of the small matrix.
c
c                bwrk holds the initial matrix which is destroyed.
c                svec has the transformation matrix.
c                note that resid is used as temporary storage in rdiag.
c
c            title='iteration = '//itoc(iter)//' diagonalizing '//
c     1            'matrix of dimension N = '//itoc(size)
c
            if(prnt(4)) then
               title='initial small matrix'
               call prntfm(title,bwrk,size,size,maxvec,maxvec,iout)
            endif
c        
            call rdiag(bwrk,eigwrk,svec,work,0.d0,iter,maxvec,
     1                 size,prnt(5))
c      
c        Step 2:
c
c           transform vec and hvec to the new basis defined by the
c           diagonalization in rdiag.  the small matrix becomes diagonal under
c           this transformation and we fill it with the eigenvalues.
c           then form the residuals and check for convergence.  the converged
c           eigenvalues are stored and the eigenpairs written to disk.  the
c           unconverged residuals are moved from their current positions in 
c           resid to the leading positions.
c
            write(iout,4) iter, size
            call frmres(vec,hvec,resid,svec,eig,eigwrk,
     1                  b,bwrk,0.d0,cnverg,maxerr,n,size,
     2                  num2do,con,uncon,maxvec,iter,
     3                  prnt(6),rtdone,code(2))
c
c     
c           check to see if all num2do roots are converged
c           if so, we are done and can print results and then proceed to 
c           the next set of roots.
c 
            if(con.eq.num2do) then
c
               cntrl='finished'
               write(iout,5) num2do
               rtdone = rtdone + num2do
               if(rtdone.eq.nroot) then
c
c                 we not only done with these roots but are finished
c                 with all the roots.
c
                  write(iout,6)
                  return
               else
c            
c                 we need to prepare for the next set of roots.
c
c                 begin next set of roots with the best available
c                 vectors.  these are the the remaining vectors coming 
c                 from the diagonalization of the small matrix.
c
                  left=size-num2do
                  left=min(left,maxvec)
                  if(left.ne.0) then
                     write(iout,7) left
                     call copy(vec(1,num2do+1),vec,n*left)  
                  elseif(used.lt.ntrial) then
                     left = ntrial - used
                     n2rd=min(num2do,left)
                     write(iout,8) n2rd
                     used = used + n2rd
                     call iosys('read real trials from rwf '//
     1                          'without rewinding',n*n2rd,
     2                           vec,0,' ')
                  else
                     write(iout,9)
                     call lnkerr('quit. no available trials')
                  endif
               endif
            else
c
c              all roots are not converged.  set the error to the largest
c              current error and either restart or continue the 
c              iteration sequence.
c
               error=min(error,maxerr)
c
c              how many new vectors could be added in principle
c
               numnew = maxvec - size
c                                             
c              how many will we add
c
               addvec = min(numnew,uncon)
c              
c              check if the number of old plus new vectors will exceed
c              the maxvec limit to determine if a restart is needed.
c
               chkno = size + addvec    
               if(chkno.gt.maxvec) then
c             
c              the maximum number of vectors is exceeded.
c
                  min2kp = con + uncon + uncon
                  size = min(min2kp,size)
                  write(iout,11) size
                  numnew = maxvec - size
                  addvec = min(numnew,uncon)
               endif
               write(iout,12) addvec
c
c              maximum number of vectors is still within the allowed
c              limits.  add vectors to the set from the unconverged
c              residuals and put them after the current vectors.
c
               if(addvec.ne.0) then
                  beg = size + 1
                  call genvec(vec(1,beg),resid,diag,eigwrk,work,
     1                        nx,ny,nz,n,dim,pre,addvec,
     2                        iter,precon,nblck,prnt(10))
               else
c
c                 the only way to continue is if there are more
c                 available trial vectors
c
                  if(used.gt.ntrial) then
                     call lnkerr('cannot contine. no more vectors')
                  else
                     left = ntrial - used
                     n2rd=min(nattim,left)
                     write(iout,8) n2rd
                     used = used + n2rd
                     call iosys('read real trials from rwf '//
     1                          'without rewinding',n*n2rd,
     2                           vec(1,beg),0,' ')
                  endif            
               endif
               orth=.false.
               call init(hx,hy,hz,diag,vec,hvec,b,bwrk,
     1                   resid,thresh,'fill',dim,n,
     2                   addvec,beg,size,nx,ny,nz,rtdone,
     3                   code(3)(1:8),maxvec,orth,prnt)
               if(prnt(11)) then
                  call tstovl(vec,n,size)
               endif
            endif
            nin = size
         enddo
         if(iter.ge.maxit) then
            write(iout,13)   
            return
         endif        
 1000 continue
      write(iout,14)
      do 20 i=1,nroot
         write(iout,15) i, eig(i)
 20   continue           
      return
 1    format(/,1x,'davidson eigenvalue solver',/,10x,
     1            'number of roots               = ',i4,/,10x,
     2            'number of roots at a time     = ',i4,/,10x,
     3            'number of trials              = ',i4,/,10x,
     4            'maximum number of iterations  = ',i4,/,10x,
     5            'maximum number of vectors     = ',i4,/,10x,
     6            'convergence criterion for rms = ',e15.8)      
 2    format(/,5x,'pass = ',i3,/,5x,
     1            'processing root = ',i4,' to root =',i4)
 3    format(/,5x,'beginning davidson iterations:',/,5x,
     1            'initial error = 'e15.8) 
 4    format (/,5x,'cycle = ',i4,2x,'no. vectors = ',i4)
 5    format(/,5x,'all',1x,i5,1x,'roots are converged')
 6    format(/,5x,'we are done')
 7    format(/,5x,'preparing for the next set of roots',
     1       /,5x,'number of added vectors from unconverged roots = ',
     2                                                             i5)
 8    format(/,5x,'there are no vectors left from the diagonalization',
     1       /,5x,'read ',i4,' trial vectors for next guesses')
 9    format(/,5x,'no new trial vectors available for next pass')
 11   format(/,5x,'maximum number of vectors exceeded',
     1        /5x,'contract back to ',1x,i4,' vectors')
 12   format(/,5x,'number of vectors added = ',i4)                   
 13   format(/,5x,'iteration limit exceeded.',/,5x,
     1            'quit and return to main')
 14   format(///,30x,'final summary',/,17x,'root',22x,'energy')
 15   format(15x,i4,20x,f15.8)

      end       














