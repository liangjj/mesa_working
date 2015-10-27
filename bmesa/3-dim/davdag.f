*deck @(#)davdag.f	1.5  7/30/91
      subroutine davdag(oper,status,file,d,v,hv,resid,n,maxvec,i3,
     1                  i4,i5,ipvt,small,root,eigvec,cvn,schmdt,prres)
c
c***begin prologue     davdag
c***date written       960302   (yymmdd)
c***revision date               (yymmdd)
c
c
c***keywords           eigenvalue, large matrices, davidson algorithm,
c                      liu algorithm, diagonalization
c***author             paul saxe (lanl), barry schneider(nsf)
c***purpose            iterative solution of an eigenvalue problem for
c                      one or a few of the lowest eigenvalues and
c                      vectors.
c***description
c
c    davdag is intended for extracting one or a few of the lowest
c    eigenvalues and vectors of a large symmetric matrix.  the major
c    computational step is typically the product of the matrix times
c    a vector. however, there are other situations where the calculation
c    can be dominated by the schmidt orthonormalization of the basis or
c    the repetetive diagonalization of the "small" matrix.  in addition,
c    there are some new procedures to solve the residual equation which
c    are also potentially time consuming.  the work required is
c    roughly proportional to the number of eigenvalues desired times
c    the dimension of the matrix squared.  this may be misleading for very
c    sparse matrices where the number of non-zero matrix elements in each row
c    in much less than the matrix dimension.
c
c    the user calls the subroutine to perform one specific task of the entire
c    davidson procedure.  that task is controlled by the first 
c    parameter, oper, in the calling sequence.  the sequence of the calls is
c    dictated by the procedure and is not arbitrary.  the description of
c    the tasks may be found in the body of the code.
c
c----------------------------------------------------------------------
c
c***references         b. liu, "the simultaneous expansion method for
c          the iterative solution of several of the lowest eigenvalues
c          and corresponding eigenvectors of large real-symmetric
c          matrices" in numerical algorithms in chemistry:
c          algebraic methods, a report of an nrcc workshop, 1978.
c
c          e. r. davidson, j. computational phys. 17, p87 (1975).
c
c***end prologue       davdag
c
c
      implicit integer (a-z)
c
      parameter (mxroot=500)
c
      real*8 d(n), v(n,*), hv(n,*), resid(n,*)
      real*8 small(maxvec,*), root(maxvec), eigvec(maxvec,*), cvn(*)
      real*8 sdot, maxdif, eroot, nrzero, norm, thresh, efixed, test
      real*8 ovrlap, cnvgnc
      dimension ipvt(n)
      character*(*) oper, status, file
      character*4 itoc
      character*3 truth
      character*80 title
      logical prnt, diag, schmdt, prres
      logical debug
      integer ptroot(mxroot)
c
      common /io/ inp,itape6
c
      data debug/.false./
      data nrzero / 1.0d-06 /
c
      save oldnum, added, totnum, nattim, numcnv
      save diag, cycle, maxdif, cnvgnc, efixed, nroots
      save ntodo, newtr, thresh
      save nrzero, ptroot, prnt, iout, lwr, upr
c
c
c**********************************************************************c
c         decide what to do based on "oper" passed in
c**********************************************************************c
c
c                      'initialize'
c
      if (oper(1:10).eq.'initialize') then
c**********************************************************************c
c              save local copies of important variables
c**********************************************************************c
c                     # of roots
         nroots=i3
c                     output file unit
         iout=i4
c                     # roots sought in each pass
         nattim=i5
c                     threshold for acceptance of a vector
         thresh=v(1,1)
c                     additive constant to total eigenvalue
         efixed=small(1,1)
c                     convergence criterion for root   
         cnvgnc=root(1)
         prnt=status.ne.'noprint'
         diag=oper(12:).eq.'with diagonals'
c*********************************************************************c
c        initialize some local variables and create some files
c*********************************************************************c
         call izero(ptroot,mxroot)
c        initialize some counters
         oldnum=0
         newtr=0
         added=0
         cycle=0
         totnum=0
         numcnv=0 
         lwr=1
         upr=nattim
c         
c*********************************************************************c
c           check for exceeding maximum number of expansion vectors
c*********************************************************************c
c                       'check matrix size'
c
      else if (oper.eq.'check matrix size') then
c
         num2ad=i3
         if (oldnum+added+num2ad.gt.maxvec) then
             status='exceeded'
         else
c*********************************************************************c
c                      update counters
c*********************************************************************c
             status='ok'
             added=added+num2ad
             totnum=totnum+num2ad
         end if
c
c*********************************************************************c
c            set up and solve the small set of equations
c            note that the number of expansion vectors
c            is equal to or greater than the number of roots
c*********************************************************************c
c         
c                            'solve'
c
      else if (oper.eq.'solve') then
         status='continue'
         if (oldnum.eq.0) then
             call rzero(small,maxvec*maxvec)
         end if
c
c*********************************************************************c
c        first, add product of trial and diagonals to product
c*********************************************************************c
c
         do 10 vector=oldnum+1,oldnum+added
            do 20 i=1,n
               hv(i,vector) = d(i)*v(i,vector) + hv(i,vector)
   20       continue
c*********************************************************************c
c                add the next row to the small matrix
c*********************************************************************c
c
            do 30 trial=1,vector
               small(vector,trial)=sdot(n,v(1,trial),1,hv(1,vector),1)
               small(trial,vector)=small(vector,trial)
   30       continue
   10    continue
c   
c*********************************************************************c
c                  print the the small matrix 
c*********************************************************************c
c
         oldnum=oldnum+added
         nnpnum=oldnum*(oldnum+1)/2
         if(debug) then
            title=' davidson:small matrix'
            call prntrm(title,small,oldnum,oldnum,maxvec,maxvec,iout)       
         endif
c
c*********************************************************************c
c                   solve the small eigenproblem
c*********************************************************************c
         call mmove(small,eigvec,oldnum,oldnum,maxvec,maxvec)
         call tred2(maxvec,oldnum,eigvec,root,resid,eigvec)
         call tql2(maxvec,oldnum,root,resid,eigvec,error)
         if(error.ne.0) then
            call lnkerr('davidson: error code from rsp')
         end if
c
c*********************************************************************c
c                   print results for this cycle
c*********************************************************************c
c
         cycle=cycle+1
         if (prnt) then
            write (iout,1) cycle, oldnum, totnum
            write(iout,2)
         end if
c
c*********************************************************************c
c               at this point we have a set of oldnum vectors
c               of which the lowest nroots vectors are what we 
c               want to converge.
c*********************************************************************c   
c
c
c*********************************************************************c
c               form the new "best" vectors.
c               copy them back to their proper vector slots.
c               at this point there are oldnum vectors which are        
c               the results of diagonalizing the small eigenproblem.
c               this number is not necessarily the number of roots
c               which we need.
c               then do the same for the products of h on the vectors.     
c*********************************************************************c
c
         call ebcxx(resid,v,eigvec,n,oldnum,oldnum,n,n,maxvec)
         call copy(resid,v,n*oldnum)
c
         call ebcxx(resid,hv,eigvec,n,oldnum,oldnum,n,n,maxvec)
         call copy(resid,hv,n*oldnum)
c
c*********************************************************************c
c                print results if requested
c*********************************************************************c
c
         if(debug) then
            title='new best vectors for cycle = '//itoc(cycle)
            call prntrm(title,v,n,oldnum,n,maxvec,iout)
         end if
c
         if(debug) then
            title='new best products for cycle = '//itoc(cycle)
            call prntrm(title,hv,n,oldnum,n,maxvec,iout)
         end if
c
c*********************************************************************c
c                we now form the residuals (h-e)c for the
c                all of the eigenvectors. 
c*********************************************************************c
c
         do 40 nroot=1,oldnum
            do 50 i=1,n
               resid(i,nroot)=hv(i,nroot)-root(nroot)*v(i,nroot)
 50         continue
            cvn(nroot)=sqrt(sdot(n,resid(1,nroot),1,resid(1,nroot),1))
 40      continue
         if(debug) then
            title=' (h-e)c '
            call prntrm(title,resid,n,oldnum,n,maxvec,iout)
         end if
c
c**********************************************************************c
c              check on convergence, and print results
c**********************************************************************c
         status='continue'
         do while(status.eq.'continue')
            pt=0
            ncnv=0
            maxdif=0.d0
            lower=lwr
            if(oldnum.ge.upr) then
               upper=upr
            else
               upper=oldnum
            endif               
            do 60 nroot=lower,upper
               if (cvn(nroot).gt.cnvgnc) then
                   pt=pt+1
                   ptroot(pt)=nroot
                   truth='no'
               else
                   ncnv=ncnv+1
                   truth='yes'
               end if
               if (cvn(nroot).gt.maxdif) then
                   maxdif=cvn(nroot)
               end if
               if (prnt) then
                  write (iout,3) nroot, root(nroot)+efixed, 
     1                                  cvn(nroot), truth
               endif
   60       continue
c
c**********************************************************************c
c          pt now tells us how may roots are unconverged and
c          the array ptroot has which roots they are.                  
c          converged vectors are frozen, counters are updated and
c          when all nattim are finished, the next group of vectors
c          are added to the list.
c**********************************************************************c
c 
            ntodo=upper-lower+1-ncnv
            i3=ntodo
c
            if (prnt) then
               endfile iout
               backspace iout
            end if
            if(ntodo.eq.0) then
               numcnv=numcnv+upper-lower+1
               if(numcnv.eq.nroots) then
c**********************************************************************c 
c                 all roots converged
c**********************************************************************c      
                  status='all roots converged'
                  return
               else
                  lwr=upr+1
                  upr=min(upr+nattim,nroots)
               endif
            else
               status='more'
            endif
         enddo               
         i4=numcnv
         status='continue'
         if (prnt) then
             endfile iout
             backspace iout
         end if
c
c**********************************************************************c
c               reform the small matrix in the new basis 
c**********************************************************************c
            do 70 nroot=1,oldnum
               do 80 trial=1,nroot
                  small(nroot,trial)
     1                               =
     2                           sdot(n,hv(1,nroot),1,v(1,trial),1)
                  small(trial,nroot) = small(nroot,trial)
c              write(iout,*) 'i',nroot,'j',trial,small(nroot,trial)
  80           continue
  70        continue
c
         added=0
         newtr=0
c         
c**********************************************************************c
c            create new trial vectors. this is done by solving the
c            equation for the error vector approximately.  the original
c            davidson procedure approximates the hamiltonian by its
c            diagonal, which is accurate for diagonally dominant matrices.
c            the approximate error vector is schmidt orthogonalized 
c            to the other vectors and then added to the trials if it is 
c            not linearly dependent on the other vectors.  the size of 
c            the norm is used as the criterion for linear dependence.      
c**********************************************************************c
c                      'new trials'
c
      else if (oper.eq.'new trials') then
c
c
         if(oldnum+ntodo.le.maxvec) then
            if(prres) then
               do 85 nroot=1,ntodo
                  pt=ptroot(nroot)
                  call copy(resid(1,pt),v(1,oldnum+nroot),n)
 85            continue
               title='residuals for cycle = '//itoc(cycle)
               call prntrm(title,v(1,oldnum+1),n,ntodo,n,ntodo,iout)
            endif
         endif      
         do 90 nroot=1,ntodo
            pt=ptroot(nroot)                  
c         
c*********************************************************************c
c           form the trial vector by:
c              1. solving the equation ( h -e ) error = residual at
c                 some approximate level. 
c              2. orthonormalizing it to the other trial vectors.
c              3. adding the new vector, if it is linearly
c                 independent to the trial set.   
c*********************************************************************c
c          
               call diaslv(resid(1,pt),root(pt),d,n)
c
c*********************************************************************c
c          orthonormalize to previous vectors
c*********************************************************************c
c
               do 200 i=1,oldnum+newtr
                  ovrlap=sdot(n,resid(1,pt),1,v(1,i),1)
                  call saxpy(n,-ovrlap,v(1,i),1,resid(1,pt),1)
                  norm=sqrt(sdot(n,resid(1,pt),1,resid(1,pt),1))
                  if (norm.lt.thresh) go to 90
  200          continue
c
               call sscal(n,1.0d+00/norm,resid(1,pt),1)
c
c
               if(schmdt) then  
                  do 300 i=1,oldnum+newtr
                     ovrlap=sdot(n,resid(1,pt),1,v(1,i),1)
                     call saxpy(n,-ovrlap,v(1,i),1,resid(1,pt),1)
                     norm=sqrt(sdot(n,resid(1,pt),1,resid(1,pt),1))
  300             continue
c
                  call sscal(n,1.0d+00/norm,resid(1,pt),1)
c
c               write(iout,*) 'second'
               endif
c               write(iout,*) 'output passed vector'
c               write(iout,*) (resid(ii,pt),ii=1,n)
               newtr=newtr+1
               call copy(resid(1,pt),v(1,oldnum+newtr),n)
   90       continue
         if(newtr.eq.0) then
            write(iout,6)
            status='none left'
         else
             status='done'                              
             write(iout,7) newtr
         endif 
         i3=newtr   
c
c                     'get vector'
c
      else if (oper.eq.'get vector') then
         if (i3.le.0.or.i3.gt.nroots) then
            status='error'
            return
         end if
c
         if (maxdif.ge.cnvgnc) then
            status='not converged'
         else
            status='ok'
         end if
c
c                 'restart'
c
      elseif(oper.eq.'restart') then
             if(i4.eq.nroots) then
                status='converged'
                return
             else
                i3=min(i4+i4,nroots+nroots,maxvec)
                i5=i3
                oldnum=0
                totnum=0
                newtr=0
                added=0
                cycle=0
                lwr=1
                upr=nattim
                write(iout,8) i3
                status='continue'
                return
             endif
c
c                 'cleanup'
c
      elseif(oper.eq.'cleanup') then
             nvecs=i3
c
c                       'finish'
c
      else if (oper.eq.'finish') then
         write(iout,9)
      else                    
       write(iout,*) status
         call lnkerr('bad call to davidson')
      end if
c
c
      return
    1 format (/,5x,' cycle = ',i4,1x,' number of vectors = ',i4,1x,
     1                     ' total number of vectors = ',i4)
    2 format(/,5x,'     root     ',3x,'     energy    ',7x,
     1            '      rms    ',3x,'   converged   ')
    3 format(9x,i3,10x,e15.8,6x,e15.8,8x,a3)
    4 format(/,1x,'residuals for cycle = ',i3)   
    5 format(/,1x,5e15.8)   
    6 format(/,15x,'there are no more trial vectors available')
    7 format(/,15x,'adding ',i3,' new vectors to the space')
    8 format(/,15x,'restarting with ',i4,' vectors')   
    9 format(/,1x,'eigenvalue calculation finished')         
      end


