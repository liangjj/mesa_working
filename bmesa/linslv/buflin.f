*deck buflin.f
c***begin prologue     buflin
c***date written       951229   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c***keywords           linear systems, variation-iteration
c***author             schneider, barry (nsf)
c***source             
c***purpose            driver for large linear systems solver using
c***references         a buffered file
c
c***routines called    
c***end prologue       buflin
      subroutine buflin(hbuf,ibuf,diag,vec,hvec,energy,rhs,dinv,
     1                 smat,svec,ipvt,thresh,cnverg,ops,n,nrhs,
     2                 nattim,mxvc,iter,lenbuf,ntot,incore)      
      implicit integer (a-z)
      real*8 hbuf, diag
      real*8 vec, hvec, svec, smat, rhs
      real*8 cnverg, thresh, rjunk, energy, dinv
      logical logkey, schmdt, incore, prnt
      character*(*) ops
      character*32 type
      character*80 dir, status
      character*8 prtflg
      dimension hbuf(lenbuf), ibuf(2,lenbuf), diag(n), rhs(n,nrhs)
      dimension smat(mxvc,mxvc,2), svec(mxvc,nrhs,2), dinv(n)
      dimension vec(n,mxvc), hvec(n,mxvc), ipvt(n), prnt(8)
      common/io/inp, iout
      write(iout,1) 
      call iosys ('read character "print flag" from rwf',-1,0,0,prtflg)      
      prnt(1)=logkey(ops,'print=linslv=all',.false.,' ')
      prnt(2)=logkey(ops,'print=linslv=vectors',.false.,' ')
      prnt(3)=logkey(ops,'print=linslv=overlaps',.false.,' ')
      prnt(4)=logkey(ops,'print=linslv=small-matrices',.false.,' ')
      prnt(5)=logkey(ops,'print=linslv=residuals',.false.,' ')
      prnt(6)=logkey(ops,'print=linslv=new-vectors',.false.,' ')
      prnt(7)=logkey(ops,'print=linslv=converged-solutions',
     1               .false.,' ')
      prnt(8)=logkey(ops,'print=linslv=best-solutions',.false.,' ')
      schmdt=logkey(ops,'linslv=two-schmidt-orthogonalizations',
     1              .false.,' ')
      type='matrix'
      north=1
      if(schmdt) then
         north=2
      endif   
      if (prtflg.eq.'minimum') then
         status='noprint'
      else
         status='print'
      end if
      ntrips=nrhs/nattim
      left=nrhs-ntrips*nattim
      if(left.ne.0) then
         ntrips=ntrips+1
      else
         left=nattim
      endif
      start=1
      do 40 trips=1,ntrips
         num2do=nattim
         if(trips.eq.ntrips) then
            num2do=left
         endif
         dir='initialize'
         call linslv(dir,type,vec,rhs(1,start),cnverg,thresh,
     1               rjunk,rjunk,rjunk,mxvc,n,iout,num2do,nvec,
     2               prnt,schmdt)
         status=''
         ntrial=num2do
         do while(status.ne.'converged')
            nold=nvec
            if(status.ne.'converged') then
               dir='generate basis functions'
c              number of trial vectors, ntrial and total number
c              of vectors, nvec returned after the orthonormalization
c              is performed.  
               call linslv(dir,status,vec,rjunk,rjunk,rjunk,
     1                     smat(1,1,2),rjunk,rjunk,ntrial,n,
     2                     mxvc,mxvc,nvec,prnt,schmdt)
               if(status.eq.'iteration limit exceeded') then
                  write(iout,2) status, nvec
                  return 
               endif
               if(ntrial.ne.0) then
c                 operate on the trial vectors with hamiltonian to prepare
c                 for projection onto vector space and solution of
c                 set of trial linear equations.   
                  call honv(hbuf,ibuf,vec(1,nold+1),hvec(1,nold+1),n,
     1                      ntrial,lenbuf,ntot,incore,prnt)
                  if(type.eq.'matrix') then
                     call gvmmul(diag,vec(1,nold+1),hvec(1,nold+1),
     1                           hvec(1,nold+1),n,ntrial)   
                  else
                     call vmmul(dinv,hvec(1,nold+1),hvec(1,nold+1),n,
     1                          ntrial)     
                  endif
                  dir='solve equations'
c                 solve the trial set of equations, check for convergence
c                 and prepare for the next set of trial vectors which are
c                 taken from the previous set of iterates.  the call ensures
c                 that the trial set never gets bigger than the maximum.
                  call linslv (dir,status,vec,hvec,smat(1,1,2),
     1                         smat(1,1,1),svec(1,1,2),svec(1,1,1),
     2                         rhs(1,start),ipvt,n,iter,ntrial,nvec,
     3                         prnt,schmdt)
                  if(status.eq.'starting to diverge') then                
                     write(iout,3) nold
                     return
                  elseif(status.eq.'continue iterations') then
c                    copy the trial set into the proper vector slots and
c                    then return to the orthogonalization step.    
                     dir='new trials'
                     call linslv(dir,status,vec(1,nvec+1),
     1                            hvec(1,nold+1),rjunk,rjunk,rjunk,
     2                            rjunk,rjunk,junk,n,junk,ntrial,junk,
     3                            prnt,schmdt)
                  endif
               else
                  write(iout,4) nold
                  return
               endif
            else      
               call linslv('final best solution',status,vec,rjunk,
     1                      rjunk,rjunk,svec(1,1,2),rjunk,
     2                      rhs(1,start),junk,n,mxvc,mxvc,nvec,
     3                      prnt,schmdt)
               write(iout,2) status, nvec
               return
            endif 
         enddo
         start=start+num2do                                         
 40   continue   
      return
 1    format(/,5x,'enter iterative linear system solver')
 2    format(/,5x,'message from linslv:',/,26x,a32,/,5x,
     1            'number of vectors used= ',i5)
 3    format(/,5x,'solution starting to diverge.  will quit',/,5x,
     1            'number of vectors used = ',i5)     
 4    format(/,5x,'no more vectors.  we are done.',/,10x,
     1            'number of vectors used = ',i5)     
      end       
