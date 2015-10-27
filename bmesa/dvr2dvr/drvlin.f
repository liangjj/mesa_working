*deck drvlin.f
c***begin prologue     drvlin
c***date written       951229   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c***keywords           linear systems, variation-iteration
c***author             schneider, barry (nsf)
c***source             
c***purpose            driver large linear systems solver.
c***references         
c
c***routines called    
c***end prologue       drvlin
      subroutine drvlin(hbuf,ibuf,diag,vec,hvec,trvec,energy,rhs,dinv,
     1                 smat,svec,ipvt,thresh,cnverg,ops,n,nrhs,
     2                 nattim,mxvc,iter,lenbuf,ntot,nvec,ntrial,
     3                 incore,smooth,grid)      
      implicit integer (a-z)
      real*8 hbuf, diag
      real*8 vec, hvec, trvec, svec, smat, rhs
      real*8 cnverg, thresh, rjunk, energy, dinv
      logical logkey, schmdt, incore, prnt
      character*(*) ops, smooth
      character*32 type, chrkey
      character*80 dir, status
      character*8 prtflg
      character*2 itoc
      dimension hbuf(lenbuf), ibuf(2,lenbuf), diag(n), rhs(n,nrhs)
      dimension smat(mxvc,mxvc,2), svec(mxvc,nrhs,2), dinv(n)
      dimension vec(n,mxvc), hvec(n,mxvc), trvec(n,ntrial)
      dimension ipvt(n), prnt(8)
      common/io/inp, iout
      call iosys('create real "solutions for grid '//
     1            itoc(grid)//'" on ham',-1,0,0,' ')
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
      type=chrkey(ops,'linslv=form','one-minus-matrix',' ')
      if(type.eq.'matrix') then
         do 10 i=1,n
            dinv(i)=1.d0
 10      continue      
      elseif(type.eq.'one-plus-matrix') then
         do 20 i=1,n
            dinv(i)=1.d0/diag(i)
 20      continue
      elseif(type.eq.'one-minus-matrix') then
         do 30 i=1,n
            dinv(i) = 1.d0/( energy - diag(i) )
 30      continue
      else
         call lnkerr('error in form of equations')
      endif         
      call vmmul(dinv,rhs,rhs,n,nrhs)   
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
         ntr=ntrial  
         call linslv(dir,type,vec,trvec,cnverg,thresh,
     1               rjunk,rjunk,rjunk,mxvc,n,iout,num2do,ntr,
     2               prnt,schmdt)
         nvec=ntr
         status=''
         ntr=num2do
         do while(status.ne.'converged')
            nold=nvec
            if(status.ne.'converged') then
               dir='generate basis functions'
c              number of trial vectors, ntr and total number
c              of vectors, nvec returned after the orthonormalization
c              is performed.  
               call linslv(dir,status,vec,rjunk,rjunk,rjunk,
     1                     smat(1,1,2),rjunk,rjunk,ntr,n,
     2                     mxvc,mxvc,nvec,prnt,schmdt)
               if(status.eq.'iteration limit exceeded') then
                  write(iout,2) status, nvec
                  call iosys('endfile "solutions for grid '//itoc(grid)
     1                       //'" on ham',0,0,0,' ')   
                  return 
               endif
               if(ntr.ne.0) then
c                 operate on the trial vectors with hamiltonian to prepare
c                 for projection onto vector space and solution of
c                 set of trial linear equations.   
                  call honv(hbuf,ibuf,diag,vec(1,nold+1),
     1                      hvec(1,nold+1),n,ntr,lenbuf,ntot,
     2                      incore,smooth,prnt,grid)
                  if(type.eq.'matrix') then
                     call gvmmul(diag,vec(1,nold+1),hvec(1,nold+1),
     1                           hvec(1,nold+1),n,ntr)   
                  else
                     call vmmul(dinv,hvec(1,nold+1),hvec(1,nold+1),n,
     1                          ntr)     
                  endif
                  dir='solve equations'
c                 solve the trial set of equations, check for convergence
c                 and prepare for the next set of trial vectors which are
c                 taken from the previous set of iterates.  the call ensures
c                 that the trial set never gets bigger than the maximum.
                  call linslv (dir,status,vec,hvec,smat(1,1,2),
     1                         smat(1,1,1),svec(1,1,2),svec(1,1,1),
     2                         rhs(1,start),ipvt,n,iter,ntr,nvec,
     3                         prnt,schmdt)
                  if(status.eq.'starting to diverge') then                
                     write(iout,3) nold
                     call iosys('endfile "solutions for grid '//
     1                           itoc(grid)//'" on ham',0,0,0,' ')   
                     return
                  elseif(status.eq.'continue iterations') then
c                    copy the trial set into the proper vector slots and
c                    then return to the orthogonalization step.    
                     dir='new trials'
                     call linslv(dir,status,vec(1,nvec+1),
     1                            hvec(1,nold+1),rjunk,rjunk,rjunk,
     2                            rjunk,rjunk,junk,n,junk,ntr,junk,
     3                            prnt,schmdt)
                  endif
               else
                  write(iout,4) nold
                  call iosys('endfile "solutions for grid '//
     1                        itoc(grid)//'" on ham',0,0,0,' ')   
                  return
               endif
            else      
               call linslv('final best solution',status,vec,rjunk,
     1                      rjunk,rjunk,svec(1,1,2),rjunk,
     2                      rhs(1,start),junk,n,mxvc,mxvc,nvec,
     3                      prnt,schmdt)
               write(iout,2) status, nvec
               call iosys('endfile "solutions for grid '//
     1                     itoc(grid)//'" on ham',0,0,0,' ')   
               return
            endif 
         enddo
         call iosys('write real "solutions for grid '//itoc(grid)
     1               //'" on ham without rewinding',num2do*n,
     2                    rhs(1,start),0,' ')
         start=start+num2do                                         
 40   continue
      call iosys('endfile "solutions for grid '//itoc(grid)
     1            //'" on ham',0,0,0,' ')   
      return
      call iosys('rewind all on ham read-and-write',0,0,0,' ')
 1    format(/,5x,'enter iterative linear system solver')
 2    format(/,5x,'message from linslv:',/,26x,a32,/,5x,
     1            'number of vectors used= ',i5)
 3    format(/,5x,'solution starting to diverge.  will quit',/,5x,
     1            'number of vectors used = ',i5)     
 4    format(/,5x,'no more vectors.  we are done.',/,10x,
     1            'number of vectors used = ',i5)     
      end       
