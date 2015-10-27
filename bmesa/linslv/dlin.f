*deck dlin.f
c***begin prologue     dlin
c***date written       951229   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c***keywords           linear systems, variation-iteration
c***author             schneider, barry (nsf)
c***source             
c***purpose            driver for large linear systems solver using
c***references         raw matrix elements
c
c***routines called    
c***end prologue       dlin
      subroutine dlin(ham,v,u,eig,vec,hvec,energy,rhs,d,dinv,
     1                 smat,svec,ipvt,t,s,thresh,cnverg,ops,n,nrhs,
     2                 nattim,mxvc,iter,precon)      
      implicit integer (a-z)
      real*8 ham, v, u, eig, vec, hvec, d, dinv, smat, svec, rhs
      real*8 cnverg, thresh, rjunk, energy, t, s
      logical logkey, schmdt, prnt
      character*(*) ops, precon
      character*32 type
      character*80 dir, status
      character*8 prtflg
      dimension ham(n,n), v(n,n), u(n,n), eig(n)
      dimension vec(n,mxvc), hvec(n,mxvc), rhs(n,nrhs), d(n)
      dimension dinv(n), smat(mxvc,mxvc,2), svec(mxvc,nrhs,2)
      dimension t(n,nrhs), s(n,nrhs), ipvt(n), prnt(8)
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
c
c     calculate the diagonal scaling matrix in the h0 representation or use
c     the raw diagonal.
c
      if (precon.eq.'diagonalize-h0') then
          do 10 i=1,n
             dinv(i) = 1.d0/ (energy - eig(i) )
 10       continue
c
c     transform the right hand side from the dvr to the h0 representation
c
         call ebtc(t,u,rhs,n,n,nrhs)
         call copy(t,rhs,n*nrhs)
c
c     scale by the diagonal preconditioner
c
         do 20 i=1,n
            do 30 j=1,nrhs
               rhs(i,j) = rhs(i,j) * dinv(i)
 30         continue
 20      continue 
c
c     transform the right hand side back to the dvr representation
c
         call ebc(t,u,rhs,n,n,nrhs)
         call copy(t,rhs,n*nrhs)  
      elseif(precon.eq.'invert-h0') then
         do 40 i=1,n
            do 50 j=1,n
               ham(i,j) = - ham(i,j)
 50         continue
            ham(i,i) = energy + ham(i,i)
 40      continue   
         call sgefa(ham,n,n,ipvt,info)
         do 60 i=1,nrhs
            call sgesl(ham,n,n,ipvt,rhs(1,i),0)
 60      continue   
      else
         call preph(ham,v,d,n)
         do 70 i=1,n
            dinv(i) = 1.d0/ ( energy - d(i) )
 70      continue
c
c     scale by the diagonal
c
         do 80 i=1,n
            do 90 j=1,nrhs
               rhs(i,j) = rhs(i,j) * dinv(i)
 90         continue
 80      continue   
      endif             
      type='one-minus-matrix'
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
      do 100 trips=1,ntrips
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
c
c                 prepare a new trial vector of the preconditioned
c                 system.
c
                  if(precon.eq.'diagonalize-h0') then 
                     call preslv(dinv,v,vec(1,nold+1),
     1                           hvec(1,nold+1),u,t,s,n,
     2                           n,ntrial,prnt)
                  elseif(precon.eq.'invert-h0') then
                     call h0slv(ham,v,vec(1,nold+1),hvec(1,nold+1),
     1                          ipvt,n,ntrial)
                  else
                     call dslv(dinv,ham,vec(1,nold+1),hvec(1,nold+1),
     1                         n,ntrial,prnt)
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
 100  continue   
      return
 1    format(/,5x,'enter iterative linear system solver')
 2    format(/,5x,'message from linslv:',/,26x,a32,/,5x,
     1            'number of vectors used= ',i5)
 3    format(/,5x,'solution starting to diverge.  will quit',/,5x,
     1            'number of vectors used = ',i5)     
 4    format(/,5x,'no more vectors.  we are done.',/,10x,
     1            'number of vectors used = ',i5)     
      end       
