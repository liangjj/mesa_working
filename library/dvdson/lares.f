*deck lares.f
c***begin prologue     lares
c***date written       980420   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c***keywords           residual calculation
c***author             schneider, barry (nsf)
c***source             
c***                   
c***references         
c
c***routines called    
c***end prologue       lares
      subroutine lares(vec,hvec,coef,rhs,energy,scale,cnverg,resid,
     1                 maxerr,t1,t2,n,m,nrhs,ncon,addvec,maxvec,
     2                 it,prnt)
      implicit integer (a-z)
      real*8 vec, hvec, coef, rhs, energy, resid, t1, t2
      real*8 scale, cnverg, sdot, maxerr, err
      character*16 status
      character*80 title
      character*4 itoc
      logical prnt
      dimension vec(n,*), hvec(n,*), coef(maxvec,nrhs), rhs(n,nrhs)
      dimension resid(n,*), t1(n,*), t2(n,*), prnt(4)
      common/io/inp, iout
c
c        first express the solution vectors in the original basis
c        and put them t1
c       
      call ebcxx(t1,vec,coef,n,m,nrhs,n,n,maxvec)
c
c        now the effect of the hamiltonian on the solution vectors
c        and put them in t2
c
      call ebcxx(t2,hvec,coef,n,m,nrhs,n,n,maxvec)
      if(prnt(1)) then
         title='information for iteration = '//itoc(it)
         write(iout,1) title
      endif               
      do 10 i=1,nrhs
         do 20 j=1,n
            resid(j,i) = rhs(j,i) - ( energy*t1(j,i) - t2(j,i) )
 20      continue
 10   continue
      if(prnt(2)) then
         title='residuals iteration = '//itoc(it)
         call prntfm(title,resid,n,nrhs,n,maxvec,iout)
      endif
      if(prnt(3)) then
         title='solution vectors iteration = '//itoc(it)
         call prntfm(title,t1,n,nrhs,n,maxvec,iout)
      endif
      if(prnt(4)) then
          title='hamiltonian on solution vectors iteration = '
     1          //itoc(it)
         call prntfm(title,t2,n,nrhs,n,maxvec,iout)
      endif
      addvec=0
      ncon=0
      maxerr=0.d0
      do 30 i=1,nrhs
         err = scale*sqrt (sdot(n,resid(1,i),1,
     1                            resid(1,i),1) )
         maxerr=max(err,maxerr)
         if(err.le.cnverg) then
            status='converged'
            ncon=ncon+1
            call iosys('write real "solution for right hand side = '
     1                 //itoc(i)//'" to ham',n,t1(1,i),0,' ')
            if(prnt(1)) then
               write(iout,2) i, err, status
            endif                        
         else
            status='unconverged'
            addvec=addvec+1
            call copy(resid(1,i),resid(1,addvec),n)
            if(prnt(1)) then
               write(iout,2) i, err, status
            endif     
         endif
 30   continue
      test=m+addvec
      if(test.gt.maxvec) then
         addvec=maxvec-m
         write(iout,3) addvec
      endif    
      return
 1    format(/,5x,a80)         
 2    format(/,5x,'solution        = ',i4,/,5x,
     1            'rms error       = ',e15.8,/,5x,
     2            'status          = ',a16)
 3    format(/,5x,'number of added vectors will exceed maxvec',/,5x,
     1            'number of vectors actually added = ',i4)                   
      end

