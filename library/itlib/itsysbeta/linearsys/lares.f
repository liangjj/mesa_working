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
      subroutine lares(vec,hvec,b,coef,rhs,energy,cnverg,resid,
     1                 maxerr,t,list,n,m,nrhs,ncon,addvec,maxvec,
     2                 it,prnt)
      implicit integer (a-z)
      real*8 vec, hvec, b, coef, rhs, energy, resid, t
      real*8 cnverg, ddot, maxerr, err, zero, one
      character*16 status
      character*80 title
      character*4 itoc
      logical prnt
      dimension vec(n,maxvec), hvec(n,maxvec)
      dimension b(maxvec,nrhs), coef(maxvec,nrhs), rhs(n,nrhs)
      dimension resid(n,nrhs), t(n,nrhs), prnt(4), list(*)
      common/io/inp, iout
      data zero, one /0.d0,1.d0/
c
c        first express the solution vectors, Psi, in the original basis
c        and put them in t.
c       
      call dgemm('n','n',n,m,m,one,vec,n,coef,maxvec,zero,t,n)      
c
c        calculate H*Psi and put in resid.
c
      call dgemm('n','n',n,m,m,one,hvec,n,coef,maxvec,zero,resid,n)
c
c
      if(prnt(1)) then
         title='information for iteration = '//itoc(it)
         write(iout,1) title
      endif
c
c        final residual
c
      do 10 i=1,nrhs
         do 20 j=1,n
            resid(j,i) = rhs(j,i) - ( energy*t(j,i) - resid(j,i) )
 20      continue
 10   continue
      if(prnt(2)) then
         title='residuals iteration = '//itoc(it)
         call prntfm(title,resid,n,nrhs,n,maxvec,iout)
      endif
      addvec=0
      ncon=0
      maxerr=0.d0
      do 30 i=1,nrhs
         err = sqrt (ddot(n,resid(1,i),1,
     1                      resid(1,i),1) )
         maxerr=max(err,maxerr)
         if(err.le.cnverg) then
            status='converged'
            ncon=ncon+1
	    write(70) list(i)
            write(70) (t(j,i), j=1,n)
            write(iout,2) list(i), err, status
         else
            status='unconverged'
            addvec=addvec+1
            call copy(rhs(1,i),rhs(1,addvec),n)
            call copy(resid(1,i),resid(1,addvec),n)
            call copy(b(1,i),b(1,addvec),maxvec)
            list(addvec)=i
            write(iout,2) list(i), err, status
         endif
 30   continue
      if(ncon.eq.nrhs) then
         write(iout,3)
      else
         write(iout,4) ncon, addvec
      endif
      return
 1    format(/,5x,a80)         
 2    format(/,5x,'solution        = ',i4,/,5x,
     1            'rms error       = ',e15.8,/,5x,
     2            'status          = ',a16)
 3    format(/,15x,'all solutions are converged.  we are done')
 4    format(/,5x,'number converged solutions   = ',i4,/,5x,
     1            'number unconverged residuals = ',i4)
      end



