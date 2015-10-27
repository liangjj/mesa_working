*deck diis.f
c***begin prologue     diis
c***date written       961209   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c***keywords           diis, pulay, extrapolation
c***author             schneider, barry (nsf)
c***source             math
c***purpose            extrapolate on a set of objects using the error
c***                   vector to extrapolate
c***                   
c***                   
c***references         peter pulay
c
c***routines called    
c***end prologue       diis
      subroutine diis(p,resid,pvec,dpvec,b,btmp,sol,norm0,error,
     1                cnverg,ipvt,iter,nobj,nerr,trunc,maxit,
     2                prnt,flag,incore)
      implicit integer (a-z)
      real*8 p, resid, pvec, dpvec, b, btmp, sol
      real*8 zero, one, mone, sdot, tstone, norm0, error, cnverg
      character*80 title
      character*3 itoc
      logical prnt, flag, incore
      dimension p(nobj), resid(nerr), pvec(nobj,*)
      dimension dpvec(nerr,*)
      dimension b(maxit+1,maxit+1), btmp(maxit+1,maxit+1)
      dimension sol(maxit+1), ipvt(maxit+1)
      common/io/inp, iout
      data zero, one, mone/ 0.d0, 1.d0, -1.d0/
      if(iter.gt.maxit) then
         write(iout,1)
         call lnkerr('exceeds maxit')
      else
         do 10 i=1,iter-1
            do 20 j=1,i
               btmp(i,j) = b(i,j)
               btmp(j,i) = b(j,i)
 20         continue
 10      continue
         if(incore) then
            do 30 i=1,iter-1
               b(i,iter) = sdot(nerr,dpvec(1,i),1,dpvec(1,iter),1)
               b(iter,i) = b(i,iter)
               btmp(i,iter) = b(i,iter)
               btmp(iter,i) = b(iter,i)
 30         continue   
            b(iter,iter) = sdot(nerr,dpvec(1,iter),1,dpvec(1,iter),1)
         else
            do 40 i=1,iter-1
               call iosys('read real "error vector iteration = '
     1                    //itoc(i)//'" from lamdat',nerr,
     2                      resid,0,' ')
               b(i,iter) = sdot(nerr,resid,1,dpvec,1)
               b(iter,i) = b(i,iter)
               btmp(i,iter) = b(i,iter)
               btmp(iter,i) = b(iter,i)
 40         continue
            b(iter,iter) = sdot(nerr,dpvec,1,dpvec,1)   
         endif
         btmp(iter,iter) = b(iter,iter)
         do 50 i=1,iter
            b(i,iter+1) = mone
            b(iter+1,i) = mone
            btmp(i,iter+1) = mone
            btmp(iter+1,i) = mone
 50      continue
         b(iter+1,iter+1) = zero   
         btmp(iter+1,iter+1) = zero   
         call rzero(sol,iter+1)
         sol(iter+1) = mone
         if(prnt) then
            title='b matrix for iteration = '//itoc(iter)
c            call prntrm(title,b,iter+1,iter+1,maxit+1,maxit+1,iout)
            call prntfm(title,b,iter+1,iter+1,maxit+1,maxit+1,iout)
         endif
         error=1.d+10
         call sgefa(btmp,maxit+1,iter+1,ipvt,info)
         if(info.ne.0) then
            call lnkerr('error in linear solve:singular matrix')
         endif
         call sgesl(btmp,maxit+1,iter+1,ipvt,sol,0)
         if(prnt) then
            title='solution vector for iteration = '//itoc(iter)
c            call prntrm(title,sol,iter+1,1,maxit+1,1,iout)
            call prntfm(title,sol,iter+1,1,maxit+1,1,iout)
         endif
         tstone = zero
         do 60 i=1,iter
            tstone = tstone +sol(i)
 60      continue
         if(incore) then
            call ebc(p,pvec,sol,nobj,iter,1)
         else
            call smul(p,pvec,sol(iter),nobj)
            do 70 i=1,iter-1
               call iosys('read real "object vector iteration = '
     1                    //itoc(i)//'" from lamdat',nobj,
     2                      pvec,0,' ')
               call saxpy(nobj,sol(i),pvec,1,p,1)
 70         continue
         endif   
         if(prnt) then
            title='extrapolated object for iteration = '//itoc(iter)
c            call prntrm(title,p,nobj,1,nobj,1,iout)
            call prntfm(title,p,nobj,1,nobj,1,iout)
         endif
         if(incore) then
            call ebc(resid,dpvec,sol,nerr,iter,1)
         else
            call smul(resid,dpvec,sol(iter),nerr)
            do 80 i=1,iter-1
               call iosys('read real "error vector iteration = '
     1                    //itoc(i)//'" from lamdat',nerr,
     2                      dpvec,0,' ')
               call saxpy(nerr,sol(i),dpvec,1,resid,1)
 80         continue
         endif
         if(prnt) then
            title='extrapolated error for iteration = '//itoc(iter)
c            call prntrm(title,resid,nerr,1,nerr,1,iout)
            call prntfm(title,resid,nerr,1,nerr,1,iout)
         endif
         error = 0.d0
         do 90 i=1,nerr
            error = max(error,abs(resid(i)))
 90      continue
c         error=error/(norm0*norm0)   
         error=error/norm0
         write(iout,2) iter, tstone, error
         if(error.le.cnverg) then
            flag=.false.
            return
         else
            if(iter.lt.trunc) then
               flag=.false.
               return
            else
               call copy(p,pvec,nobj)
               flag=.true.
               write(iout,3)
               return
            endif            
         endif
      endif
      return
 1    format(/,5x,'iteration exceeds ',i3,' quit')
 2    format(/,5x,'iteration  = ',i4,/,5x,
     1            'lamda      = ',e15.8,/,5x,
     2            'rms error  = ',e15.8)    
 3    format(/,5x,'restart')
      end       
