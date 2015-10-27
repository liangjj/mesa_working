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
      subroutine diis(p,resid,pvec,dpvec,b,btmp,sol,error,
     1                cnverg,ipvt,iter,nobj,nerr,trunc,maxit,prnt,flag)
      implicit integer (a-z)
      complex*16 p, resid, pvec, dpvec, b, btmp, sol
      complex*16 cdotc, zeroc, onec, monec, tstone
      real*8 error, cnverg
      character*80 title
      character*2 itoc
      logical prnt, flag
      dimension p(nobj), resid(nerr), pvec(nobj,maxit+1)
      dimension dpvec(nerr,maxit+1)
      dimension b(maxit+1,maxit+1), btmp(maxit+1,maxit+1)
      dimension sol(maxit+1), ipvt(maxit+1)
      common/io/inp, iout
      data zeroc, onec, monec / (0.d0,0.d0), (1.d0,0.d0), (-1.d0,0.d0) /
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
         do 30 i=1,iter-1
            b(i,iter) = cdotc(nerr,dpvec(1,i),1,dpvec(1,iter),1)
            b(iter,i) = cdotc(nerr,dpvec(1,iter),1,dpvec(1,i),1)
            btmp(i,iter) = b(i,iter)
            btmp(iter,i) = b(iter,i)
 30      continue   
         b(iter,iter) = cdotc(nerr,dpvec(1,iter),1,dpvec(1,iter),1)
         btmp(iter,iter) = b(iter,iter)
         do 40 i=1,iter
            b(i,iter+1) = monec
            b(iter+1,i) = monec
            btmp(i,iter+1) = monec
            btmp(iter+1,i) = monec
 40      continue
         b(iter+1,iter+1) = zeroc   
         btmp(iter+1,iter+1) = zeroc   
         call czero(sol,iter+1)
         sol(iter+1) = monec
         if(prnt) then
            title='b matrix for iteration = '//itoc(iter)
            call prntcm(title,b,iter+1,iter+1,maxit+1,maxit+1,iout)
         endif
         error=1.d+10
         call cgefa(btmp,maxit+1,iter+1,ipvt,info)
         if(info.ne.0) then
            call lnkerr('error in linear solve:singular matrix')
         endif
         call cgesl(btmp,maxit+1,iter+1,ipvt,sol,0)
         if(prnt) then
            title='solution vector for iteration = '//itoc(iter)
            call prntcm(title,sol,iter+1,1,maxit+1,1,iout)
         endif
         tstone = zeroc
         do 50 i=1,iter
            tstone = tstone +sol(i)
 50      continue
         call cebc(p,pvec,sol,nobj,iter,1)
         if(prnt) then
            title='extrapolated object for iteration = '//itoc(iter)
            call prntcm(title,p,nobj,1,nobj,1,iout)
         endif
         call cebc(resid,dpvec,sol,nerr,iter,1)
         error = 0.d0
         do 60 i=1,nerr
            error = max(error,abs(resid(i)))
 60      continue
         write(iout,2) iter, tstone, error
         if(error.le.cnverg) then
            flag=.false.
            return
         else
            if(iter.lt.trunc) then
               flag=.false.
               return
            else
               call cc2opy(p,pvec(1,1),nobj)
               flag=.true.
               write(iout,3)
               return
            endif            
         endif
      endif
      return
 1    format(/,5x,'iteration exceeds ',i3,' quit')
 2    format(/,5x,'iteration  = ',i4,/,5x,
     1            'lamda      = ',2e15.8,/,5x,
     2            'rms error  = ',e15.8)    
 3    format(/,5x,'restart')
      end       
