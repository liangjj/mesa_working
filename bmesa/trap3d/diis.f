*deck diis.f
c***begin prologue     diis
c***date written       961209   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c***keywords           diis, pulay, extrapolation
c***author             schneider, barry (nsf)
c***source             math
c***purpose            extrapolate on a set of objects using the error
c***                   vector to extrapolate without peripheral storage.
c***                   
c***                   
c***references         peter pulay
c
c***routines called    
c***end prologue       diis
      subroutine diis(pcur,pold,ham01,vnl,psi,hpsi,
     1                b,btmp,hh,pp,ph,hp,sol,error,ipvt,iter,n,
     2                maxit,prnt)
      implicit integer (a-z)
      real*8 pcur, pold, ham01, vnl, psi, hpsi
      real*8 b, btmp, hh, pp, ph, hp, sol, error
      real*8 zero, one, mone, tstone
      character*80 title
      character*3 itoc
      logical prnt
      dimension pcur(n*n), pold(n*n)
      dimension ham01(n,n), vnl(n,*), psi(n,*), hpsi(n,*)
      dimension b(maxit+1,maxit+1), btmp(maxit+1,maxit+1)
      dimension sol(maxit+1), ipvt(maxit+1)
      dimension hh(maxit,maxit), pp(maxit,maxit), ph(maxit,maxit)
      dimension hp(maxit,maxit)
      common/io/inp, iout
      data zero, one, mone/ 0.d0, 1.d0, -1.d0/
      nobj=n*n
      nerr=n*n
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
c                         th
c              calculate i  error matrix
c
            b(i,iter) = 2.d0*( hh(iter,i)*pp(iter,i) 
     1                                  - 
     2                         hp(iter,i)*ph(iter,i) )
            b(iter,i) = b(i,iter)
            btmp(i,iter) = b(i,iter)
            btmp(iter,i) = b(iter,i)
 30      continue
          b(iter,iter) = 2.d0*( hh(iter,iter)*pp(iter,iter) 
     1                                   - 
     2                         hp(iter,iter)*ph(iter,iter) )
         btmp(iter,iter) = b(iter,iter)
         do 40 i=1,iter
            b(i,iter+1) = mone
            b(iter+1,i) = mone
            btmp(i,iter+1) = mone
            btmp(iter+1,i) = mone
 40      continue
         b(iter+1,iter+1) = zero   
         btmp(iter+1,iter+1) = zero   
         call rzero(sol,iter+1)
         sol(iter+1) = mone
         if(prnt) then
            title='b matrix for iteration = '//itoc(iter)
c            call prntrm(title,b,iter+1,iter+1,maxit+1,maxit+1,iout)
            call prntfm(title,b,iter+1,iter+1,maxit+1,maxit+1,iout)
         endif
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
         do 50 i=1,iter
            tstone = tstone +sol(i)
 50      continue
         call newfck(pcur,ham01,vnl,sol,n,iter)
         call copy(pcur,pold,n*n)
         if(prnt) then
            title='extrapolated object for iteration = '//itoc(iter)
            call prntfm(title,pcur,nobj,1,nobj,1,iout)
         endif
         call newerr(psi,hpsi,sol,error,n,iter)
         write(iout,2) iter, tstone
      endif
      return
 1    format(/,5x,'iteration exceeds ',i3,' quit')
 2    format(/,5x,'iteration  = ',i4,/,5x,
     1            'lamda      = ',e15.8)
      end       
