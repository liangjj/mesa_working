*deck rdiis.f
c***begin prologue     rdiis
c***date written       960925   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c***keywords           diis, extrapolate,least-squares
c***author             komornici, andrew (polyatomics)
c***source             math
c***purpose            extrapolate on a set of objects
c***description        objects(obj) are extrpolated using the error vector
c***                   (xerr) and the diis methodology of pulay.
c
c                      pvec(n,maxit)    iterate vectors.
c                      dlpvec(n,maxit)  error vectors.  
c                      vec(n)           extrapolated vector.                 
c                      ipvt             index vector for sgefa
c                      iter             current iteration.
c                      n                length of object vector
c                      maxit            maximum number of iterations
c
c***references         pulay, chem. phys. letts., 73, 393(1980)
c
c***routines called    sdot, sgefa, sgesl
c***end prologue       rdiis
      subroutine rdiis(pvec,dlpvec,vec,b,btmp,sol,ipvt,iter,n,maxit)
      implicit integer (a-z)
      real*8 pvec, dlpvec, vec, b, btmp, sol
      real*8 zero, one, mone
      real*8 sdot
      real*8 sum, error, xerror
      dimension pvec(n,maxit+1), dlpvec(n,maxit+1), vec(n)
      dimension b(maxit+1,maxit+1), sol(maxit+1), btmp(maxit+1,maxit+1)
      dimension ipvt(maxit+1)
      parameter (zero=0.0d0, one=1.0d0, mone=-1.d0)
      common/io/inp, iout
c
c        *****  update iteration counters:                      *****
      if(iter.lt.1) then
         write(iout,1)
         return
      elseif(iter.eq.1) then
         write(iout,2)
         do 10 i=1,n
            dlpvec(i,1) = pvec(i,2) - pvec(i,1)
 10      continue
         b(1,1) = sdot(n,dlpvec(1,1),1,dlpvec(1,1),1)
         return
      elseif(iter.gt.1.and.iter.le.maxit) then
         error = sqrt(sdot(n,dlpvec(1,iter),1,dlpvec(1,iter),1))
c
c        *****  update the current pulay b matrix:              *****
         do 20 ii=1,iter
            b(ii,iter) = sdot(n,dlpvec(1,ii),1,dlpvec(1,iter),1)
            b(iter,ii) = b(ii,iter)
 20      continue
c        *****  transfer into temp array:                       *****
         do 30 i=1,iter
            do 40 j=1,i
               btmp(i,j) = b(i,j)
               btmp(j,i) = b(j,i)
 40         continue
 30      continue   
c
c        *****  put in remaining row and column:                *****
c        *****  construct right-hand-side:                      *****
         do 50 i=1,iter
            btmp(i,iter+1) = mone
            btmp(iter+1,i) = mone
 50      continue   
         btmp(iter+1,iter+1) = zero
         call rzero(sol,iter+1)
         sol(iter+1) = mone      
c        *****  solve set of simultaneous equations:            *****
         call sgefa(btmp,maxit+1,iter+1,ipvt,info)
         if(info.ne.0) then
            call lnkerr('error in solution of diis linear equations')
         endif
         call sgesl(btmp,maxit+1,iter+1,ipvt,sol,0)
c        *****  check that solution vector adds to one          *****
         sum = zero
         do 60 i=1,iter
            sum = sum + sol(i)
 60      continue
         write(iout,3) iter, sum
c        *****  extrapolate the objects:                        *****
         call rzero(vec,n)
         call ebcxx(vec,dlpvec,sol,n,iter,1,n,n,n)
         xerror = sqrt(sdot(n,vec,1,vec,1))
c        *****  check on magnitude of new error vector:         *****
c        *****  if this test is satisfied, extrapolation        *****
c        *****  may have gone south.                            *****
         if(xerror.gt.error) then
            write(iout,4)  iter, error, xerror
            return
         else
            call copy(vec,dlpvec(1,iter+1),n)
            call rzero(vec,n)
            call ebcxx(vec,pvec,sol,n,iter,1,n,n,n)
         endif
      elseif(iter.gt.maxit) then
         call lnkerr('number iterations exceeded in diis')
      endif
      return 
 1    format(/,10x,'iteration outside range.  return')
 2    format(/,10x,'begin rdiis extrapolation')
 3    format(/,10x,'iteration = ',i3,1x,'coefficient sum = ',e15.8)
 4    format(10x,'end of rdiis extrapolation:  ',/,
     &       10x,'iter: ',i3,/,
     &       10x,'norm of input  error vector:  ',1pe12.4,/,
     &       10x,'norm of output error vector:  ',1pe12.4,/ )
      end
