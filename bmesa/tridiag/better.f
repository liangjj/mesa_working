*deck better.f
c***begin prologue     better
c***date written       950704   (yymmdd)
c***revision date               (yymmdd)
c***keywords           finite difference, band, eigenvalues, numerov
c***author             schneider, barry(nsf)
c***source             
c***purpose            calculate eigenvalues/eigenvectors of numerov
c***                   tridiagonal matrix from guesses using finite 
c***                   difference.
c***
c***description        the finite difference eigenvalues are used to start
c***                   off a root search procedure based on the numerov
c***                   algorithm.  the latter procedure gives much better
c***                   eigenvalues and eigenvectors for the same numerical
c***                   effort.  however you must be close enough using the
c***                   finite difference to bracket the roots.  crude
c***                   finite difference guesses can lead to no convergence
c***                   for the higher roots.
c***references
c
c***routines called
c***end prologue      better
      subroutine better(band,r,v,eig,vec,s,wts,ipiv,work,stp,nroots,n,
     1                  order,bw,maxit,bound,err,bc,type,prnt)
      implicit integer (a-z)
      dimension band(n,-bw:bw), r(n), v(n), eig(n), vec(n,*) ,work(n,3)
      dimension ipiv(n,2), wts(n), s(*)
      real*8 band, r, v, eig, vec, work, stp,  bound, err, s
      real*8 lower, upper, e, eps, val, test, rcond, prd, eold
      real*8 check, wts
      character*(*) type, bc
      logical prnt
      common /io/ inp, iout
      type='standard-finite-difference'
      write(iout,1) type, bc
      do 10 i=2,nroots+1
c        establish upper and lower bounds for search and number of
c        values to locate zero of determinant.
         lower=eig(i)-eig(i)*bound
         upper=eig(i)+eig(i)*bound
         eps=abs((upper-lower)/20)
         if(upper.lt.lower) then
            eps=-eps
         endif
         write(iout,2) i, lower, upper        
         e=lower
         call condit(e,rcond,band,r,v,ipiv,work,stp,n,order,bw,
     1               bc,type)
         val=sign(1.d0,rcond)
         prd=1.d0
         if(eps.gt.0.d0) then
            do while (e.ge.lower.and.e.le.upper.and.prd.gt.0.d0)
               eold=e
               e=e+eps
               call condit(e,rcond,band,r,v,ipiv,work,stp,n,order,bw,
     1                     bc,type)
               test=sign(1.d0,rcond)
               prd=test*val
            enddo
         else
            do while (e.le.lower.and.e.ge.upper.and.prd.gt.0.d0)
               eold=e
               e=e+eps
               call condit(e,rcond,band,r,v,ipiv,work,stp,n,order,bw,
     1                     bc,type)
               test=sign(1.d0,rcond)
               prd=test*val
            enddo
         endif
         if(prd.lt.0.d0) then
c           ok.  we found a sign change now zero in on the eigenvalue,
            lower=eold
            upper=e
            call fzero(lower,upper,lower,err,err,maxit,iflag,band,
     1                 r,v,ipiv,work,stp,n,order,bw,bc,type)
            eig(i)=lower
c           find corresponding eigenvector
 99         do 20 j=1,n
               work(j,1)=v(j)-lower
 20         continue   
            if (type.eq.'standard-finite-difference'
     1                     .or.
     2          type.eq.'standard-numerov') then
                call fdiff(band,work(1,1),stp,n,order,bw,type)
                call bndcnd(band,n,order,bw,bc)
            else
                call gfdiff(band,r,work(1,1),n,order,bw,type)
            endif
            vec(1,i)=0.d0
            vec(2,i)=1.d0
            vec(3,i)=-band(2,0)/band(2,1)
            do 30 j=3,n-2
               vec(j+1,i) = - ( band(j,0)*vec(j,i) + 
     1                          band(j,-1)*vec(j-1,i) )/band(j,1)
 30         continue
            check=band(n-1,-1)*vec(n-2,i)+band(n-1,0)*vec(n-1,i)
            write(iout,*) 'check = ',check   
            if(bc.eq.'zero-function') then
               vec(n,i)=0.d0
            else
               vec(n,i) = ( 4.d0*vec(n-1,i)-vec(n-2,i) )/3.d0
c                vec(n,i)=vec(n-1,i)
            endif
c            call sgtrmv(band(2,0),band(2,-1),band(2,1),vec(2,i),
c     1                  work(2,1),n-2)
c            write(iout,*) (work(j,1),j=2,n-1)
         else
c            no sign change found.  quit and proceed to next root.             
             write(iout,3) i
         endif
 10   continue   
      if(prnt) then
        call orth(s,vec(1,2),wts,r,work,nroots,n)
      endif
      return
 1    format('entering correction subroutine',/,5x,
     1       'type correction = ',a32/,5x,'boundary-condition = ',a32)
 2    format(/,'correction. fitting root = ',i3,/,
     1         'lower bound = ',e15.8,2x,'upper bound = ',e15.8)
 3    format(/,'no sign change in this interval for root = ',i3)
 4    format(/,'solution for root = ',i3,(/,5e15.8))
 5    format('orthonormality of vector = ',i3,1x,' vector = ',i3,
     1       ' is ',e15.8)
      end

