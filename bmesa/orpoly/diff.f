*deck diff.f
c***begin prologue     diff
c***date written       951229   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c***keywords           chebychev polynomials
c***author             schneider, barry (nsf)
c***source             math
c***purpose            solve differential equation by expanding in
c***                   orthogonal polynomials.
c***                   
c***references         
c
c***routines called    
c***end prologue       diff
      subroutine diff(pn,dpn,ddpn,x,wts,ham,eig,work,der,n,npts,
     1                pottyp,prnt)
      implicit integer (a-z)
      real*8 pn, dpn, ddpn, x, ham, eig, work, wts, der
      character*80 title
      character*(*) pottyp
      logical prnt
      dimension pn(npts,0:n-1), dpn(npts,0:n-1), ddpn(npts,0:n-1)
      dimension x(npts), ham(n,n), eig(n), work(n), wts(n), der(2)
      common/io/inp, iout 
      do 10 i=1,n
         do 20 j=1,i
            ham(i,j)=0.d0
            do 30 k=1,npts
               ham(i,j) = ham(i,j) - .5d0*wts(k)*pn(k,i-1)*ddpn(k,j-1)
 30         continue
            ham(i,j) = ham(i,j) +.5d0 *
     1           ( pn(npts,i-1)*( dpn(npts,j-1) -der(2)*pn(npts,j-1) )
     2                          - 
     3             pn(1,i-1)*( dpn(1,j-1) - der(1)*pn(1,j-1) ) )
            ham(j,i)=ham(i,j)
 20      continue   
 10   continue
      if (pottyp.eq.'one') then
          do 40 i=1,n
             do 50 j=1,i
                do 60 k=1,npts
                   ham(i,j) = ham(i,j) - pn(k,i-1)*wts(k)*exp(-x(k))*
     1                                   pn(k,j-1)
 60             continue
                ham(j,i) = ham(i,j)
 50          continue
 40       continue   
      elseif(pottyp.eq.'one') then
          do 70 i=1,n
             do 80 j=1,i
                do 90 k=1,npts
                   ham(i,j) = ham(i,j) + pn(k,i-1)*wts(k)*pn(k,j-1)
 90             continue
                ham(j,i) = ham(i,j)
 80          continue
 70       continue   
      elseif(pottyp.eq.'coulomb') then
          do 100 i=1,n
             do 110 j=1,i
                do 120 k=1,npts
                   ham(i,j) = ham(i,j) - pn(k,i-1)*wts(k)*pn(k,j-1)/x(k)
 120            continue
                ham(j,i) = ham(i,j)
 110         continue
 100      continue   
      endif
      if (prnt) then
          title='hamiltonian'   
          call prntrm(title,ham,n,n,n,n,iout)
      endif
      call tred2(n,n,ham,eig,work,ham)
      call tql2(n,n,eig,work,ham,ierr)
      title='eigenvalues'
      call prntrm(title,eig,n,1,n,1,iout)
      return
 1    format('x = ',e15.8,1x,'fexact = ',e15.8,1x,'fcheb = ',e15.8)
      end       

