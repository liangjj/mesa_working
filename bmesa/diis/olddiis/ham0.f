*deck ham0.f
c***begin prologue     ham0
c***date written       960723   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c***keywords           harmonic oscillator
c***author             schneider, barry (nsf)
c***source             
c***purpose            harmonic oscillator eigenvalues and vectors.
c***                   
c***description        construct hamiltonian in polynomial basis.  
c***references         
c
c***routines called    
c***end prologue       ham0
      subroutine ham0(p,dp,ddp,wt,eigr,ham,omega,n,npts,pottyp,
     1                tfermi,prn)
      implicit integer (a-z)
      real*8 p, dp, ddp, wt, eigr, ham, omega, fac
      logical tfermi, prn
      character*80 title
      character*(*) pottyp
      dimension p(npts,0:n-1), dp(npts,0:n-1), ddp(npts,0:n-1)
      dimension wt(npts), eigr(n), ham(n,n) 
      common/io/inp, iout
      call rzero(ham,n*n)
      if(.not.tfermi) then
         do 10 i=1,n
            do 20 j=1,i
               do 30 k=1,npts
                  ham(i,j) = ham(i,j) - wt(k)*p(k,i-1)*ddp(k,j-1)
 30            continue
               ham(i,j) = ham(i,j) + p(npts,i-1)*dp(npts,j-1) 
     1                             - p(1,i-1)*dp(1,j-1)
               ham(j,i) = ham(i,j)
 20         continue   
 10      continue
      endif
      if(pottyp.eq.'harmonic-oscillator') then
         fac = .25d0*omega*omega
c         fac=.25d0
         do 40 i=1,n
            ham(i,i) = ham(i,i) + fac*eigr(i)*eigr(i)
 40      continue
      elseif(pottyp.eq.'exponential') then
         do 50 i=1,n
            ham(i,i) = ham(i,i) - 2.0d0*exp(-eigr(i))
 50      continue
      elseif(pottyp.eq.'well') then
         do 60 i=1,n
            ham(i,i) = ham(i,i) - 2.0d0
 60      continue
      endif
      if (prn) then
          title='hamiltonian'
          call prntrm(title,ham,n,n,n,n,iout)
      endif     
      return
      end       

