*deck ham0.f
c***begin prologue     ham0
c***date written       960723   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c***keywords           harmonic oscillator
c***author             schneider, barry (nsf)
c***source             
c***purpose            harmonic oscillator eigenvalues and vectors.
c***                   
c***description        diagonalize radial harmonic oscillator in polynomial
c***                   basis.  
c***references         
c
c***routines called    
c***end prologue       ham0
      subroutine ham0(p,dp,ddp,wt,eigr,hamil,v,n,npts,pottyp,prn)
      implicit integer (a-z)
      real*8 p, dp, ddp, wt, eigr, hamil, v
      logical prn
      character*80 title
      character*(*) pottyp
      dimension p(npts,0:n-1), dp(npts,0:n-1), ddp(npts,0:n-1)
      dimension wt(npts), eigr(n), hamil(n,n), v(n,n) 
      common/io/inp, iout 
      call iosys('write real "dvr space eigenvalues" to lamdat',
     1            n,eigr,0,' ')
      call iosys('write real "dvr space functions" to lamdat',
     1            n*npts,p,0,' ')
      call iosys('write real "first derivative of dvr space '//
     1           'functions" to lamdat',n*npts,dp,0,' ')
      call iosys('write real "second derivative of dvr space '//
     1           'functions" to lamdat',n*npts,ddp,0,' ')
      call rzero(v,n*n)
      do 10 i=1,n
         do 20 j=1,i
            do 30 k=1,npts
               v(i,j) = v(i,j) - wt(k)*p(k,i-1)*ddp(k,j-1)
 30         continue
            v(i,j) = v(i,j) + p(npts,i-1)*dp(npts,j-1) 
     1                      - p(1,i-1)*dp(1,j-1)
            v(j,i) = v(i,j)
 20      continue   
 10   continue
      if(pottyp.eq.'harmonic-oscillator') then
         do 40 i=1,n
            v(i,i) = v(i,i) + .25d0*eigr(i)*eigr(i)
 40      continue
      elseif(pottyp.eq.'exponential') then
         do 50 i=1,n
            v(i,i) = v(i,i) - 2.0d0*exp(-eigr(i))
 50      continue
      elseif(pottyp.eq.'well') then
         do 60 i=1,n
            v(i,i) = v(i,i) - 2.0d0
 60      continue
      endif
      call copy(v,hamil,n*n)   
      if (prn) then
          title='unperturbed hamiltonian'
          call prntrm(title,hamil,n,n,n,n,iout)
      endif     
      return
      end       
