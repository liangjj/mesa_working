*deck ham1.f
c***begin prologue     ham1
c***date written       960723   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c***keywords           time development
c***author             schneider, barry (nsf)
c***source             
c***purpose            find time-dependent wavefunction by expansion.
c***                   
c***description        simple time-dependent hamiltonian using polynomial
c***                   basis.  
c***references         
c
c***routines called    
c***end prologue       ham1
      subroutine ham1(p,dp,ddp,wt,eigt,hamil,n,npts,type,prn)
      implicit integer (a-z)
      real*8 p, dp, ddp, wt, eigt
      complex*16 hamil, eye
      logical prn
      character*80 title
      character*(*) type
      dimension p(npts,0:n-1), dp(npts,0:n-1), ddp(npts,0:n-1)
      dimension wt(npts), eigt(n), hamil(n,n)
      common/io/inp, iout
      data eye/(0.d0,1.d0)/ 
      call iosys('write real "dvr time eigenvalues" to lamdat',
     1            n,eigt,0,' ')
      call iosys('write real "dvr time functions" to lamdat',
     1            n*npts,p,0,' ')
      call iosys('write real "first derivative of dvr time '//
     1           'functions" to lamdat',n*npts,dp,0,' ')
      call iosys('write real "second derivative of dvr time '//
     1           'functions" to lamdat',n*npts,ddp,0,' ')
      call czero(hamil,n*n)
      do 10 i=1,n
         do 20 j=1,n
            do 30 k=1,npts
               hamil(i,j) = hamil(i,j) + wt(k)*p(k,i-1)*dp(k,j-1)
 30         continue
c            hamil(j,i) = - hamil(i,j)
             hamil(i,j) = eye*hamil(i,j) 
c            hamil(j,i) = eye*hamil(j,i)  
 20      continue   
 10   continue
      if(type.eq.'t') then
         do 40 i=1,n
            hamil(i,i) = hamil(i,i) - eigt(i)
 40      continue
      elseif(type.eq.'1') then
         do 50 i=1,n
            hamil(i,i) = hamil(i,i) + 1.d0
 50      continue
      endif
      if (prn) then
          title='time-dependent hamiltonian'
          call prntcm(title,hamil,n,n,n,n,iout)
      endif     
      return
      end       
