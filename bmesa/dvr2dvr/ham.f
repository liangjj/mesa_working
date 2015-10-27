*deck ham.f
c***begin prologue     ham
c***date written       960723   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c***keywords           simle one-dimensional hamiltonian
c***author             schneider, barry (nsf)
c***source             
c***purpose            matrix elements of hamiltonian.
c***                   
c***references         
c
c***routines called    
c***end prologue       ham
      subroutine ham(p,dp,ddp,wt,eig,hamil,n,npts,pottyp,grid,prn)
      implicit integer (a-z)
      real*8 p, dp, ddp, wt, eig, hamil
      logical prn
      character*80 title
      character*2 itoc
      character*(*) pottyp
      dimension p(npts,0:n-1), dp(npts,0:n-1), ddp(npts,0:n-1)
      dimension wt(npts), eig(n), hamil(n,n) 
      common/io/inp, iout
      call rzero(hamil,n*n)
      do 10 i=1,n
         do 20 j=1,i
            do 30 k=1,npts
               hamil(i,j) = hamil(i,j) - wt(k)*p(k,i-1)*ddp(k,j-1)
 30         continue
            hamil(i,j) = hamil(i,j) + p(npts,i-1)*dp(npts,j-1) 
     1                      - p(1,i-1)*dp(1,j-1)
            hamil(j,i) = hamil(i,j)
 20      continue   
 10   continue
      call vscale(hamil,hamil,.5d0,n*n)
      if(pottyp.eq.'harmonic-oscillator') then
         do 40 i=1,n
            hamil(i,i) = hamil(i,i) + .5d0*eig(i)*eig(i)
 40      continue
      elseif(pottyp.eq.'exponential') then
         do 50 i=1,n
            hamil(i,i) = hamil(i,i) - exp(-eig(i))
 50      continue
      elseif(pottyp.eq.'well') then
         do 60 i=1,n
            hamil(i,i) = hamil(i,i) - 1.0d0
 60      continue
      endif
      call iosys('write real "hamiltonian for grid '//itoc(grid)
     1           //'" to ham',n*n,hamil,0,' ')     
      if (prn) then
          title='hamiltonian'
          call prntrm(title,hamil,n,n,n,n,iout)
      endif     
      return
      end       



