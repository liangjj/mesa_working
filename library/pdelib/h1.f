*deck h1.f
c***begin prologue     h1
c***date written       960723   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c***keywords           time development
c***author             schneider, barry (nsf)
c***source             
c***purpose            find time-dependent wavefunction by expansion.
c***                   
c***description        simple time-dependent hamiltonian using polynomial
c***                   basis.  
c         
c                           d
c                      hbar _    
c                           dt
c***references         
c
c***routines called    
c***end prologue       h1
      subroutine h1(p,dp,wt,hamt,hbar,n,prn)
      implicit integer (a-z)
      real*8 p, dp, wt, hbar, hamt
      logical prn
      character*80 title
      dimension p(n,n), dp(n,n)
      dimension wt(n), hamt(n,n)
      common/io/inp, iout
      call rzero(hamt,n*n)
      do 10 i=1,n
         do 20 j=1,n
            hamt(i,j) = wt(i)*p(i,i)*dp(i,j)
 20      continue   
 10   continue
      call sscal(n*n,hbar,hamt,1)
      if (prn) then
          title='time-dependent hamiltonian'
          call prntrm(title,hamt,n,n,n,n,iout)
      endif     
      return
      end       









