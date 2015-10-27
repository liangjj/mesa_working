*deck errmat.f
c***begin prologue     errmat
c***date written       951229   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c***keywords           non-linear schroedinger equation, gross-pitaevski
c***author             schneider, barry (nsf)
c***source             tstdiis
c***purpose            form current approximation to the diis error vector
c***                   
c***                   
c***references         
c
c***routines called    
c***end prologue       errmat
      subroutine errmat(f,evec,rho,error,n,iter,prnt)
      implicit integer (a-z)
      real*8 f, evec, rho
      real*8 error
      character*3 itoc
      character*80 title
      logical prnt
      dimension f(n,n), rho(n,n), evec(n,n)
      common/io/inp, iout
c
c     calculate error matrix = f*rho - rho*f = [f,rho]
c
      call ebc(evec,f,rho,n,n,n)
      call ambc(evec,rho,f,n,n,n)
      error=0.d0
      do 10 i=1,n
         do 20 j=1,i
            error = max(error,abs(evec(i,j)))
 20      continue
 10   continue   
      if(prnt) then
         write(iout,1) iter
         title='error matrix iteration = '//itoc(iter)
         call prntfm(title,evec,n,n,n,n,iout)
      endif                  
      return 
   1  format(/,5x,'error matrix formation iteration =', i3)
      end       
