*deck cmpare.f
c***begin prologue     cmpare
c***date written       960723   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c***keywords           compare approximate and exact model wavefunction
c***author             schneider, barry (nsf)
c***source             
c***purpose            
c***                   
c***references         
c
c***routines called    
c***end prologue       cmpare
      subroutine cmpare(coef,pt,t,t0,psi0,nt,type)
      implicit integer (a-z)
      real*8 pt, t, t0, coef
      complex*16 exact, approx, eye, psi0
      character*(*) type
      dimension pt(nt,nt), t(nt)
      dimension coef(nt,2)
      common/io/inp, iout
      data eye/(0.d0,1.d0)/
      if(type.ne.'t') then
         call lnkerr('error in time perturbation')
      else 
         write(iout,1)
         do 10 i=1,nt
            exact = ( exp(-eye*(t(i)*t(i) - t0*t0)*.5d0) )*psi0
            approx = ( coef(i,1) + eye*coef(i,2) )*pt(i,i) + psi0
            write(iout,2) t(i), approx, exact
 10      continue
      endif
  1   format(/,1x,'comparison between approximate and exact solution ',
     1       /,1x,'for linear potential in time')
  2   format(/,1x,'time        = ',e15.8,/,1x,
     1            'approximate = ',e15.8,1x,e15.8,/,1x, 
     2            'exact       = ',e15.8,1x,e15.8)                   
      return
      end       
