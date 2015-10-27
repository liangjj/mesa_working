*deck fobj.f
c***begin prologue     fobj
c***date written       930623   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c***keywords           energy grid
c***author             schneider, barry (nsf)
c***source             
c***purpose            make the objective function for a newton-raphson
c***                   
c
c***references         
c
c***routines called    
c***end prologue      fobj
      subroutine fobj(f,df,rmat,drmat,kappa,n,der)
c
      implicit integer (a-z)
      real*8  f, df, rmat, drmat, kappa
      logical der 
      common/io/inp, iout
c
      f = 1.d0 + kappa * rmat
      if(der) then
         df = kappa * drmat - rmat / kappa       
      endif
      return
      end
