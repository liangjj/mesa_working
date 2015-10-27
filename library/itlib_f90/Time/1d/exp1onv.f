*deck exp1onv.f
c***begin prologue     exp1onv
c***date written       011013   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c***keywords           time development
c***author             schneider, barry (nsf)
c***source             
c***purpose            exponential of potential on vector
c***                   
c***references         
c
c***routines called    
c***end prologue       exp1onv
      subroutine exp1onv(vecin,vecout,v,t,n,nc,nvc)
      implicit integer (a-z)
      real*8 v, t
      complex*16 vecin, vecout
      dimension vecin(n,nc,nvc), vecout(n,nc,nvc), v(n,nc,nc)
      common/io/inp, iout
      call apbc(vecout,h1,vecin,n1,n1,nt*nc*2*nvc)
      do 10 ic=1,nc
         do 20 i=1,nvc
            call apbct(vecout(1,1,ic,1,i),vecin(1,1,ic,2,i),ht,n1,nt,nt)   
            call ambct(vecout(1,1,ic,2,i),vecin(1,1,ic,1,i),ht,n1,nt,nt)   
 20      continue   
 10   continue   
      call vneg(vecout,vecout,n1*nt*nc*2*nvc)
      return
      end       

