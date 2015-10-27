*deck h1tonv.f
c***begin prologue     h1tonv
c***date written       960723   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c***keywords           time development
c***author             schneider, barry (nsf)
c***source             
c***purpose            hamiltonian times space*time vector
c***                   
c***references         
c
c***routines called    
c***end prologue       h1tonv
      subroutine h1tonv(h1,ht,vecout,vecin,n1,nt,nc,nvc)
      implicit integer (a-z)
      real*8 h1, ht, vecout, vecin
      dimension h1(n1,n1), ht(nt,nt), vecin(n1,nt,nc,2,nvc)
      dimension vecout(n1,nt,nc,2,nvc)
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

