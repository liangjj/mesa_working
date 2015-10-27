*deck h2tonv.f
c***begin prologue     h2tonv
c***date written       960723   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c***keywords           time development
c***author             schneider, barry (nsf)
c***source             
c***purpose            space hamiltonian times space*time vector
c***                   
c***references         
c
c***routines called    
c***end prologue       h2tonv
      subroutine h2tonv(h1,h2,ht,vecout,vecin,n1,n2,nt,nc,nvc)
      implicit integer (a-z)
      real*8 h1, h2, ht, vecout, vecin
      dimension h1(n1,n1), h2(n2,n2), ht(nt,nt)
      dimension vecin(n2,n1,nt,nc,2,nvc), vecout(n2,n1,nt,nc,2,nvc)
      common/io/inp, iout
      call apbc(vecout,h2,vecin,n2,n2,n1*nt*nc*2*nvc)
      do 10 i=1,nt
         do 20 ic=1,nc
            do 30 j=1,nvc
               call apbct(vecout(1,1,i,ic,1,j),vecin(1,1,i,ic,1,j),
     1                    h1,n2,n1,n1)
               call apbct(vecout(1,1,i,ic,2,j),vecin(1,1,i,ic,2,j),
     1                    h1,n2,n1,n1)
 30         continue   
 20      continue
 10   continue
      do 40 ic=1,nc
         do 50 j=1,nvc
            call apbct(vecout(1,1,1,ic,1,j),vecin(1,1,1,ic,2,j),
     1                 ht,n2*n1,nt,nt)
            call ambct(vecout(1,1,1,ic,2,j),vecin(1,1,1,ic,1,j),
     1                 ht,n2*n1,nt,nt)
 50      continue   
 40   continue
      call vneg(vecout,vecout,n2*n1*nt*nc*2*nvc)   
      return
      end       

