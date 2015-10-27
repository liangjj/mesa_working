*deck h3tonv.f
c***begin prologue     h3tonv
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
c***end prologue       h3tonv
      subroutine h3tonv(h1,h2,h3,ht,vecout,vecin,n1,n2,n3,nt,nc,nvc)
      implicit integer (a-z)
      real*8 h1, h2, h3, ht, vecout, vecin
      dimension h1(n1,n1), h2(n2,n2), h3(n3,n3)
      dimension vecin(n3,n2,n1,nt,nc,2,nvc)
      dimension vecout(n3,n2,n1,nt,nc,2,nvc)
      common/io/inp, iout
      call apbc(vecout,h3,vecin,n3,n3,n2*n1*nt*nc*2*nvc)
      do 10 i=1,n1
         do 20 j=1,nt
            do 30 ic=1,nc
               do 40 k=1,nvc
                  call apbct(vecout(1,1,i,j,ic,1,k),
     1                       vecin (1,1,i,j,ic,1,k),h2,n3,n2,n2)
                  call apbct(vecout(1,1,i,j,ic,2,k),
     1                       vecin (1,1,i,j,ic,2,k),h2,n3,n2,n2)
 40            continue   
 30         continue   
 20      continue
 10   continue
      do 50 j=1,nt
         do 60 ic=1,nc
            do 70 k=1,nvc
               call apbct(vecout(1,1,1,j,ic,1,k),
     1                    vecin (1,1,1,j,ic,1,k),h1,n3*n2,n1,n1)
               call apbct(vecout(1,1,1,j,ic,2,k),
     1                    vecin (1,1,1,j,ic,2,k),h1,n3*n2,n1,n1)
 70         continue
 60      continue
 50   continue   
      do 100 ic=1,nc
         do 200 k=1,nvc
            call apbct(vecout(1,1,1,1,ic,1,k),
     1                 vecin (1,1,1,1,ic,2,k),ht,n3*n2*n1,nt,nt)
            call ambct(vecout(1,1,1,1,ic,2,k),
     1                 vecin (1,1,1,1,ic,1,k),ht,n3*n2*n1,nt,nt)
 200     continue   
 100  continue   
      call vneg(vecout,vecout,n3*n2*n1*nt*nc*2*nvc)
      return
      end       

