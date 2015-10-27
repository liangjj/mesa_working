*deck vnlse.f
c***begin prologue     vnlse
c***date written       960723   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c***keywords           copy
c***author             schneider, barry (nsf)
c***source             
c***purpose            mean time dependent hartree field
c***                   
c***                   
c***references         
c
c***routines called    
c***end prologue       vnlse
      subroutine vnlse(v,psi,gamma,p1,p2,p3,ni,nj,n,nd,nt,nc,dim)
      implicit integer (a-z)
      real*8 v, psi, gamma, p1, p2, p3
      real*8 prdj, prdk, prdl
      dimension nd(dim)
      dimension v(n,nt,nc,nc), psi(n,nt,nc,2)
      dimension p1(nd(1),nd(1)), p2(nd(2),nd(2)), q3(nd(3),nd(3))
      common/io/inp, iout
      if(dim.eq.1) then
         do 10 i=1,nt
            do 20 j=1,nd(1)
               v(j,i,ni,nj) = v(j,i,ni,nj) + 
     1                  gamma * ( psi(j,i,ni,1)*psi(j,i,ni,1) 
     2                          + psi(j,i,ni,2)*psi(j,i,ni,2)
     3                          + psi(j,i,nj,1)*psi(j,i,nj,1)
     4                          + psi(j,i,nj,2)*psi(j,i,nj,2) )
     3                          * p1(j,j)*p1(j,j)
 20         continue   
 10      continue
      elseif(dim.eq.2) then
         do 30 i=1,nt
            cnt=0
            do 40 j=1,nd(1)
               prdj = p1(j,j)*p1(j,j)
               do 50 k=1,nd(2)
                  prdk = p2(k,k)*p2(k,k)
                  cnt=cnt+1 
                  v(cnt,i,ni,nj) = v(cnt,i,ni,nj) + gamma * 
     1                       ( psi(cnt,i,ni,1)*psi(cnt,i,ni,1)
     2                      +  psi(cnt,i,ni,2)*psi(cnt,i,ni,2)
     3                      +  psi(cnt,i,nj,1)*psi(cnt,i,nj,1)
     4                      +  psi(cnt,i,nj,2)*psi(cnt,i,nj,2) )
     5                      *  prdj * prdk
 50            continue
 40         continue   
 30      continue   
      elseif(dim.eq.3) then
         do 60 i=1,nt
            cnt=0 
            do 70 j=1,nd(1)
               prdj = p1(j,j)*p1(j,j)
               do 80 k=1,nd(2)
                  prdk = p2(k,k)*p2(k,k)
                  do 90 l=1,nd(3)
                     prdl = p3(l,l)*p3(l,l)
                     cnt = cnt + 1
                     v(cnt,i,ni,nj) = v(cnt,i,ni,nj) + gamma    * 
     1                          ( psi(cnt,i,ni,1)*psi(cnt,i,ni,1)
     2                          + psi(cnt,i,ni,2)*psi(cnt,i,ni,2)
     3                          + psi(cnt,i,nj,1)*psi(cnt,i,nj,1)
     4                          + psi(cnt,i,nj,2)*psi(cnt,i,nj,2) )
     3                          * prdj * prdk * prdl
 90               continue
 80            continue
 70         continue
 60      continue
      else
         call lnkerr('dimension error')  
      endif
      return
      end       


