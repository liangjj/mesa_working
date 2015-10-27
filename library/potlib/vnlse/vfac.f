*deck vfac.f
c***begin prologue     vfac
c***date written       960723   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c***keywords           copy
c***author             schneider, barry (nsf)
c***source             
c***purpose            gaussian time factors
c***                   
c***                   
c***references         
c
c***routines called    
c***end prologue       vfac
      subroutine vfac(v,fac,q1,q2,q3,ni,nj,n,nd,nt,nc,dim)
      implicit integer (a-z)
      real*8 v, fac, q1, q2, q3
      dimension nd(dim)
      dimension v(n,nt,nc,nc), fac(nt), q1(nd(1)), q2(nd(2)), q3(nd(3))
      common/io/inp, iout
      if(dim.eq.1) then
         do 10 i=1,nt
            do 20 j=1,nd(1)
               v(j,i,ni,nj) = v(j,i,ni,nj) + fac(i) * q1(j)
 20         continue   
 10      continue
      elseif(dim.eq.2) then
         do 30 i=1,nt
            cnt=0
            do 40 j=1,nd(1)
               do 50 k=1,nd(2)
                  cnt=cnt+1 
                  v(cnt,i,ni,nj) = v(cnt,i,ni,nj) + fac(i) * 
     1                                        sqrt ( q1(j) * q1(j)
     2                                                     + 
     3                                               q2(k) * q2(k) )
 50            continue
 40         continue   
 30      continue   
      elseif(dim.eq.3) then
         do 60 i=1,nt
            cnt=0 
            do 70 j=1,nd(1)
               do 80 k=1,nd(2)
                  do 90 l=1,nd(3)
                     cnt = cnt + 1
                     v(cnt,i,ni,nj) = v(cnt,i,ni,nj) 
     1                                + fac(i) * sqrt( q1(j)*q1(j) + 
     2                                                 q2(k)*q2(k) + 
     3                                                 q3(l)*q3(l) )
 90               continue
 80            continue
 70         continue
 60      continue
      else
         call lnkerr('dimension error')  
      endif
      return
      end       


