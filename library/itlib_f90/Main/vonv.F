*deck vonv.f
c***begin prologue     vonv
c***date written       960723   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c***keywords           time development
c***author             schneider, barry (nsf)
c***source             
c***routines called    
c***end prologue       vonv
      subroutine vonv(vecout,v,vecin,n,nc,m)
      implicit integer (a-z)
      real*8 vecout, v, vecin, one
      dimension vecout(n,nc,m), vecin(n,nc,m), v(n,nc,nc)
      common/io/inp, iout
      call rzero(vecout,n*nc*m)
      do 10 i=1,n
         do 20 j=1,m
            do 30 k=1,nc
               do 40 l=1,nc
                  vecout(i,k,j) = vecout(i,k,j) + v(i,k,l)*vecin(i,l,j) 
 40            continue
 30         continue
 20      continue
 10   continue   
      return
      end       
