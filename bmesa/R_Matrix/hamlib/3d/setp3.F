*deck setp3.f
c***begin prologue     setp3
c***date written       951229   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c***keywords           two-dim
c***author             schneider, barry (nsf)
c***source             3-dim
c***purpose            form three dimensional index array.
c***                   
c***                   
c***references         
c
c***routines called    
c***end prologue       setp3
      subroutine setp3(ind,n,prnt)
      implicit integer (a-z)
      logical prnt
      dimension n(3)
      dimension ind(n(3),n(2),n(1))
      common/io/inp, iout
      cnt=0
      do 10 i=1,n(1)
         do 20 j=1,n(2)
            do 30 k=1,n(3)
               cnt=cnt+1
               ind(k,j,i)=i
 30         continue
 20      continue
 10   continue 
      if(prnt) then
         write(iout,1)
         cnt=0
         do 40 i=1,n(1)
            do 50 j=1,n(2)
               do 60 k=1,n(3)
                  cnt=cnt+1
                  write(iout,2) k, j, i, ind(k,j,i)
 60            continue
 50         continue
 40      continue 
      endif
      return
 1    format(/,1x,'hamiltonian indices')
 2    format('k = ',i3,1x,'j = ',i3,1x,'i  = ',i3,1x,'index = ',i5)
      end       

