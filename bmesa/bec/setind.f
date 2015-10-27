*deck setind.f
c***begin prologue     setind
c***date written       951229   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c***keywords           two-dim
c***author             schneider, barry (nsf)
c***source             3-dim
c***purpose            form two and three dimensional index arrays.
c***                   
c***                   
c***references         
c
c***routines called    
c***end prologue       setind
      subroutine setind(ind,nply,n,dim)
      implicit integer (a-z)
      dimension ind(n,*), nply(dim)
      common/io/inp, iout 
      if(dim.eq.2) then
         cnt=0 
         do 10 i=1,nply(1)
            ii = nply(2)*(i-1)
            do 20 j=1,nply(2)
               cnt=cnt+1
               ind(cnt,1)=i
               ind(cnt,2)=j
               ind(cnt,3) = ii + j
 20         continue
 10      continue   
      elseif(dim.eq.3) then
         nsq=nply(2)*nply(3)
         cnt=0
         do 30 i=1,nply(1)
            ii = nsq*(i-1)
            do 40 j=1,nply(2)
               jj = nply(3)*(j-1)
               do 50 k=1,nply(3)
                  cnt=cnt+1
                  ind(cnt,1)=i
                  ind(cnt,2)=j
                  ind(cnt,3)=k
                  ind(cnt,4) = ii + jj + k
 50            continue
 40         continue
 30      continue   
      endif
      return
      end       

