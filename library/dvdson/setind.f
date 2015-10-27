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
      subroutine setind(ind,nply,n,ng,dim,ngrid)
      implicit integer (a-z)
      dimension nply(ngrid,3)
      dimension ind(n,*)
      common/io/inp, iout
      if(dim.eq.1) then
         cnt=0
         do 10 i=1,nply(ng,1)
            cnt=cnt+1
            ind(cnt,1)=i
            ind(cnt,2)=cnt
            ind(cnt,3)=i
            ind(cnt,4)=cnt
 10      continue   
      elseif(dim.eq.2) then
         cnt=0 
         do 20 i=1,nply(ng,1)
            do 30 j=1,nply(ng,2)
               cnt=cnt+1
               ind(cnt,1)=i
               ind(cnt,2)=j
               ind(cnt,3) = cnt
               ind(cnt,4)=i
               ind(cnt,5)=j
               ind(cnt,6) = cnt               
 30         continue
 20      continue
      elseif(dim.eq.3) then
         cnt=0
         do 40 i=1,nply(ng,1)
            do 50 j=1,nply(ng,2)
               do 60 k=1,nply(ng,3)
                  cnt=cnt+1
                  ind(cnt,1)=i
                  ind(cnt,2)=j
                  ind(cnt,3)=k
                  ind(cnt,4) = cnt
                  ind(cnt,5)=i
                  ind(cnt,6)=j
                  ind(cnt,7)=k
                  ind(cnt,8) = cnt                  
 60            continue
 50         continue
 40      continue   
      endif
      return
      end       

