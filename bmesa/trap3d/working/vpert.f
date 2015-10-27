*deck vpert.f
c***begin prologue     vpert
c***date written       970797   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c***keywords           3-d schroedinger equation
c***author             schneider, barry (nsf)
c***source             trap3d
c***purpose            interaction potential in dvr representation
c***                   
c***                   
c***                   
c***references         
c
c***routines called    
c***end prologue       vpert
      subroutine vpert(v,v1,v2,v3,n,n1,n2,n3,numdim)
      implicit integer (a-z)
      real*8 v, v1, v2, v3
      dimension v1(n1), v2(n2), v3(n3), v(n)
      if(numdim.eq.1) then
         do 10 i=1,n
            v(i) = v(i) - v1(i)
 10      continue
      elseif(numdim.eq.2) then
         count=0
         do 20 i=1,n1
            do 30 j=1,n2
               count=count+1
               v(count) = v(count) - v1(i) - v2(j)
 30         continue
 20      continue
      elseif(numdim.eq.3) then
         count=0
         do 40 i=1,n1
            do 50 j=1,n2
               do 60 k=1,n3
                  count=count+1
                  v(count) = v(count) - v1(i) - v2(j) - v3(k)
 60            continue
 50         continue   
 40      continue
      endif   
      return
      end       
