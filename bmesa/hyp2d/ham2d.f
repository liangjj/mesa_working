*deck ham2d.f
c***begin prologue     ham2d
c***date written       000710   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c***keywords           two-dim
c***author             schneider, barry (nsf)
c***source             3-dim
c***purpose            non zero hamiltonian elements and indices for 
c***                   two dimensional packed hamiltonian. 
c***                   
c***                   
c***references         
c
c***routines called    
c***end prologue       ham2d
      subroutine ham2d(ham,v,ind,h1,h2,q1,dscale,n)
      implicit integer (a-z)
      real*8 ham, v, h1, h2, q1, pre, dscale
      dimension n(3)
      dimension ham(n(3),n(3)), ind(n(2),n(1)), v(n(3))
      dimension h1(n(1),n(1)), h2(n(2),n(2)), q1(n(1))
      common/io/inp, iout
c
c 
      do 10 i=1,n(1)
         pre = dscale/(q1(i)*q1(i))
         do 20 j=1,n(2) 
            do 30 k=1,j-1
               if(h2(j,k).ne.0.d0) then
                  itot=ind(j,i)
                  jtot=ind(k,i)
                  ham(itot,jtot) = ham(itot,jtot) + pre*h2(j,k)
               endif
 30         continue
 20      continue
 10   continue      
      do 40 i=1,n(2)
         do 50 j=1,n(1)
            do 60 k=1,j-1
               if(h1(j,k).ne.0.d0) then                
                  itot=ind(i,j)
                  jtot=ind(i,k)
                  ham(itot,jtot) = ham(itot,jtot) + h1(i,j)  
               endif
 60         continue   
 50      continue   
 40   continue
c
c     do the diagonals
c
      cnt=0
      do 100 i=1,n(1)
         pre = .25d0*dscale/(q1(i)*q1(i))
         do 200 j=1,n(2)
            cnt=cnt+1
            ham(cnt,cnt) = ham(cnt,cnt) + h1(i,i) + h2(j,j) 
     1                                   + v(cnt) + pre
 200     continue
 100  continue   
      return
      end       










