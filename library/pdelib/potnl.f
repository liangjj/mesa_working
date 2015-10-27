*deck potnl.f
c***begin prologue     potnl
c***date written       970813   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c***keywords           non-linear schroedinger equation, gross-pitaevski
c***author             schneider, barry (nsf)
c***source             trap3d
c***purpose            form current approximation to non-linear
c***                   potential.
c***                   
c***                   
c***references         
c
c***routines called    
c***end prologue       potnl
      subroutine potnl(vnl,eig1,eig2,eig3,p1,p2,p3,psi,gamma,
     1                 n,n1,n2,n3,nd1,nd2,nd3,dim)
      implicit integer (a-z)
      real*8 vnl, eig1, eig2, eig3, p1, p2, p3, psi
      real*8 gamma
      dimension vnl(n), eig1(n1), eig2(n2), eig3(n3)
      dimension p1(nd1,n1), p2(nd2,n2)
      dimension p3(nd3,n3), psi(n)
      common/io/inp, iout
      if(dim.eq.1) then
         do 10 i=1,n
            vnl(i)=p1(i,i)*psi(i)/eig1(i)
            vnl(i)=vnl(i)*vnl(i)
 10      continue
      elseif(dim.eq.2) then
         count=0
         do 20 i=1,n1
            do 30 j=1,n2
               count=count+1
               vnl(count) = p1(i,i)*p2(j,j)*psi(count)
               vnl(count)=vnl(count)*vnl(count)
 30         continue   
 20      continue   
      elseif(dim.eq.3) then
         count=0
         do 40 i=1,n1
            do 50 j=1,n2
               do 60 k=1,n3
                  count=count+1
                  vnl(count) = p1(i,i)*p2(j,j)*p3(k,k)*psi(count)
                  vnl(count)=vnl(count)*vnl(count)
 60            continue   
 50         continue   
 40      continue   
      endif
      call smul(vnl,vnl,gamma,n) 
      return 
      end       





