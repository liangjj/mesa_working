*deck two
      subroutine two(vcij,ec,energy,rbox)
      implicit integer (a-z)
      common /io/ inp, iout
      real*8 vcij, ec, a, b, c, root, vec, sum, dum
      real*8 energy, rbox, ksq, k, chi, psi
      dimension vcij(3), root(2), vec(2,2), dum(2), chi(2)
      dimension ec(2), psi(2)
      dum(1)=vcij(1)+ec(1)
      dum(2)=vcij(3)+ec(2)
      b=-(dum(1)+dum(2))
      c=dum(1)*dum(2)-vcij(2)*vcij(2)
      root(1)=-b+sqrt(b*b-4.d0*c)
      root(2)=-b-sqrt(b*b-4.d0*c)
      root(1)=.5d0*root(1)
      root(2)=.5d0*root(2)
      do 10 i=1,2
         a=dum(1)-root(i)
         vec(1,i)=1.d0
         vec(2,i)=-a/vcij(2)
         sum=vec(1,i)*vec(1,i)+vec(2,i)*vec(2,i)
         sum=1.d0/sqrt(sum)
         vec(1,i)=sum*vec(1,i)
         vec(2,i)=sum*vec(2,i)
         write(iout,1) vec(1,i), vec(2,i)
 10   continue
      do 20 i=1,2
         ksq=energy-root(i)
         if (ksq.ge.0.d0) then
             k=sqrt(2.d0*ksq)
             chi(i)=sin(k*rbox)
         else
             ksq=root(i)-energy
             k=sqrt(2.d0*ksq)
             chi(i)=sinh(k*rbox)
         endif
 20   continue
      psi(1)=vec(1,1)*chi(1)+vec(1,2)*chi(2)
      psi(2)=vec(2,1)*chi(1)+vec(2,2)*chi(2)
 1    format(/,1x,'vector',1x,2e15.8) 
      return
      end



