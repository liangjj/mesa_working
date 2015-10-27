*deck ke2exp.f 
c***begin prologue     ke2exp
c***date written       960718   (yymmdd)
c***revision date               (yymmdd)
c***keywords           time, dvr, split-operator
c***                   
c***author             schneider, b. i.(nsf)
c***source             ke2exp
c***purpose            split operator propagation.
c***                    
c***references       
c
c***routines called    iosys, util and mdutil
c***end prologue       ke2exp

      subroutine ke2exp(psi,eig1,eig2,u1,u2,t,work,n1,n2)
      implicit integer (a-z)
      real*8 eig1, eig2, u1, u2, t
      complex*16 psi
      complex*16 temp, work, eye
      dimension psi(n2,n1)
      dimension eig1(n1), eig2(n2)
      dimension u1(n1,n1), u2(n2,n2)
      dimension work(n2,n1)
      common/io/inp, iout     
      data eye / (0.d0,1.d0) /
c
c     transform to the diagonal representation of KE1.
c
      call ecbc(work,psi,u1,n2,n1,n1)
c
c        multiply by the diagonal KE1 in the transformed representation
c
      do 10 i=1,n1
         temp = exp(-eye*eig1(i)*t)
         do 20 j=1,n2
            work(j,i) = temp*work(j,i) 
 20      continue
 10   continue   
c
c        transform back to the original representation
c
      call ecbct(psi,work,u1,n2,n1,n1)
c
c     do the same for KE2
c
      call ebtcc(work,u2,psi,n2,n2,n1)
      do 30 j=1,n2
         temp = exp(-eye*eig2(j)*t)
         do 40 i=1,n1
            work(j,i) = temp*work(j,i)
 40      continue
 30   continue   
      call ebcc(psi,u2,work,n2,n2,n1)
      return
      end





















