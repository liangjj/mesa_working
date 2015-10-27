*deck ke1exp.f 
c***begin prologue     ke1exp
c***date written       960718   (yymmdd)
c***revision date               (yymmdd)
c***keywords           time, dvr, split-operator
c***                   
c***author             schneider, b. i.(nsf)
c***source             ke1exp
c***purpose            split operator propagation.
c***                    
c***references       
c
c***routines called    iosys, util and mdutil
c***end prologue       ke1exp

      subroutine ke1exp(psi,eig1,u1,t,work,n)
      implicit integer (a-z)
      real*8 eig1, u1, t
      complex*16 psi
      complex*16 work, eye
      dimension psi(n)
      dimension eig1(n)
      dimension u1(n,n)
      dimension work(*)
      common/io/inp, iout     
      data eye / (0.d0,1.d0) /
c
c     transform to the diagonal representation of KE.
c
      call ebtcc(work,u1,psi,n,n,1)
c
c        multiply by the diagonal KE in the transformed representation
c
      do 10 i=1,n
         work(i) = work(i)*exp(-eye*eig1(i)*t)
 10   continue   
c
c        transform back to the original representation
c
      call ebcc(psi,u1,work,n,n,1)
      return
      end





















