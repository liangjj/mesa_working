*deck ke3exp.f 
c***begin prologue     ke3exp
c***date written       960718   (yymmdd)
c***revision date               (yymmdd)
c***keywords           time, dvr, split-operator
c***                   
c***author             schneider, b. i.(nsf)
c***source             ke3exp
c***purpose            split operator propagation.
c***                    
c***references       
c
c***routines called    iosys, util and mdutil
c***end prologue       ke3exp

      subroutine ke3exp(psi,eig1,eig2,eig3,u1,u2,u3,t,
     1                  work,n1,n2,n3)
      implicit integer (a-z)
      real*8 eig1, eig2, eig3, u1, u2, u3, t
      complex*16 psi
      complex*16 temp, work, eye
      dimension psi(n3,n2,n1)
      dimension eig1(n1), eig2(n2), eig3(n3)
      dimension u1(n1,n1), u2(n2,n2), u3(n3,n3)
      dimension work(n3,n2,n1)
      common/io/inp, iout     
      data eye / (0.d0,1.d0) /
c
c     transform to the diagonal representation of KE1.
c
      call ecbc(work,psi,u1,n3*n2,n1,n1)
c
c        multiply by the diagonal KE1 in the transformed representation
c
      do 10 i=1,n1
         temp = exp(-eye*eig1(i)*t)
         do 20 j=1,n2
            do 30 k=1,n3
               work(k,j,i) = temp*work(k,j,i)
 30         continue
 20      continue
 10   continue   
c
c        transform back to the original representation
c
      call ecbct(psi,work,u1,n3*n2,n1,n1)
c
c     do the same for KE2
c
      do 40 i=1,n1
         call ecbc(work(1,1,i),psi(1,1,i),u2,n3,n2,n2)
 40   continue
      do 50 j=1,n2
         temp = exp(-eye*eig2(i)*t)
         do 60 i=1,n1
            do 70 k=1,n3
               work(k,j,i) = temp*work(k,j,i)
 70         continue
 60      continue
 50   continue   
      do 80 i=1,n1
         call ecbct(psi,work(1,1,i),u2,n3,n2,n2)
 80   continue   
c
c     and finally the same for KE3
c
      call ebtcc(work,u3,work,n3,n3,n2*n1)
      do 90 k=1,n3
         temp = exp(-eye*eig3(k)*t)
         do 100 j=1,n2
            do 200 i=1,n1
               work(k,j,i) = temp*work(k,j,i)
 200        continue
 100     continue
 90   continue   
      call ebcc(psi,u3,work,n3,n3,n2*n1)
c      
      return
      end





















