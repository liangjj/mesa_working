*deck hqpfill
c***prologue           hqpfill
c***date written       960615   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c***keywords           hamiltonian
c***author             schneider, barry (nsf)
c***source             3-dim
c***purpose            fill a hamiltonian using the
c***                   non zero elements and their indices.
c***                   
c***                   
c***references         
c
c***routines called    
c***end prologue       hqpfill
      subroutine hqpfill(ibuf,rbuf,h,n,npvec,number)
      implicit integer (a-z)
      real*8 rbuf, h
      dimension rbuf(number), ibuf(2,number)
      dimension h(n,npvec)
      common/io/inp, iout 
      do 1000 nel=1,number
            i=ibuf(1,nel)
            j=ibuf(2,nel)
            h(i,j)=h(i,j) + rbuf(nel)
 1000 continue   
      return
      end       
