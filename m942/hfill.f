*deck hfill
c***prologue           hfill
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
c***end prologue       hfill
      subroutine hfill(ibuf,rbuf,hfull,htri,n,number,type)
      implicit integer (a-z)
      real*8 rbuf, hfull, htri
      character*(*) type
      dimension rbuf(number), ibuf(2,number)
      dimension hfull(n,n), htri(n*(n+1)/2)
      common/io/inp, iout 
      if(type.eq.'triangle') then
         do 1000 nel=1,number
            i=ibuf(1,nel)
            j=ibuf(2,nel)
            ii= i*(i-1)/2 + j          
            htri(ii)=htri(ii) + rbuf(nel)
 1000    continue   
      else
         do 2000 nel=1,number
            i=ibuf(1,nel)
            j=ibuf(2,nel)
            hfull(i,j)=hfull(i,j) + rbuf(nel)
 2000    continue   
      endif   
      return
      end       
