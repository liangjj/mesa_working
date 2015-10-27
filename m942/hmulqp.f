*deck hmulqp
c***prologue           hmulqp
c***date written       960615   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c***keywords           hamiltonian
c***author             schneider, barry (nsf)
c***source             3-dim
c***purpose            matrix vector multiply
c***                   hamiltonian and their indices.
c***                   
c***                   
c***references         
c
c***routines called    
c***end prologue       hmulqp
      subroutine hmulqp(ibuf,rbuf,pvec,sqp,nwks,nwksq,npvec,
     #                  number,lenbuf)
      implicit integer (a-z)
      real*8 rbuf, sqp
      real*8 pvec
      dimension rbuf(lenbuf), ibuf(2,lenbuf)
      dimension pvec(nwks,npvec)
      dimension sqp(nwksq,npvec)
      do 1000 nel=1,number
         i=ibuf(1,nel)
         j=ibuf(2,nel)
         if(i.gt.nwksq.and.j.le.nwksq) then
            do 2000 k=1,npvec
               sqp(j,k) = sqp(j,k) + rbuf(nel)*pvec(i,k)
 2000       continue
         end if
 1000 continue
      return
      end       


















