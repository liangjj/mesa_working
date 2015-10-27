*deck hmul
c***prologue           hmul
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
c***end prologue       hmul
      subroutine hmul(ibuf,rbuf,pvec,spp,sqp,nwks,nwksp,nwksq,
     #                npvec,number,lenbuf)
      implicit integer (a-z)
      real*8 rbuf, diag, spp, sqp, tmp
      real*8 pvec
      dimension rbuf(lenbuf), ibuf(2,lenbuf)
      dimension pvec(nwks,npvec)
      dimension spp(nwksp,npvec), sqp(nwksq,npvec)
      common/io/inp, iout 
      do 1000 nel=1,number
         i=ibuf(1,nel)
         j=ibuf(2,nel)
         if(i.gt.nwksq.and.j.le.nwksq) then
            do 2000 k=1,npvec
               sqp(j,k) = sqp(j,k) + rbuf(nel)*pvec(i,k)
 2000       continue
         else if(i.gt.nwksq.and.j.gt.nwksq) then
            ii = i - nwksq
            jj = j - nwksq
            do 3000 k=1,npvec
               spp(ii,k) = spp(ii,k) + rbuf(nel)*pvec(j,k)
               spp(jj,k) = spp(jj,k) + rbuf(nel)*pvec(i,k)
 3000       continue
         end if
 1000 continue
      return
      end       


















