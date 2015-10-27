*deck eshft
c***prologue           eshft
c***date written       960615   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c***keywords           hamiltonian
c***author             schneider, barry (nsf)
c***source             3-dim
c***purpose            form E - H
c***                   
c***                   
c***                   
c***references         
c
c***routines called    
c***end prologue       eshft
      subroutine eshft(hfull,htri,energy,n,ntri,type)
      implicit integer (a-z)
      real*8 hfull, htri, energy
      character*(*) type
      dimension hfull(n,n), htri(ntri)
      common/io/inp, iout 
      if(type.eq.'triangle') then
         call vneg(htri,htri,ntri)
         ii=1
         htri(ii) = energy + htri(ii)
         do 1000 i=2,n
            ii=ii+i
            htri(ii) = energy + htri(ii)
 1000    continue   
      else
         call vneg(hfull,hfull,n*n) 
         do 2000 i=1,n
            hfull(i,i) = energy + hfull(i,i)
 2000    continue   
      endif   
      return
      end       
