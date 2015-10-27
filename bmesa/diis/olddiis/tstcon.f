*deck tstcon.f
c***begin prologue     tstcon
c***date written       960723   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c***keywords           
c***author             schneider, barry (nsf)
c***source             
c***purpose            
c***                   
c***description          
c***references         
c
c***routines called    
c***end prologue       tstcon
      subroutine tstcon(v1,v2,scr,ham,psi,error,type,nobj,n)
      implicit integer (a-z)
      real*8 v1, v2, scr, ham, psi, error, sdot
      character*(*) type
      dimension v1(*), v2(*), scr(*), ham (n,n), psi(*)
      common/io/ inp, iout
      if(type.eq.'fock matrix') then
         ij=0
         do 10 i=1,n
            do 20 j=1,i
               ij = ij+1 
               v2(ij) = ham(i,j)  
 20         continue
 10      continue
      elseif(type.eq.'wavefunction') then
          call copy(psi,v2,n)
c          write(iout,*) (v1(ii),ii=1,nobj)
c          write(iout,*) (v2(ii),ii=1,nobj)
      elseif(type.eq.'hartree') then
          call copy(scr,v2,n)
      endif          
      do 30 i=1,nobj
         scr(i) = v1(i) - v2(i)
 30   continue
      error=sqrt( sdot(nobj,scr,1,scr,1) )   
      write(iout,1) error
      return
 1    format(/1x,'error from tstcon = ',e15.8)
      end       

