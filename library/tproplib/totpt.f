*deck totpt.f
c***begin prologue     totpt
c***date written       951229   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c***keywords           
c***author             schneider, barry (nsf)
c***source             
c***purpose            
c***                   
c***references         
c
c***routines called    
c***end prologue       totpt
      subroutine totpt(npt,nptsr,nreg,qi,prnt)
      implicit integer (a-z)
      logical prnt
      dimension nptsr(*)
      common/io/inp, iout
      if(prnt) then
         write(iout,1) qi
      endif
      npt=0
      do 10 i=1,nreg
         if(prnt) then 
            write(iout,2) i, nptsr(i)
         endif
         npt=npt+nptsr(i)
 10   continue
      if(prnt) then
         write(iout,3) npt   
      endif
      return
 1    format(/,15x,'point information for dimension = ',i1)
 2    format(/,1x,'region = ',i2,1x,'number of points = ',i4)
 3    format(/,15x,'total number of points = ',i4)
      end       
