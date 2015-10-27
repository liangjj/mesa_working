*deck trimm
c***begin prologue     trimm
c***date written       920209   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c***keywords           trim a matrix
c***author             schneider, barry (nsf)
c***source             mylib 
c***                                         
c***purpose            trim the leading or trailing parts of a matrix
c***                   
c***                  
c***references         none
c
c***routines called    
c***end prologue       trimm
      subroutine trimm(min,mout,mscr,keep,nkeep,n)
      implicit integer (a-z)
      real*8 min, mout, mscr
      dimension min(n,n), mout(nkeep,nkeep), mscr(nkeep,nkeep)
      dimension keep(nkeep)
      do 10 i=1,nkeep
         do 20 j=1,nkeep
            mscr(i,j)=min(keep(i),keep(j))         
 20      continue   
 10   continue
      call copy(mscr,mout,nkeep*nkeep) 
      return
      end
