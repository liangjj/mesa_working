*deck blkout.f
c***begin prologue     blkout
c***date written       960723   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c***keywords           time development
c***author             schneider, barry (nsf)
c***source             
c***purpose            output trial vectors.
c***                   
c***references         
c
c***routines called    
c***end prologue       blkout
      subroutine blkout(vecin,work,cnt,n,m)
      implicit integer (a-z)
      real*8 vecin, work
      dimension vecin(m,m), work(n)
      common/io/inp, iout
      do 10 i=1,m
         call rzero(work,n)            
         call copy(vecin(1,i),work(cnt+1),msize)
         call iosys('write real vectors to ham without rewinding',
     1               n,work,0,' ')
 10   continue   
      return
      end       
