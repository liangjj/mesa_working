*deck gwadd.f
c***begin prologue     gwadd
c***date written       951229   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c***keywords           
c***author             schneider, barry (nsf)
c***source             
c***purpose            memory outlay for a particular grid.
c***                   
c***references         
c
c***routines called    
c***end prologue       gwadd
 1    subroutine gwadd(q,wt,eigc,wtc,dim,words)
      implicit integer (a-z)
      common/io/inp, iout
      ibeg=words
      q=ibeg
      wt=q+dim
      eigc=wt+dim
      wtc=eigc+dim
      ibeg=wtc+dim
      words=ibeg
      return
      end       
