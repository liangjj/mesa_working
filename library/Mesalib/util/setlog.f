*deck setlog.f
c***begin prologue     setlog
c***date written       951229   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c***keywords           
c***author             schneider, barry (nsf)
c***source             
c***purpose            set logical variables depending on input.
c***                   
c***references         
c
c***routines called    
c***end prologue       setlog
      subroutine setlog(prn,key,card,n)
      implicit integer (a-z)
      logical prn, logkey
      character*(*) key, card
      dimension prn(n), key(n) 
      common/io/inp, iout
      do 10 i=1,n
         prn(i)=logkey(card,key(i),.false.,' ')
 10   continue   
      return      
      end       






