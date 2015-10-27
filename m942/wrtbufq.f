*deck wrtbufq.f
c***begin prologue     wrtbufq
c***date written       960615   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c***keywords           wrtbufqiltonian
c***author             schneider, barry (nsf)
c***source             3-dim
c***purpose            non-zero elements of hqq
c***                   
c***description                           
c***                   
c***references         
c
c***routines called    
c***end prologue       wrtbufq
      subroutine wrtbufq(ibuf,rbuf,ibufq,rbufq,header,nwksq,num,lenbuf,
     #                   n,nonz)
      implicit integer (a-z)
      real*8 rbuf, rbufq
      dimension rbuf(lenbuf), ibuf(2,lenbuf)
      dimension rbufq(lenbuf), ibufq(2,lenbuf)
      character*(*) header
      common/io/inp, iout 
      do 1000 nel=1,num
         if(ibuf(1,nel).le.nwksq) then
            n=n+1  
            if(n.gt.lenbuf) then
               nonz=nonz+lenbuf
               call iosys('write integer '//header//
     #                    ' to hamiltonian '//
     #                    'without rewinding',2*lenbuf,ibufq,0,' ') 
               call iosys('write integer '//header//
     #                    ' to hamiltonian '//
     #                    'without rewinding',wptoin(lenbuf),
     #                     rbufq,0,' ')        
               n=1
            endif
            ibufq(1,n)=ibuf(1,nel)
            ibufq(2,n)=ibuf(2,nel)
            rbufq(n)=rbuf(nel)
         endif
 1000 continue   
      return
      end       


















