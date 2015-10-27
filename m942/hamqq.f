*deck hamqq.f
c***begin prologue     hamqq
c***date written       960615   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c***keywords           hamqqiltonian
c***author             schneider, barry (nsf)
c***source             3-dim
c***purpose            non-zero elements of hqq
c***                   
c***description                           
c***                   
c***references         
c
c***routines called    
c***end prologue       hamqq
      subroutine hamqq(ibuf,rbuf,ibufq,rbufq,diagq,header,nwksq,
     #                 lenbuf,nonz,nodsk,prnt)
      implicit integer (a-z)
      real*8 rbuf, rbufq, diagq
      logical nodsk, prnt
      dimension rbuf(lenbuf), ibuf(2,lenbuf)
      dimension rbufq(lenbuf), ibufq(2,lenbuf)
      dimension diagq(nwksq)
      character*(*) header(10)
      common/io/inp, iout 
      call iosys('create integer '//header(9)//' on hamiltonian',
     #            -1,0,0,' ')
      n=0
      nonzqq=0
      if(nodsk) then
         call rdham(ibuf,rbuf,header(1),lenbuf,nonz) 
         call wrtbufq(ibuf,rbuf,ibufq,rbufq,header(9),nwksq,nonz,
     #                lenbuf,n,nonzqq)
         if(n.gt.0) then
            call lstwrt(ibufq,rbufq,header(9),lenbuf,n)
            nonzqq=nonzqq+n
         end if
      else
         trips = nonz/lenbuf
         left = nonz - trips*lenbuf
         do 1000 trp=1,trips
            call iosys('read integer '//header(1)//
     #                 ' from hamiltonian without rewinding',
     #                   2*lenbuf,ibuf,0,' ')
            call iosys('read integer '//header(1)//
     #                 ' from hamiltonian without rewinding',
     #                   wptoin(lenbuf),rbuf,
     #                   0,' ')
            call wrtbufq(ibuf,rbuf,ibufq,rbufq,header(9),nwksq,
     #                   lenbuf,lenbuf,n,nonzqq)
 1000    continue   
         if(left.ne.0) then
            call iosys('read integer '//header(1)//
     #                 ' from hamiltonian without rewinding',
     #                   2*left,ibuf,0,' ')
            call iosys('read integer '//header(1)//
     #                 ' from hamiltonian without rewinding',
     #                   wptoin(left),rbuf,0,' ')
            call wrtbufq(ibuf,rbuf,ibufq,rbufq,header(9),
     #                   nwksq,left,lenbuf,n,nonzqq)
            if(n.gt.0) then
               call lstwrt(ibufq,rbufq,header(9),lenbuf,n)
               nonzqq=nonzqq+n
            end if
         end if
      endif
      call iosys('endfile '//header(9)//
     #           ' on hamiltonian',0,0,0,' ')
      call iosys('write integer '//header(10)//
     #           ' to hamiltonian',1,nonzqq,0,' ')
      write(iout,1) nonzqq
      if(prnt) then
         call wrthqq(ibufq,hbufq,diagq,header(8),nwksq,lenbuf,nonzqq)
      endif 
      return
 1    format(/,1x,'number of non-zero q-space '
     1            'matrix elements = ',i5 )
 2    format(/,1x,a32)
      end       
