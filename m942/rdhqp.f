*deck rdhqp.f
c***begin prologue     rdhqp
c***date written       960615   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c***keywords           rdhqpiltonian
c***author             schneider, barry (nsf)
c***source             3-dim
c***purpose            fill hqp using its non-zero elements
c***                   
c***description                           
c***                   
c***references         
c
c***routines called    
c***end prologue       rdhqp
      subroutine rdhqp(ibuf,rbuf,hqp,header,n,npvec,lenbuf,
     #                 nonz,nodsk,prnt)
      implicit integer (a-z)
      real*8 rbuf, hqp
      logical nodsk, prnt
      character*(*) header
      character*80 title
      dimension rbuf(lenbuf), ibuf(2,lenbuf)
      dimension hqp(n,npvec)
      common/io/inp, iout 
      call rzero(hqp,n*npvec)
      if(nonz.gt.0) then
         if(nodsk) then
            call iosys('read integer '//header//' from '//
     #                 'hamiltonian without rewinding',
     #                  2*nonz,ibuf,0,' ')
            call iosys('read integer '//header//' from '//
     #                 'hamiltonian without rewinding',
     #                  wptoin(nonz),rbuf,0,' ')
            call hqpfill(ibuf,rbuf,hqp,n,npvec,nonz)
         else
            trips = nonz/lenbuf
            left = nonz - trips*lenbuf
            do 1000 trp=1,trips
               call iosys('read integer '//header//' from '//
     #                    'hamiltonian without rewinding',
     #                     2*lenbuf,ibuf,0,' ')
               call iosys('read integer '//header//' from '//
     #                    'hamiltonian without rewinding',
     #                     wptoin(lenbuf),rbuf,0,' ')
               call hqpfill(ibuf,rbuf,hqp,n,npvec,lenbuf)
 1000       continue   
            if(left.ne.0) then
               call iosys('read integer '//header//' from '//
     #                    'hamiltonian without rewinding',
     #                     2*left,ibuf,0,' ')
               call iosys('read integer '//header//' from '//
     #                    'hamiltonian without rewinding',
     #                     wptoin(left),rbuf,0,' ')
               call hqpfill(ibuf,rbuf,hqp,n,npvec,left)
            end if
         endif
      endif
      if(prnt) then
         title='H_QP Matrix Elements'
         call prntrm(title,hqp,n,npvec,n,npvec,iout)
      endif
      call iosys('rewind '//header//' on hamiltonian '//
     #           'read-and-write',0,0,0,' ')
      return
      end       
