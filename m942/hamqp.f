*deck hamqp.f
c***begin prologue     hamqp
c***date written       960615   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c***keywords           hamiltonian
c***author             schneider, barry (nsf)
c***source             3-dim
c***purpose            out of core matrix vector multiply
c***                   using non-zero elements of matrix.
c***description                           
c***                   
c***references         
c
c***routines called    
c***end prologue       hamqp
      subroutine hamqp(ibuf,rbuf,hqp,pvec,header,nwksp,nwksq,nwks,
     #                 npvec,lenbuf,nonz,nodsk,prnt)
      implicit integer (a-z)
      real*8 rbuf, hqp, pvec
      character*80 title
      character*(*) header
      logical nodsk, prnt
      dimension rbuf(lenbuf), ibuf(2,lenbuf)
      dimension hqp(nwksq,npvec)
      dimension pvec(nwks,npvec), header(10)
      common/io/inp, iout 
      data null / 1.d-14 /
      if(nodsk) then
         call hmulqp(ibuf,rbuf,pvec,hqp,nwks,nwksq,npvec,nonz,lenbuf)
      else
         call iosys('rewind '//header(1)//' on hamiltonian '//
     #              'read-and-write',0,0,0,' ')
         trips = nonz/lenbuf
         left = nonz - trips*lenbuf
         do 1000 trp=1,trips
            call iosys('read integer '//header(1)//' from '//
     #                 'hamiltonian without rewinding',
     #                  2*lenbuf,ibuf,0,' ')
            call iosys('read integer '//header(1)//' from '//
     #                 'hamiltonian without rewinding',
     #                  wptoin(lenbuf),rbuf,0,' ')
            call hmulqp(ibuf,rbuf,pvec,hqp,nwks,nwksp,nwksq,
     #                  npvec,lenbuf,lenbuf)
 1000    continue
         if(left.ne.0) then
            call iosys('read integer '//header(1)//' from '//
     #                 'hamiltonian without rewinding',
     #                  2*left,ibuf,0,' ')
            call iosys('read integer '//header(1)//' from '//
     #                 'hamiltonian without rewinding',
     #                  wptoin(left),rbuf,0,' ')
            call hmulqp(ibuf,rbuf,pvec,hqp,nwks,nwksq,npvec,left,lenbuf)
         end if
      end if
      call iosys('create integer '//header(5)//' on hamiltonian',
     #            -1,0,0,' ')
      n=0
      nonzqp=0
      do 2000 i=1,nwksq
         do 3000 j=1,npvec
            if(abs(hqp(i,j)).gt.null) then
               n=n+1
               if(n.gt.lenbuf) then
                  nonzqp=nonzqp+lenbuf
                  call iosys('write integer '//header(5)//
     #                       ' to hamiltonian '//
     #                       'without rewinding',2*lenbuf,ibuf,0,' ') 
                  call iosys('write integer '//header(5)//
     #                       ' to hamiltonian '//
     #                       'without rewinding',wptoin(lenbuf),
     #                        rbuf,0,' ')        
                  n=1
               endif
               ibuf(1,n)=i
               ibuf(2,n)=j
               rbuf(n)=hqp(i,j)
            endif
 3000    continue         
 2000 continue
      if (n.gt.0) then
          nonzqp=nonzqp+n
          call iosys('write integer '//header(5)//
     #               ' to hamiltonian '//
     #               'without rewinding',2*n,ibuf,0,' ') 
          call iosys('write integer '//header(5)//
     #               ' to hamiltonian '//
     #               'without rewinding',wptoin(n),
     #                rbuf,0,' ')        
      endif
      call iosys('endfile '//header(5)//
     #           ' on hamiltonian',0,0,0,' ')
      call iosys('write integer '//header(6)//' to hamiltonian',
     #            1,nonzqp,0,' ')
      call iosys('write real H_QP on hamiltonian without rewinding',
     #            npvec*nwksq,hqp,0,' ')
      if(prnt) then
         title='Contracted H_QP'
         call prntrm(title,hqp,nwksq,npvec,nwksq,npvec,iout)
      endif
      write(iout,1) nonzqp
      return
 1    format(/,1x,'number of non-zero contracted qp-space '
     1            'matrix elements = ',i5 )
      end       


















