*deck wrthqq.f
c***begin prologue     wrthqq
c***date written       960615   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c***keywords           hamiltonian
c***author             schneider, barry (nsf)
c***source             3-dim
c***purpose            print non-zero matrix elements
c***                   of hamiltonian and their indices.
c***                   
c***                   
c***references         
c
c***routines called    
c***end prologue       wrthqq
      subroutine wrthqq(ibuf,hbuf,diag,header,nwksq,lenbuf,nonz)
      implicit integer (a-z)
      real*8 hbuf, diag
      character*80 title
      character*(*) header
      dimension hbuf(lenbuf), ibuf(2,lenbuf), diag(*)
      dimension header(2)
      common/io/inp, iout 
      call iosys('read real '//header(1)//' from hamiltonian',
     #            nwksq,diag,0,' ')
      title='diagonal hqq elements'
      call prntrm(title,diag,nwksq,1,nwksq,1,iout)
      if(nonz.ne.0) then
         write(iout,1)
         write(iout,2)
         trips = nonz/lenbuf
         left = nonz - trips*lenbuf
         do 1000 trp=1,trips
            call iosys('read integer '//header(2)//' from '//
     #                 'hamiltonian without rewinding',2*lenbuf,
     #                  ibuf,0,' ') 
            call iosys('read integer '//header(2)//' from '//
     #                 'hamiltonian without rewinding',wptoin(lenbuf),
     #                  hbuf,0,' ')        
            do 2000 i=1,lenbuf
               write(iout,1) ibuf(1,i), ibuf(2,i), hbuf(i)
 2000       continue   
 1000    continue   
         if(left.ne.0) then
            call iosys('read integer '//header(2)//' from '//
     #                 'hamiltonian without rewinding',2*left,ibuf,
     #                  0,' ') 
            call iosys('read integer '//header(2)//' from '//
     #                 'hamiltonian without rewinding',wptoin(left),
     #                  rbuf,0,' ')        
            do 3000 i=1,left
               write(iout,1) ibuf(1,i), ibuf(2,i), hbuf(i)
 3000       continue   
         endif
      endif
      return
 1    format(/,1x,'non-zero matrix elements')
 2    format(/,5x,'    i   ',4x,'   j   ',3x,'  matrix element  ')

      end       


















