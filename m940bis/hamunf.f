*deck hamunf.f
c***begin prologue     hamunf
c***date written       960615   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c***keywords           hamiltonian
c***author             schneider, barry (nsf)
c***source             3-dim
c***purpose            non-zero matrix elements of
c***                   hamiltonian and their indices to
c***                   unformatted file.
c***                   
c***                   
c***references         
c
c***routines called    
c***end prologue       hamunf
      subroutine hamunf(ibuf,rbuf,diag,lenbuf,nwks,title)
      implicit integer (a-z)
      real*8 rbuf, diag
      character*32 title
      dimension rbuf(lenbuf), ibuf(2,lenbuf), diag(*)
      dimension title(3)
      common/io/inp, iout 
      call iosys('rewind all on hamiltonian',0,0,0,' ')
      open(unit=80,file='hamiltonian',err=99,
     1     form='unformatted',status='unknown')
      rewind(80)
      call iosys('read real '//title(3)//' from hamiltonian',nwks,
     1            diag,0,' ')
      call iosys('read integer '//title(2)//' from hamiltonian',
     1            1,nonz,0,' ')
      write(80) nwks, lenbuf, nonz
      write(80) (diag(i),i=1,nwks)
      trips = nonz/lenbuf
      left = nonz - trips*lenbuf
      last=lenbuf 
      if(left.ne.0) then
         trips=trips+1
         last=left
      endif
      do 20 trp=1,trips
         number=lenbuf
         if(trp.eq.trips) then
            number=last
         endif
         call iosys('read integer '//title(1)//' from hamiltonian '//
     1              'without rewinding',2*number,ibuf,0,' ')
         write(80) (ibuf(1,k),k=1,number)
         write(80) (ibuf(2,k),k=1,number)
         call iosys('read integer '//title(1)//' from hamiltonian '//
     1              'without rewinding',wptoin(number),rbuf,0,' ')
         write(80) (rbuf(k),k=1,number)
 20   continue
      rewind(80)   
      close(80)
      call iosys('rewind all on hamiltonian read-and-write',0,0,0,' ')
 99   continue
      return
      end       


















