*deck ham.f
c***begin prologue     ham
c***date written       960615   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c***keywords           hamiltonian
c***author             schneider, barry (nsf)
c***source             3-dim
c***purpose            non-zero matrix elements of partitioned
c***                   hamiltonian and their indices.
c***                   
c***                   
c***references         
c
c***routines called    
c***end prologue       ham
      subroutine ham(ibuf,rbuf,diag,h,lenbuf,nwks,title,prnt)
      implicit integer (a-z)
      real*8 rbuf, diag, h
      character*80 tit
      logical prnt
      character*32 title
      dimension rbuf(lenbuf), ibuf(2,lenbuf), diag(*)
      dimension h(nwks,nwks)
      dimension title(3)
      common/io/inp, iout 
      call iosys('rewind all on hamiltonian',0,0,0,' ')
      open(unit=80,file='hamiltonian',err=99,
     1     form='unformatted',status='unknown')
      rewind(80)
      call iosys('read real '//title(3)//' from hamiltonian',nwks,
     1            diag,0,' ')
      call rzero(h,nwks*nwks)
      do 10 i=1,nwks
         h(i,i)=diag(i)
 10   continue
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
         do 30 nel=1,number
            i=ibuf(1,nel)
            j=ibuf(2,nel)
            h(i,j) =  h(i,j) + rbuf(nel) 
 30      continue   
 20   continue
      rewind(80)
      close(80)
      call iosys('rewind all on hamiltonian read-and-write',0,0,0,' ')
      do 40 i=1,nwks
         do 50 j=1,i
            h(j,i)=h(i,j)
 50      continue
 40   continue
      call iosys('write real h to hamiltonian',nwks*nwks,h,0,' ')
      if(prnt) then
         tit='full hamiltonian'
         call prntrm(tit,h,nwks,nwks,nwks,nwks,iout)
      endif
 99   continue
      return
      end       


















