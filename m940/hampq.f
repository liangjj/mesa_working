*hampp.f
deck hampq.f
c***begin prologue     hampq
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
c***end prologue       hampq
      subroutine hampq(ibuf,rbuf,ipqbuf,hpqbuf,hpq,lenbuf,
     1                 nwksp,nwksq,conpq,title,prnt)
      implicit integer (a-z)
      real*8 rbuf, hpqbuf, hpq
      logical conpq, prnt
      character*80 tit
      character*32 title
      dimension rbuf(lenbuf), ibuf(2,lenbuf), hpq(nwksp,nwksq)
      dimension hpqbuf(lenbuf), ipqbuf(2,lenbuf), title(2)
      common/io/inp, iout 
      call iosys('rewind all on hamiltonian',0,0,0,' ')
      open(unit=60,file='hampq',err=99,
     1     form='unformatted',status='unknown')
      rewind(60)
      nwks=nwksp+nwksq
      nonzpq=0
      write(60) nwksp, nwksq, nwks, lenbuf, nonzpq
      call iosys('read integer "number of elements" from hamiltonian',
     1            1,nonz,0,' ')
      call iosys('create integer '//title(1)//' on hamiltonian',
     1           -1,0,0,' ')
      if(conpq) then
         call rzero(hpq,nwksp*nwksq)
      endif
      trips = nonz/lenbuf
      left = nonz - trips*lenbuf
      last=lenbuf 
      if(left.ne.0) then
         trips=trips+1
         last=left
      endif
      n=0
      nonzpq=0
      do 20 trp=1,trips
         number=lenbuf
         if(trp.eq.trips) then
            number=last
         endif
         call iosys('read integer buffers from hamiltonian '//
     1              'without rewinding',2*number,ibuf,0,' ')
         call iosys('read integer buffers from hamiltonian '//
     1              'without rewinding',wptoin(number),rbuf,0,' ')
         do 30 nel=1,number
            i=ibuf(1,nel)
            j=ibuf(2,nel)
            if(i.gt.nwksq.and.j.le.nwksq) then
               n=n+1
               if (n.gt.lenbuf) then
                   nonzpq=nonzpq+lenbuf
                   call iosys('write integer '//title(1)//
     $                        ' to hamiltonian without rewinding',
     $                         2*lenbuf,ipqbuf,0,' ')
                   write(60) (ipqbuf(1,k),k=1,lenbuf)
                   write(60) (ipqbuf(2,k),k=1,lenbuf)
                   call iosys('write integer '//title(1)//
     $                        ' to hamiltonian without rewinding',
     $                         wptoin(lenbuf),hpqbuf,0,' ')
                   write(60) (hpqbuf(k),k=1,lenbuf)
                   n=1
               endif
               ipqbuf(1,n)=i-nwksq
               ipqbuf(2,n)=j
               hpqbuf(n)=rbuf(nel)
               if(conpq) then
                  hpq(ipqbuf(1,n),ipqbuf(2,n)) = 
     1            hpq(ipqbuf(1,n),ipqbuf(2,n)) + hpqbuf(n)
               endif 
            endif
 30      continue   
 20   continue   
      if (n.gt.0) then
          nonzpq=nonzpq+n
          call iosys('write integer '//title(1)//' to hamiltonian'//
     $               ' without rewinding',2*n,ipqbuf,0,' ')
          write(60) (ipqbuf(1,k),k=1,n)
          write(60) (ipqbuf(2,k),k=1,n)
          call iosys('write integer '//title(1)//' to hamiltonian'//
     $               ' without rewinding',wptoin(n),hpqbuf,0,' ')
          write(60) (hpqbuf(k),k=1,n)
      endif
      rewind(60)
      write(60) nwksp, nwksq, nwks, lenbuf, nonzpq
      close(60)
      call iosys('endfile '//title(1)//' on hamiltonian',0,0,0,' ')
      call iosys('write integer '//title(2)//' to hamiltonian',
     1            1,nonzpq,0,' ')
      call iosys('rewind all on hamiltonian read-and-write',0,0,0,' ')
      write(iout,1) nonzpq
      call rdham(hpqbuf,ipqbuf,hpqbuf,lenbuf,0,title)
      if(conpq) then
         call iosys('write real h(pq) to hamiltonian',nwksp*nwksq,
     1               hpq,0,' ')
         if(prnt) then
            tit='pq-space hamiltonian'
            call prntrm(tit,hpq,nwksp,nwksq,nwksp,nwksq,iout)
         endif
      endif
 99   continue
      return
 1    format(/,1x,'number of non-zero pq-space matrix elements = ',i5 )
      end       


















