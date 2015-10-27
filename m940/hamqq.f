*deck hamqq.f
c***begin prologue     hamqq
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
c***end prologue       hamqq
      subroutine hamqq(ibuf,rbuf,diag,iqbuf,hqbuf,diagq,
     1                 hqq,lenbuf,nwksp,nwksq,conqq,title,prnt)
      implicit integer (a-z)
      real*8 rbuf, diag, hqbuf, diagq, hqq
      character*32 title
      character*80 tit
      logical conqq, prnt
      dimension rbuf(lenbuf), ibuf(2,lenbuf), diag(*)
      dimension hqbuf(lenbuf), iqbuf(2,lenbuf), diagq(*)
      dimension hqq(nwksq,nwksq)
      dimension title(3)
      common/io/inp, iout 
      call iosys('rewind all on hamiltonian',0,0,0,' ')
      nwks=nwksp+nwksq
      call iosys('read real diagonals from hamiltonian',nwks,diag,0,' ')
      call copy(diag,diagq,nwksq)
      if(conqq) then
         call rzero(hqq,nwksq*nwksq)
         do 10 i=1,nwksq
            hqq(i,i)=diagq(i)
 10      continue
      endif   
      open(unit=70,file='hamqq',err=99,
     1     form='unformatted',status='unknown')
      rewind(70)
      nonzq=0
      write(70) nwksp, nwksq, nwks, lenbuf, nonzq
      write(70) (diagq(i),i=1,nwksq)
      call iosys('write real '//title(3)//' to hamiltonian',nwksq,
     1            diagq,0,' ')   
      call iosys('read integer "number of elements" from hamiltonian',
     1            1,nonz,0,' ')
      call iosys('create integer '//title(1)//' on hamiltonian',
     1           -1,0,0,' ')
      trips = nonz/lenbuf
      left = nonz - trips*lenbuf
      last=lenbuf 
      if(left.ne.0) then
         trips=trips+1
         last=left
      endif
      n=0
      nonzq=0
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
            if(i.le.nwksq.and.j.le.nwksq) then
               n=n+1
               if (n.gt.lenbuf) then
                   nonzq=nonzq+lenbuf
                   write(70) (iqbuf(1,k),k=1,lenbuf)
                   write(70) (iqbuf(2,k),k=1,lenbuf)
                   call iosys('write integer '//title(1)//
     $                        ' to hamiltonian without rewinding',
     $                         2*lenbuf,iqbuf,0,' ')
                   call iosys('write integer '//title(1)//
     $                        ' to hamiltonian without rewinding',
     $                         wptoin(lenbuf),hqbuf,0,' ')
                   write(70) (hqbuf(k),k=1,lenbuf)
                   n=1
               endif
               iqbuf(1,n)=i
               iqbuf(2,n)=j
               hqbuf(n)=rbuf(nel)
               if(conqq) then
                  hqq(iqbuf(1,n),iqbuf(2,n)) = 
     1            hqq(iqbuf(1,n),iqbuf(2,n)) + hqbuf(n) 
               endif  
            endif
 30      continue   
 20   continue   
      if (n.gt.0) then
          nonzq=nonzq+n
          call iosys('write integer '//title(1)//' to hamiltonian'//
     $               ' without rewinding',2*n,iqbuf,0,' ')
          write(70) (iqbuf(1,k),k=1,n)
          write(70) (iqbuf(2,k),k=1,n)
          call iosys('write integer '//title(1)//' to hamiltonian'//
     $               ' without rewinding',wptoin(n),hqbuf,0,' ')
          write(70) (hqbuf(k),k=1,n)
      endif
      rewind(70)
      write(70) nwksp, nwksq, nwks, lenbuf, nonzq
      close(70)
      call iosys('endfile '//title(1)//' on hamiltonian',0,0,0,' ')
      call iosys('write integer '//title(2)//' to hamiltonian',
     1            1,nonzq,0,' ')
      call iosys('rewind all on hamiltonian read-and-write',0,0,0,' ')
      write(iout,1) nonzq
      if(prnt) then
         call rdham(hqbuf,iqbuf,diagq,lenbuf,nwksq,title)
      endif
      if(conqq) then
         do 50 i=1,nwksq
            do 60 j=1,i
               hqq(j,i)=hqq(i,j)
 60         continue
 50      continue
         call iosys('write real h(qq) to hamiltonian',nwksq*nwksq,
     1               hqq,0,' ')
         if(prnt) then
            tit='q-space hamiltonian'
            call prntrm(tit,hqq,nwksq,nwksq,nwksq,nwksq,iout)
         endif
      endif 
 99   continue
      return
 1    format(/,1x,'number of non-zero q-space matrix elements = ',i5 )
      end       
