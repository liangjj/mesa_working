*deck hampp.f
c***begin prologue     hampp
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
c***end prologue       hampp
      subroutine hampp(ibuf,rbuf,diag,ipbuf,hpbuf,diagp,
     1                 hpp,lenbuf,nwksp,nwksq,conpp,title,prnt)
      implicit integer (a-z)
      real*8 rbuf, diag, hpbuf, diagp, hpp
      character*32 title
      character*80 tit
      logical conpp, prnt
      dimension rbuf(lenbuf), ibuf(2,lenbuf), diag(*)
      dimension hpbuf(lenbuf), ipbuf(2,lenbuf), diagp(*)
      dimension hpp(nwksp,nwksp)
      dimension title(3)
      common/io/inp, iout 
      call iosys('rewind all on hamiltonian',0,0,0,' ')
      open(unit=90,file='hampp',err=99,
     1     form='unformatted',status='unknown')
      rewind(90)
      nwks=nwksp+nwksq
      nonzp=0
      write(90) nwksp, nwksq, nwks, lenbuf, nonzp
      call iosys('read real diagonals from hamiltonian',nwks,diag,0,' ')
      count=0
      do 10 i=nwksq+1,nwks
         count=count+1
         diagp(count)=diag(i)
 10   continue
      write(90) (diagp(i),i=1,nwksp)
      if(conpp) then
         call rzero(hpp,nwksp*nwksp)
         do 20 i=1,nwksp
            hpp(i,i)=diagp(i)
 20      continue
      endif
      call iosys('write real '//title(3)//' to hamiltonian',nwksp,
     1            diagp,0,' ')   
      call iosys('create integer '//title(1)//' on hamiltonian',
     1           -1,0,0,' ')
      call iosys('read integer "number of elements" from hamiltonian',
     1            1,nonz,0,' ')
      trips = nonz/lenbuf
      left = nonz - trips*lenbuf
      last=lenbuf 
      if(left.ne.0) then
         trips=trips+1
         last=left
      endif
      n=0
      nonzp=0
      do 30 trp=1,trips
         number=lenbuf
         if(trp.eq.trips) then
            number=last
         endif
         call iosys('read integer buffers from hamiltonian '//
     1              'without rewinding',2*number,ibuf,0,' ')
         call iosys('read integer buffers from hamiltonian '//
     1              'without rewinding',wptoin(number),rbuf,0,' ')
         do 40 nel=1,number
            i=ibuf(1,nel)
            j=ibuf(2,nel)
            if(i.gt.nwksq.and.j.gt.nwksq) then
               n=n+1
               if (n.gt.lenbuf) then
                   nonzp=nonzp+lenbuf
                   call iosys('write integer '//title(1)//
     $                        ' to hamiltonian without rewinding',
     $                          2*lenbuf,ipbuf,0,' ')
                   write(90) (ipbuf(1,k),k=1,lenbuf)
                   write(90) (ipbuf(2,k),k=1,lenbuf)
                   call iosys('write integer '//title(1)//
     $                        ' to hamiltonian without rewinding',
     $                          wptoin(lenbuf),hpbuf,0,' ')
                   write(90) (hpbuf(k),k=1,lenbuf)
                   n=1
               endif
               ipbuf(1,n)=i-nwksq
               ipbuf(2,n)=j-nwksq
               hpbuf(n)=rbuf(nel)
               if(conpp) then
                  hpp(ipbuf(1,n),ipbuf(2,n)) = 
     1            hpp(ipbuf(1,n),ipbuf(2,n)) + hpbuf(n) 
               endif  
            endif
 40      continue   
 30   continue   
      if (n.gt.0) then
          nonzp=nonzp+n
          call iosys('write integer '//title(1)//
     1               ' to hamiltonian without rewinding',
     2                 2*n,ipbuf,0,' ')
          write(90) (ipbuf(1,k),k=1,n)
          write(90) (ipbuf(2,k),i=k,n)
          call iosys('write integer '//title(1)//
     1               ' to hamiltonian without rewinding',
     2                 wptoin(n),hpbuf,0,' ')
          write(90) (hpbuf(k),k=1,n)
      endif
      rewind(90)
      write(90) nwksp, nwksq, nwks, lenbuf, nonzp
      close(90)
      call iosys('endfile '//title(1)//' on hamiltonian',0,0,0,' ')
      call iosys('write integer '//title(2)//
     1           ' to hamiltonian',1,nonzp,0,' ')
      call iosys('rewind all on hamiltonian read-and-write',0,0,0,' ')
      write(iout,1) nonzp
      if(prnt) then
         call rdham(hpbuf,ipbuf,diagp,lenbuf,nwksp,title)
      endif
      if(conpp) then
         do 50 i=1,nwksp
            do 60 j=1,i
               hpp(j,i)=hpp(i,j)
 60         continue
 50      continue
         call iosys('write real h(pp) to hamiltonian',nwksp*nwksp,
     1               hpp,0,' ')
         if(prnt) then
            tit='p-space hamiltonian'
            call prntrm(tit,hpp,nwksp,nwksp,nwksp,nwksp,iout)
         endif
      endif 
 99   continue
      return
 1    format(/,1x,'number of non-zero p-space matrix elements = ',i10 )
      end       


















