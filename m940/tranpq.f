*deck tranpq.f
c***begin prologue     tranpq
c***date written       960615   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c***keywords           hamiltonian
c***author             schneider, barry (nsf)
c***source             3-dim
c***purpose            transform partitioned hamiltonian 
c***                   from primitive to contracted p space.
c***                   
c***                   
c***references         
c
c***routines called    
c***end prologue       tranpq
      subroutine tranpq(ibuf,hbuf,pvec,convec,lenbuf,nwksp,nwksq,ntot,
     1                  headr,prnt)
c
      implicit integer (a-z)
c
      real*8 hbuf, pvec, convec, null
      character*(*) headr
      character*80 title
      character*40 bufr
      logical nodsk, prnt
      dimension ibuf(2,lenbuf), hbuf(lenbuf)
      dimension pvec(nwksp,*), convec(ntot,nwksq), headr(2)
      dimension bufr(2)
      common /io/ inp,iout
      data nodsk/.false./
      data null/1.d-14/
c
      call iosys('rewind all on hamiltonian',0,0,0,' ')
      call iosys('read integer '//headr(2)//' from hamiltonian',
     1            1,nonz,0,' ')
      call hpqmul(ibuf,hbuf,lenbuf,headr(1),nonz,pvec,convec,nwksp,
     1            nwksq,ntot,nodsk)
      if(prnt) then
         title='contracted pq-space hamiltonian'
         call prntrm(title,convec,ntot,nwksq,ntot,nwksq,iout)
      endif
      bufr(1)='"contracted pq space buffers"'
      bufr(2)='"number of contracted pq space elements"'
      call iosys('open scratch as scratch',0,0,0,' ')
      call iosys('create integer '//bufr(1)//' on scratch',
     1           -1,0,0,' ')
      n=0
      nonzpq=0
      do 10 i=1,ntot
         do 20 j=1,nwksq
            if(abs(convec(i,j)).gt.null) then
               n=n+1
               if(n.gt.lenbuf) then
                  nonzpq=nonzpq+lenbuf
                  call iosys('write integer '//bufr(1)//
     1                       ' to scratch '//
     2                       'without rewinding',2*lenbuf,ibuf,0,' ')          
                  call iosys('write integer '//bufr(1)//
     1                       ' to scratch '//
     2                       'without rewinding',wptoin(lenbuf),
     3                        hbuf,0,' ')        
                  n=1
               endif
               ibuf(1,n)=i
               ibuf(2,n)=j
               hbuf(n)=convec(i,j)
            endif
 20      continue         
 10   continue
      if (n.gt.0) then
          nonzpq=nonzpq+n
          call iosys('write integer '//bufr(1)//
     1               ' to scratch without rewinding',
     2                 2*n,ibuf,0,' ')   
          call iosys('write integer '//bufr(1)//
     1               ' to scratch without rewinding',
     2                 wptoin(n),hbuf,0,' ')
      endif
      call iosys('endfile '//bufr(1)//' on scratch',0,0,0,' ')
      call iosys('write integer '//bufr(2)//' to scratch',
     1            1,nonzpq,0,' ')
      call iosys('rewind all on hamiltonian read-and-write',0,0,0,' ')
      call iosys('rewind all on scratch read-and-write',0,0,0,' ')
      write(iout,1) nonzpq
      return
 1    format(/,1x,'number of non-zero contracted pq-space '
     1            'matrix elements = ',i5 )
      end




