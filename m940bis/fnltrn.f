*deck fnltrn.f
c***begin prologue     fnltran
c***date written       000921   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c***keywords           hamiltonian
c***author             schneider, barry (nsf)
c***source             3-dim
c***purpose            transform partitioned hamiltonian 
c***                   from primitive to contracted q space.
c***                   
c***                   
c***references         
c
c***routines called    
c***end prologue       fnltrn
      subroutine fnltrn(ibuf,hbuf,qvec,convec,lenbuf,nwksp,nwksq,
     1                  qroots,ntot,prnt)
c
      implicit integer (a-z)
c
      real*8 hbuf, qvec, convec, null
      character*80 title
      character*40 bufr
      character*4 itoc
      logical prnt
      logical nodsk
      dimension ibuf(2,lenbuf), hbuf(lenbuf)
      dimension qvec(nwksq,qroots), convec(ntot,qroots)
      dimension bufr(2)
      common /io/ inp,iout
      data nodsk/.false./
      data null/1.d-14/
c
      do 10 i=1,qroots
         call iosys('read real "q-ci root '//itoc(i)//'" from rwf',
     1               nwksq,qvec(1,i),0,' ')
 10   continue   
      call rzero(convec,ntot*qroots)
      call iosys('rewind all on scratch',0,0,0,' ')
      bufr(1)='"contracted pq space buffers"'
      bufr(2)='"number of contracted pq space elements"'
      call iosys('read integer '//bufr(2)//' from scratch',
     1            1,nel,0,' ')
      trips=nel/lenbuf
      left=nel-trips*lenbuf
      do 20 nt=1,trips
         call iosys('read integer '//bufr(1)//' from scratch '//
     1              'without rewinding',2*lenbuf,ibuf,0,' ')
         call iosys('read integer '//bufr(1)//' from scratch '//
     1              'without rewinding',wptoin(lenbuf),hbuf,0,' ')
         do 30 i=1,lenbuf
            ii=ibuf(1,i)
            jj=ibuf(2,i)
            do 40 j=1,qroots
               convec(ii,j) = convec(ii,j) + hbuf(i)*qvec(jj,j)
 40         continue
 30      continue
 20   continue
      if(left.ne.0) then
         call iosys('read integer '//bufr(1)//' from scratch '//
     1              'without rewinding',2*left,ibuf,0,' ')
         call iosys('read integer '//bufr(1)//' from scratch '//
     1              'without rewinding',wptoin(left),hbuf,0,' ')
         do 50 i=1,left
            ii=ibuf(1,i)
            jj=ibuf(2,i)
            do 60 j=1,qroots
               convec(ii,j) = convec(ii,j) + hbuf(i)*qvec(jj,j)
 60         continue
 50      continue
      endif   
      if(prnt) then
         title='contracted p-and-q-space hamiltonian'
         call prntrm(title,convec,ntot,qroots,ntot,qroots,iout)
      endif
      call iosys('destroy scratch',0,0,0,' ')
      call iosys('create integer '//bufr(1)//' on hamiltonian',
     1           -1,0,0,' ')
      n=0
      nonz=0
      do 100 i=1,ntot
         do 200 j=1,qroots
            if(abs(convec(i,j)).gt.null) then
               n=n+1
               if(n.gt.lenbuf) then
                  nonz=nonz+lenbuf
                  call iosys('write integer '//bufr(1)//
     1                       ' to hamiltonian '//
     2                       'without rewinding',2*lenbuf,ibuf,0,' ')          
                  call iosys('write integer '//bufr(1)//
     1                       ' to hamiltonian '//
     2                       'without rewinding',wptoin(lenbuf),
     3                        hbuf,0,' ')        
                  n=1
               endif
               ibuf(1,n)=i
               ibuf(2,n)=j
               hbuf(n)=convec(i,j)
            endif
 200     continue         
 100  continue
      if (n.gt.0) then
          nonz=nonz+n
          call iosys('write integer '//bufr(1)//
     1               ' to hamiltonian without rewinding',
     2                 2*n,ibuf,0,' ')   
          call iosys('write integer '//bufr(1)//
     1               ' to hamiltonian without rewinding',
     2                 wptoin(n),hbuf,0,' ')
      endif
      call iosys('endfile '//bufr(1)//' on hamiltonian',0,0,0,' ')
      call iosys('write integer '//bufr(2)//' to hamiltonian',
     1            1,nonz,0,' ')
      call iosys('rewind all on hamiltonian read-and-write',0,0,0,' ')
      write(iout,1) nonz
      return
 1    format(/,1x,'number of non-zero contracted pq-space '
     1            'matrix elements = ',i5 )
      end
