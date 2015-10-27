*deck tranpp.f
c***begin prologue     tranpp
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
c***end prologue       tranpp
      subroutine tranpp(ibuf,hbuf,diag,pvec,convec,scr,lenbuf,
     1                  nwks,nwksp,ntot,headr,prnt,dsk)
c
      implicit integer (a-z)
c
      real*8 hbuf, diag, pvec, convec, scr, null
      character*(*) headr, dsk
      character*80 title
      character*40 bufr
      character*3 answer
      logical nodsk, prnt
      dimension ibuf(2,lenbuf), hbuf(lenbuf), diag(nwksp)
      dimension pvec(nwksp,*), convec(ntot,*), headr(3)
      dimension scr(nwksp,*), bufr(3)
      common /io/ inp,iout
      data nodsk/.false./
      data null/1.d-14/
c
      call iosys('does "p-space vectors" exist on '//dsk,0,0,0,answer)
      offset=nwks-nwksp+1
      if(answer.eq.'yes') then
         do 10 i=1,ntot
            call iosys('read real "p-space vectors" from '//dsk//
     1                 ' without rewinding',nwks,scr,0,' ')
            call copy(scr(offset,1),pvec(1,i),nwksp)
 10      continue   
      else
         call rzero(pvec,ntot*nwksp)
         do 20 i=1,ntot
            pvec(i,i)=1.d0
 20      continue
      endif   
c
      title='new vector set'
      call prntrm(title,pvec,nwksp,ntot,nwksp,ntot,iout)
      call rzero(scr,ntot*nwksp)
      call iosys('rewind all on hamiltonian',0,0,0,' ')
      call iosys('read integer '//headr(2)//' from hamiltonian',
     1            1,nonz,0,' ')
      call iosys('read real '//headr(3)//' from hamiltonian',nwksp,
     1            diag,0,' ')  
      call hmult(ibuf,hbuf,lenbuf,headr(1),nonz,pvec,scr,nwksp,
     1           ntot,nodsk)
      do 30 i=1,nwksp
         do 40 j=1,ntot
            scr(i,j) = scr(i,j) + diag(i)*pvec(i,j)
 40      continue
 30   continue      
      call ebtc(convec,pvec,scr,ntot,nwksp,ntot)
      if(prnt) then
         title='transformed contracted p-space hamiltonian'
         call prntrm(title,convec,ntot,ntot,ntot,ntot,iout)
      endif
      do 50 i=1,ntot
         diag(i) = convec(i,i)
 50   continue   
      bufr(1)='"contracted p space buffers"'
      bufr(2)='"number of contracted p space elements"'
      bufr(3)='"contracted p space diagonals"'
      call iosys('write real '//bufr(3)//' to hamiltonian', 
     1            ntot,diag,0,' ')
      call iosys('create integer '//bufr(1)//' on hamiltonian',
     1           -1,0,0,' ')
      n=0
      nonzp=0
      do 60 i=1,ntot
         do 70 j=1,i-1
            if(abs(convec(i,j)).gt.null) then
               n=n+1
               if(n.gt.lenbuf) then
                  nonzp=nonzp+lenbuf
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
 70      continue         
 60   continue
      if (n.gt.0) then
          nonzp=nonzp+n
          call iosys('write integer '//bufr(1)//
     1               ' to hamiltonian without rewinding',
     2                 2*n,ibuf,0,' ')   
          call iosys('write integer '//bufr(1)//
     1               ' to hamiltonian without rewinding',
     2                 wptoin(n),hbuf,0,' ')
      endif
      call iosys('endfile '//bufr(1)//' on hamiltonian',0,0,0,' ')
      call iosys('write integer '//bufr(2)//' to hamiltonian',
     1            1,nonzp,0,' ')
      call iosys('rewind all on hamiltonian read-and-write',0,0,0,' ')
      write(iout,1) nonzp
      return
 1    format(/,1x,'number of non-zero contracted p-space '
     1            'matrix elements = ',i5 )
      end




