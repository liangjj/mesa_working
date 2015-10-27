*deck @(#)iomesa.f	5.1  11/6/94
      subroutine iomesa(z,a,data,lenbin,type)
c
c***begin prologue     iomesa
      implicit integer (a-z)
c
      logical logkey,prnt
      real*8 z(*)
      real*8 pi
      character*4 itoc
      character*(*) type
      dimension a(*) 
c
      common /io/ inp,iout
      data pi /3.141592653d+00/
      data prnt /.true./
c
c
c
      dtapp=lenbin
      npass=data/dtapp
      left=data-npass*dtapp
      dtlast=dtapp
      if(left.gt.0) then
         dtlast=left
         npass=npass+1
      endif
      if(type.eq.'integer') then
         call izero(a,dtapp)
         do 10 i=1,dtapp
            a(i)=i
 10      continue   
      elseif(type.eq.'real') then
         call rzero(z,dtapp)
         do 20 i=1,dtapp
            z(i)=pi/i
 20      continue
      else
         call lnkerr('error in data type')
      endif   
      write(iout,1) npass, data
      total=0
      do 30 i=1,npass
         dtowr=dtapp
         if(i.eq.npass) then
            dtowr=dtlast
         endif
         write(iout,2) i, dtowr
         call iosys('write '//type//' array on testunit without '//
     #                 'rewinding',dtowr,a,0,' ')
         if(type.eq.'integer') then 
            write(iout,3) ( a(j),j=1,10)
         else
            write(iout,4) ( z(j),j=1,10)
         endif
         total=total+dtowr
         write(iout,5) total
 30   continue   
      write(iout,6) total
      call iosys('rewind all on testunit',0,0,0,' ')
      total=0
      do 40 i=1,npass
         dtowr=dtapp
         if(i.eq.npass) then
            dtowr=dtlast
         endif
         write(iout,7) i, dtowr
         call iosys('read '//type//' array from testunit without '//
     #              'rewinding',dtowr,a,0,' ')
         if(type.eq.'integer') then
            write(iout,3) ( a(j),j=1,10)
         else
            write(iout,4) ( z(j),j=1,10)
         endif
         total=total+dtowr
         write(iout,8) total
 40   continue   
      write(iout,9) total
      return
 1    format(/,10x,'iomesa:'
     #       /,1x,'total number of passes        = ',i10,
     #       /,1x,'total number of data elements = ',i10,//)
 2    format(/,10x,'iomesa:'
     #       //,1x,'pass = ',i10,2x,'elements to write = ',i10)
 3    format(/,10x,'iomesa:'
     #       /,1x,'first 10 elements',
     #       (/,1x,5i15))
 4    format(/,10x,'iomesa:'
     #       /,1x,'first 10 elements',
     #       (/,1x,5e15.8))
 5    format(/,10x,'iomesa:'
     #       /,1x,'elements written = ',i10)
 6    format(/,10x,'iomesa:','final total number of elements written'
     #                       ' = ',i10)
 7    format(/,10x,'iomesa:'
     #       //,1x,'pass = ',i10,2x,'elements to read = ',i10)
 8    format(/,10x,'iomesa:'
     #       /,1x,'elements read = ',i10)
 9    format(/,10x,'iomesa:','final total number of elements read'
     #                       ' = ',i10)
c
      end
