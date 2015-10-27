*deck @(#)iostd.f	5.1  11/6/94
      subroutine iostd(z,a,data,reclen)
c
c***begin prologue     iostd
      implicit integer (a-z)
c
      real*8 z
      real*8 pi
      dimension z(*), a(*)
c
      common /io/ inp,iout
      data pi /3.141592653d+00/
c
c
c
c
c     reclen is in 4 byte words make into 8 byte words
c
      dtapp=reclen
      reclen=2*reclen 
c     
c     open a direct access file having reclen as its default record length.
c
      open(unit=99,access='direct',status='scratch',
     #     recl=reclen,iostat=iostat)
c
c      just to put something in to z
c
      do 10 i=1,dtapp
         z(i)=pi/i
 10      continue
c
c     given the data and the record length, dtermine how many passes
c     are needed to write the data.
c   
      npass=data/dtapp
      left=data-npass*dtapp      
      dtlast=dtapp
      if(left.gt.0) then
         npass=npass+1
         dtlast=left
      endif
      write(iout,1) npass, data
      total=0
      do 20 i=1,npass
         dtowr=dtapp
         if(i.eq.npass) then
            dtowr=dtlast
         endif
         write(iout,2) i, dtowr
         write(unit=99,rec=i,iostat=ierr) (z(j),j=1,dtowr)
         if(ierr.ne.0) then
            write(iout,3) dtowr
            write(iout,4)
         endif
         write(iout,3) dtowr
         total=total+dtowr
         write(iout,5) total
 20   continue   
      rewind 99
      total=0
      do 30 i=1,npass
         dtowr=dtapp
         if(i.eq.npass) then
            dtowr=dtlast
         endif
         write(iout,6) i, dtowr
         read(unit=99,rec=i,iostat=ierr) (z(j),j=1,dtowr)
         if(ierr.ne.0) then
            write(iout,7) dtowr
            write(iout,8)
            write(iout,9) (z(j),j=1,100)
         endif
         write(iout,7) dtowr
         total=total+dtowr
         write(iout,5) total
 30   continue 
c
      return
 1    format(/,10x,'iostd:'
     #       /,1x,'total number of passes        = ',i10,
     #       /,1x,'total number of data elements = ',i10,//)
 2    format(/,10x,'iostd:'
     #       //,1x,'pass = ',i10,2x,'elements to write = ',i10)
 3    format(/,10x,'iostd:'
     #       /,1x,'elements written = ',i10)
 4    format(/,1x,'write error')
 5    format(/,10x,'iostd:','final total number of elements = ',i10)
 6    format(/,10x,'iostd:'
     #       //,1x,'pass = ',i10,2x,'elements to read = ',i10)
 7    format(/,10x,'iostd:'
     #       /,1x,'elements read = ',i10)
 8    format(/,1x,'read error')
 9    format(5e15.8)
      end
