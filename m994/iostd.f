*deck @(#)iostd.f	5.1  11/6/94
      subroutine iostd(z,a,data,lenbin,reclen)
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
c     reclen is in 4 byte words 
c
c     
c     open a direct access file having reclen as its default record length.
c
      open(unit=99,access='direct',status='scratch',
     #     recl=reclen,iostat=iostat)
c
c
c
c     given the data, the bin size and record length, determine how many passes
c     are needed to write the data.
c   
      numbin=data/lenbin
      left = data - numbin*lenbin
      lstbin=lenbin
      if(left.gt.0) then
         lstbin=left
         numbin=numbin+1
      endif 
      count=0
      total=0
      recno=0      
      write(iout,*) 'Test of Standard Write and Read Routines'
      write(iout,*) 'Beginning the Write Phase'
      write(iout,*) 'Write ',data,' 4 byte words'
      write(iout,*)
      write(iout,*) 'Number of Bins Used = ',numbin
      write(iout,*)
      do 10 bin=1,numbin
         el2wrt=lenbin
         if(bin.eq.numbin) then
            el2wrt=lstbin
         endif
         do 20 i=1,el2wrt
            count=count+1
            a(i)=count
 20      continue   
         write(iout,*) '      bin number to write = ',bin
         write(iout,*) '      elements to write   = ',el2wrt
         npass=el2wrt/reclen
         lastwd=reclen
         left = el2wrt - npass*reclen      
         if(left.gt.0) then
            npass=npass+1
            lastwd=left
         endif
         write(iout,*) '      number of passes    = ',npass 
         write(iout,*) '      number of last pass = ',lastwd
         first=1
         last=reclen
         do 30 j=1,npass-1
            recno=recno+1
            write(unit=99,rec=recno,iostat=ierr) (a(k),k=first,last)
            if(ierr.ne.0) then
               write(iout,*) '     record number = ',recno
               write(iout,*) '     first         =',first
               write(iout,*) '     last          =',last
c               write(iout,*) (a(k),k=first,first+9)
c               write(iout,*) (a(k),k=last-9,last)
               write(iout,*)
               write(iout,*) 'ERROR = ',ierr,'CONTINUE'
            endif
            first=last+1
            last=last+reclen
            total=total+reclen 
 30      continue   
         write(iout,*) 'total written = ',total
         recno=recno+1
         last = last - reclen + left
         write(unit=99,rec=recno,iostat=ierr) (a(k),k=first,last)
         if(ierr.ne.0) then
            write(iout,*) '     record number = ',recno
            write(iout,*) '     first         =',first
            write(iout,*) '     last          =',last
c            write(iout,*) (a(k),k=first,first+9)
c            write(iout,*) (a(k),k=last-9,last)
         endif
         total = total + last - first + 1
 10   continue   
      write(iout,*) 'final total written = ',total
c
c     read in back
c
      write(iout,*) 'Beginning the Read Phase'
      total=0
      recno=0
      do 40 bin=1,numbin
         el2wrt=lenbin
         if(bin.eq.numbin) then
            el2wrt=lstbin
         endif
         write(iout,*) '      bin number to read = ',bin
         write(iout,*) '      elements to read   = ',el2wrt
         npass=el2wrt/reclen
         write(iout,*) '      number of passes    = ',npass 
         lastwd=reclen
         left = el2wrt - npass*reclen      
         if(left.gt.0) then
            npass=npass+1
            lastwd=left
         endif
         first=1
         last=reclen
         do 50 j=1,npass-1
            recno=recno+1
            read(unit=99,rec=recno,iostat=ierr) (a(k),k=first,last)
            if(ierr.ne.0) then
               write(iout,*) '     record number = ',recno
               write(iout,*) '     first         =',first
               write(iout,*) '     last          =',last
c               write(iout,*) (a(k),k=first,first+9)
c               write(iout,*) (a(k),k=last-9,last)
               write(iout,*)
               write(iout,*) 'ERROR = ',ierr,'CONTINUE'
            endif
            first=last+1
            last=last+reclen
            total=total+reclen 
 50      continue   
         write(iout,*) 'total read = ',total
         recno=recno+1
         last = last - reclen + left
         read(unit=99,rec=recno,iostat=ierr) (a(k),k=first,last)
         if(ierr.ne.0) then
            write(iout,*) '     record number = ',recno
            write(iout,*) '     first         =',first
            write(iout,*) '     last          =',last
c            write(iout,*) (a(k),k=first,first+9)
c            write(iout,*) (a(k),k=last-9,last)
         endif
         total = total + last - first + 1
 40   continue   
      write(iout,*) 'final total read = ',total
c
      return
      end
