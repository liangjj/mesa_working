*deck @(#)readran.f	5.1  11/6/94
      subroutine readran(z,a,data,lenbin,reclen)
c
c***begin prologue     readran
      implicit integer (a-z)
c
      real*8 z
      real*8 pi
      dimension z(*), a(*)
c
      common /io/ inp,iout
      data pi /3.141592653d+00/
      write(iout,*) 'Rereading the File'
      write(iout,*)
      numbin=data/reclen
      left = data - numbin*reclen
      lstbin=reclen
      if(left.gt.0) then
         lstbin=left
         numbin=numbin+1
      endif 
      do 10 bin=1,numbin
         el2wrt=reclen
         if(bin.eq.numbin) then
            el2wrt=lstbin
         endif
         write(iout,*) '      bin number to read = ',bin
         write(iout,*) '      elements to read   = ',el2wrt
         read(unit=99,rec=bin,iostat=ierr) (a(k),k=1,el2wrt)
         write(iout,*) (a(k),k=1,100)
 10   continue   
c
      return
      end
