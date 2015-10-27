*deck @(#)pm991.f	1.3  7/30/91
      subroutine pm991
c***begin prologue     m991
c***date written       020731   (yymmdd)
c***revision date               (yymmdd)
c
c***keywords           m9911, link 991, io test
c***author             Schneider, Barry (NSF)
c***source             @(#)m911.f
c***purpose            test i/o under new compaq compiler/linux
c***description
c
c***references
c
c***routines called
c***end prologue       pm811
c
      implicit integer (a-z)
c
      character*4096 ops
      character*8 typdat, filtyp
      character*128 chrkey, name
      logical logkey
      integer need, ngot
      integer a
      real*8 z
      pointer(p,a(1)), (p,z(1))
c
      common /io/     inp,iout
c
c
 1000 format(1x,'m991:io test',//)
 1001 format(5x,'memory use        ',11x,i9)
cdir$ fastmd
c
c     ----- recover the options string -----
c
      write(iout,1000)      
      ops=' '
      call iosys('read character options from rwf',-1,0,0,ops)
      name=chrkey(ops,'unit-name','nameun',' ')
      filtyp=chrkey(ops,'file-type','new',' ')
      typdat=chrkey(ops,'data-type','integer',' ')
      len=length(typdat)
      data=intkey(ops,'number-of-data-elements',1000000,' ')
      lenbin=intkey(ops,'bin-size',100000,' ')
      lenbin=min(lenbin,data)
      call getmem(0,p,ngot,'m991:begin',0)
      call iosys('read integer mxcore from rwf',1,maxcor,0,' ')
      val=1
      if(typdat(1:len).eq.'real') then
         need=wpadti(val+lenbin)
      else
         need=val+lenbin
      endif
      write(iout,*) '    available core              = ',maxcor
      write(iout,*) '    length of bin               = ',lenbin
      write(iout,*) '    no. data elements to write  = ',data
      write(iout,*) '    memory used                 = ',need
      call getmem(need,p,ngot,'m991',0)
      write(iout,*) '    creating a '//typdat(1:len)//' file'
      len1=length(name)
      write(iout,*) '    file name = ',name(1:len1)
      write(iout,*) '    file type = ',filtyp
      call iosys('open testunit as '//filtyp,0,0,0,name)
      call iosys('create '//typdat(1:len)//' array on testunit',
     #            data,0,0,' ')
      call iomesa(z(val),a(val),data,lenbin,typdat(1:len))
      call getmem(-ngot,p,idum,'m991',idum)
c
c     ----- and exit with grace -----
c
      call chainx(0)
c
c
      stop
      end
