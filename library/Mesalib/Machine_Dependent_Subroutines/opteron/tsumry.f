*deck @(#)tsumry.f	5.1  11/6/94
      subroutine tsumry      
      implicit integer(a-z)
c
      real*8 cpuchg,iochg,syschg,zero,totcpu,totio,totsys
c
      parameter (zero=0.0d+00)
      parameter (mxlnk=200)
      integer nenter(mxlnk),lnknos(mxlnk),coruse(mxlnk)
      real*8 times(3,mxlnk)
      common/io/inp,iout
c
 1000 format(1x,'summary:')
 1010 format(5x,'link   entries               charges          core use
     $      ')
 1020 format(5x,'                     cpu       sys      io         ')
 1050 format(5x,i4,5x,i5,3(2x,f8.1),i10)
 1060 format(13x,'total:',3(2x,f8.1),i10)
c
c
c     read in the common blocks which oversee the management of
c     mesa.
c
c        'links' contains global information regarding link timings
c        and memory usage.
c        nenter is the number of times a link has been executed,
c        lnknos contains the link number, coruse the memory requested,
c        and times(3,i) contains the cpu,system,and (unused) charges,
c        respectively.
c
c
c
      call iosys('read integer "links:nenter" from rwf',
     $           -1,nenter,0,' ')
      call iosys('read integer "links:lnknos" from rwf',
     $           -1,lnknos,0,' ')
      call iosys('read integer "links:coruse" from rwf',
     $           -1,coruse,0,' ')
      call iosys('read real "links:times" from rwf',
     $           -1,times,0,' ')
c
c     now process it.
      write(iout,1000)
      write(iout,1010)
      write(iout,1020)
      totcpu=zero
      totsys=zero
      totio=zero
      maxmem=0
      do 10 i=1,mxlnk
         link=lnknos(i)
         if(link.ne.0) then
            ntimes=nenter(i)
            cpuchg=times(1,i)
            syschg=times(2,i)
            iochg=times(3,i)
            cormax=coruse(i)
            totcpu=totcpu+cpuchg
            totsys=totsys+syschg
            totio=totio+iochg
            memmax=max(memmax,cormax)
            write(iout,1050) link,ntimes,cpuchg,syschg,iochg,cormax
         endif
   10 continue
c
      write(iout,1060) totcpu,totsys,totio,memmax
c
      return
      end
