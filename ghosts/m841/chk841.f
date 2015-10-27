*deck @(#)chk841.f	1.1  11/30/90
      subroutine chk841(z,a,maxcor)
c
      implicit integer(a-z)
c
      real*8 z(maxcor),en2,sdot
c
      common /io/ inp,iout
c
      call iosys('read integer "number of grouped 2pdm elements" '//
     $     'from rwf',1,n2pdm,0,' ')
c
      call iosys('rewind "group ordered 2e-ints" on rwf',0,0,0,' ')
      call iosys('rewind "group ordered 2pdm" on rwf',0,0,0,' ')
c
      maxc=min(n2pdm,(maxcor-1)/2)
c
      dm=1
      ints=dm+maxc
c
      en2=0.d0
c
      do 1 i=1,n2pdm,maxc
c
       lnread=min(maxc,n2pdm-i+1)
       call iosys('read real "group ordered 2pdm" from rwf '
     #                    //'without rewinding',lnread,z(dm),0,' ')
       call iosys('read real "group ordered 2e-ints" from rwf '
     #                    //'without rewinding',lnread,z(ints),0,' ')
c
      en2=en2+sdot(lnread,z(dm),1,z(ints),1)
c
  1   continue
c
      write(iout,2) en2
 2    format(/,5x,' m841: two-electron energy ',g20.12)
c
      return
      end
