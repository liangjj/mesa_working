*deck @(#)rhotst.f	1.1 9/8/91
c***begin prologue     rhotst
c***date written       920614   (yymmdd)
c***revision date               (yymmdd)
c***keywords           density
c***author             schneider, barry(nsf)
c***source             @(#)m7000
c***purpose            integrate transition densities
c***
c***
c***references
c
c***routines called    util
c***end prologue
      subroutine rhotst(rhor,typa,rhoc,typb,valr,valc,npr,n)
      implicit integer (a-z)
      real *8 rhor, valr
      complex *16 rhoc, valc
      character *(*) typa, typb
      character *128 filnme
      dimension rhor(n), rhoc(n), valr(npr), valc(npr)
      common /io/ inp, iout
      filnme='"'//typa//'-'//typb//' densities"'
      if (typa.eq.'bound'.and.typb.eq.'bound') then
          call rzero(valr,npr)
          do 10 i=1,npr
             call iosys ('read real '//filnme//' from '//
     1                   'atomci without rewinding',n,
     2                    rhor,0,' ')
             do 20 j=1,n
                valr(i)=valr(i)+rhor(j)
   20        continue
   10     continue
          write(iout,100) ( valr(i),i=1,npr)
      else
          call czero(valc,npr)
          do 30 i=1,npr
             call iosys ('read real '//filnme//' from '//
     1                   'atomci without rewinding',2*n,
     2                    rhor,0,' ')
             do 40 j=1,n
                valc(i)=valc(i)+rhoc(j)
   40        continue
   30     continue
          write(iout,200) ( valc(i),i=1,npr)
      endif
      return
  100 format((5x,4e15.8))
  200 format(5x,2('(',e15.8,','e15.8,')'))
      end


