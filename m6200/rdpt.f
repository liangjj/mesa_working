*deck rdpt.f
c***begin prologue     rdpt
c***date written       930502   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c***keywords           rdpt, link m6200
c***author             schneider, barry (nsf)
c***source             m6200
c***purpose            read in radial points
c***references         
c
c***routines called
c***end prologue       rdpt
      subroutine rdpt (rpt,nshell,nr,numshl,nrmax,str,iat)
      implicit integer (a-z)
      real*8 rpt
      character*(*) str
      dimension rpt(nrmax,numshl), nr(numshl)
      common /io/ inp, iout
      call iosys ('rewind "radial points '//str//
     1            '" on lamdat read-and-write',0,0,0,' ')  
      do 10 i=1,nshell
         call iosys ('read real "radial points '//str//
     1               '" from lamdat without rewinding',nr(i),
     2                  rpt(1,i),0,' ')
   10 continue        
      return
      end


