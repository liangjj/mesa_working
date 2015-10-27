*deck rdciv.f
c***begin prologue     rdciv
c***date written       980420   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c***keywords           read ci vectors
c***author             schneider, barry (nsf)
c***source             
c***                   
c***references         
c
c***routines called    
c***end prologue       rdciv
      subroutine rdciv(vec,n,nroot)
      implicit integer (a-z)
      real*8 energy, vec
      dimension vec(n,*)
      common/io/inp, iout
      rewind(80)
      do 10 i=1,nroot
         read(80) num, energy
         write(iout,1) num, energy
         read(80) (vec(j,i), j=1,n)
 10   continue
      return
 1    format(/,1x,'root   = ',i4,/,5x,
     1            'energy = ',e15.8)
      end

