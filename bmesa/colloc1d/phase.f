deck phase.f
c***begin prologue     phase
c***date written       910128   (yymmdd)
c***revision date               (yymmdd)
c***keywords           
c***author             schneider, barry(lanl)
c***source             @(#)m6020
c***purpose            calculate phase shift and display solution
c***
c***
c***references
c
c***routines called    util
c***end prologue
      subroutine phase(psi,r,tvec,energy,y,n,prnt,nobc)
      implicit integer (a-z)
      dimension psi(n), tvec(n), r(n), dmat(n,n)
      real*8 psi, tvec, energy, y, amplt, pshft, k, r
      character*80 title
      logical prnt, nobc
      common /io/ inp, iout
      if(.not.nobc) then
         psi(n)=0.d0
         do 10 i=1,n-1
            psi(n)=psi(n)+tvec(i)*psi(i)
   10    continue
      endif         
      if (prnt) then
          title='             r             solution vector'
          write(iout,1) title
          do 20 i=1,n
             write(iout,2) r(i), psi(i)
   20     continue   
      endif
      amplt=psi(n)/y
      write(iout,*) '     the ampltude of the irregular solution = '
     1                 ,amplt
      pshft=atan(amplt)
      write(iout,*) '     the phase shift = ',-pshft
      return
    1 format(a80)
    2 format(5x,e15.8,5x,e15.8)   
      end
