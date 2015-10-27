*deck @(#)prntv.f	1.1 9/7/91
c***begin prologue     prntv
c***date written       890412   (yymmdd)
c***revision date               (yymmdd)
c***keywords           prntv, link 6003
c***authors            schneider, barry (lanl)
c***                   
c***source             m6002
c***purpose            print static potential
c***references       
c
c***routines called    
c***end prologue       prntv
      subroutine prntv(pot,grid,nstri,npt,pass)
      implicit integer (a-z)
      real *8 pot, grid
      character *80 title
      character *3 itoc
      dimension pot(npt,nstri), grid(4,npt)
      common /io/ inp, iout
      write (iout,10) pass
      title='static potential'
      do 20 i=1,nstri
         title='triangle label-'//itoc(i)
         write (iout,30) title
         write (iout,40)
         do 50 j=1,npt
            write (iout,60) grid(1,j), grid(2,j), grid(3,j), pot(j,i)
   50    continue
   20 continue
      return
   10 format(/,5x,'potential matrix for pass',1x,i3)
   30 format(5x,a80)
   40 format(/,11x,'x',13x,'y',13x,'z',12x,'potential')
   60 format(3x,f12.6,2x,f12.6,2x,f12.6,2x,f15.8)
      end
