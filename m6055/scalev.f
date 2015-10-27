*deck @(#)scalev.f	1.1 9/8/91
c***begin prologue     scalev
c***date written       890512   (yymmdd)
c***revision date               (yymmdd)
c***keywords
c***author             schneider, barry (lanl)
c***source             m6005
c***purpose            multiply potential by integration weights
c***                   handles complex potentials
c***references
c
c***routines called   none
c***end prologue      scalev
      subroutine scalev(v,grid,npnts,nstri)
      implicit integer (a-z)
      real *8 grid, v
      dimension v(npnts,nstri), grid(4,npnts)
      do 10 i=1,nstri
         do 20 grpt=1,npnts
            v(grpt,i)=v(grpt,i)*grid(4,grpt)
   20    continue
   10 continue
      return
      end
