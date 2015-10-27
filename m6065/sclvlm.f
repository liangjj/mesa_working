*deck sclvlm
c***begin prologue     sclvlm
c***date written       890512   (yymmdd)
c***revision date               (yymmdd)
c***keywords
c***author             schneider, barry (nsf)
c***source             m6060
c***purpose            multiply complex vlamdas by integration weights
c***references
c
c***routines called   none
c***end prologue      sclvlm
      subroutine sclvlm(vlamda,grid,npnts,nolam)
      implicit integer (a-z)
      real *8 grid, sq2
      complex *16 vlamda
      dimension vlamda(npnts,nolam), grid(4,npnts)
      sq2=sqrt(2.d0)
      do 10 i=1,nolam
         do 20 grpt=1,npnts
            vlamda(grpt,i)=sq2*vlamda(grpt,i)*grid(4,grpt)
   20    continue
   10 continue
      return
      end

