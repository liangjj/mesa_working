*deck @(#)tstcon.f	5.1  11/6/94
      subroutine tstcon(b,grad,rmsd,mdim,nrhs,stopj,icnvrg)
      implicit real*8(a-h,o-z)
      real*8 b(mdim,nrhs),grad(mdim,nrhs),rmsd(nrhs)
      common /io/ inp,iout
c
      call rzero(rmsd,nrhs)
c
      testj=50.d0*stopj
c
      icnvrg=1
c
      do 16 j=1,nrhs
         do 15 k=1,mdim
            rmsd(j)=rmsd(j)+(b(k,j)-grad(k,j))**2
 15      continue
c
         rmsd(j)=sqrt(rmsd(j))/float(mdim)
c
         if(rmsd(j).gt.testj) then
            icnvrg=0
         endif
 16   continue
c
      write(iout,20)(rmsd(j),j=1,nrhs)
  20  format(/,'  rmsd ',/,5(2x,e15.7))
c
      return
      end
