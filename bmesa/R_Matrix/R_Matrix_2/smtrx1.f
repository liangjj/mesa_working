*deck smtrx1.f
c***begin prologue     smtrx1
c***date written       940213   (yymmdd)
c***revision date               (yymmdd)
c***keywords           s-matrix
c***author             schneider, barry (nsf)
c***source             
c***purpose            
c***description        calculate the s-matrix from the r-matrix for 1 channel
c***references       
c
c***routines called
c***end prologue       smtrx1
      subroutine smtrx1(smat,rmat,energy,rbox)
      implicit integer (a-z)
      common /io/ inp, iout
      real*8 rmat, energy, rbox, k
      complex*16 smat, eye
      character*16 fptoc
      character*80 title
      data eye/(0.d0,1.d0)/
      k=sqrt(2.d0*energy)
      smat= (1.d0 + eye*k*rmat)*exp(-2.d0*eye*k*rbox)
     1            /( 1.d0 - eye*k*rmat)
c 
      title='S-matrix: energy = '//fptoc(energy)
      call prntcm(title,smat,1,1,1,1,iout)
      return
      end
