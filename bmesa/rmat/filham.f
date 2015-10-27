*deck filham
c***begin prologue     filham
c***date written       940213   (yymmdd)
c***revision date               (yymmdd)
c***keywords           hamiltonian, matrix
c***author             schneider, barry (nsf)
c***source             
c***purpose            fill hamiltonian matrix
c***description        
c***references       
c
c***routines called
c***end prologue       filham
      subroutine filham(ham,t,v,ec,nbf,nc,n)
      implicit integer (a-z)
      common /io/ inp, iout
      real*8 ham, t, v, ec
      dimension t(*), v(*), ham(n,n), ec(nc), nbf(nc)
      call rzero(ham,n*n)
      tcnt=0
      vcnt=0
      cntci=0
      do 10 ci=1,nc
         call fildag(ham(cntci+1,cntci+1),t(tcnt+1),ec(ci),nbf(ci),n)
         cntcj=0
         do 20 cj=1,ci
            call filodg(ham(cntci+1,cntcj+1),ham(cntcj+1,cntci+1),
     1                  v(vcnt+1),ci,cj,nbf(ci),nbf(cj),n)
            vcnt=vcnt+nbf(ci)*nbf(cj)
            cntcj=cntcj+nbf(cj)
   20    continue
         tcnt=tcnt+nbf(ci)*nbf(cj)
         cntci=cntci+nbf(ci)               
   10 continue
      return
      end
