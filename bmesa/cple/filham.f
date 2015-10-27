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
      subroutine filham(ham,t,v,ec,ni,nj,nc,nclosd,n)
      implicit integer (a-z)
      common /io/ inp, iout
      real*8 ham, t, v, ec
      dimension t(*), v(*), ham(n,n), ec(nc), ni(nc), nj(nc)
      call rzero(ham,n*n)
      tcnt=0
      vcnt=0
      cntcic=0
      cntcio=nclosd
      do 10 ci=1,nc
         mi=ni(ci)+nj(ci) 
         call fildag(ham(cntcic+1,cntcic+1),ham(cntcic+1,cntcio+1),
     1               ham(cntcio+1,cntcic+1),ham(cntcio+1,cntcio+1),
     2               t(tcnt+1),ec(ci),ni(ci),nj(ci),mi,n)
         cntcjc=0
         cntcjo=nclosd
         do 20 cj=1,ci
            mj=ni(cj)+nj(cj)
            call filodg(ham(cntcic+1,cntcjc+1),ham(cntcjc+1,cntcic+1),
     1                  ham(cntcic+1,cntcjo+1),ham(cntcjc+1,cntcio+1),
     2                  ham(cntcio+1,cntcjc+1),ham(cntcjo+1,cntcic+1),
     3                  ham(cntcio+1,cntcjo+1),ham(cntcjo+1,cntcio+1),
     4                  v(vcnt+1),ci,cj,ni(ci),nj(ci),ni(cj),nj(cj),
     5                  mi,mj,n)
            cntcjc=cntcjc+ni(cj)
            cntcjo=cntcjo+nj(cj)
            vcnt=vcnt+mi*mj
   20    continue
         cntcic=cntcic+ni(ci)
         cntcio=cntcio+nj(ci)
         tcnt=tcnt+mi*mi               
   10 continue
      return
      end
