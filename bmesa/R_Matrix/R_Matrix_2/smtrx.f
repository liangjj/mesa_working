*deck smtrx.f
c***begin prologue     smtrx
c***date written       940213   (yymmdd)
c***revision date               (yymmdd)
c***keywords           s-matrix
c***author             schneider, barry (nsf)
c***source             
c***purpose            
c***description        calculate the s-matrix from the k-matrix
c***references       
c
c***routines called
c***end prologue       smtrx
      subroutine smtrx(smat,kmat,energy,coef,work,ipvt,
     1                 nchan,nopen,lwork,prn)
      implicit integer (a-z)
      common /io/ inp, iout
      real*8 kmat, energy
      complex*16 smat, eye, coef, work
      character*16 fptoc
      character*80 title
      dimension smat(nopen,nopen), kmat(nchan,nopen), coef(nopen,nopen)
      dimension work(*), ipvt(nopen)
      data eye/(0.d0,1.d0)/
      call czero(coef,nopen*nopen)
      call czero(smat,nopen*nopen)
      do 10 chni=1,nopen
         coef(chni,chni)=1.d0
         smat(chni,chni)=1.d0
 10   continue   
      do 20 chni=1,nopen
         do 30 chnj=1,nopen
            coef(chni,chnj) = coef(chni,chnj) - eye*kmat(chni,chnj)
            smat(chni,chnj) = smat(chni,chnj) + eye*kmat(chni,chnj)
 30      continue
 20   continue   
c 
      if(prn) then
         title='coefficient matrix'   
         call prntcm(title,coef,nopen,nopen,nopen,nopen,iout)
         title='right hand side'   
         call prntcm(title,smat,nopen,nopen,nopen,nopen,iout)
      endif
c
c     solve for S-matrix
c
      call zgesv(nopen,nopen,coef,nopen,ipvt,smat,nopen,info)
      if(info.ne.0) then
         call lnkerr('error in solving for s-matrix')
      endif
      if(prn) then
         title='S-matrix: energy = '//fptoc(energy)
         call prntcm(title,smat,nopen,nopen,nopen,nopen,iout)
      endif
      return
      end
