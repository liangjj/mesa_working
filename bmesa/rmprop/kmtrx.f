*deck kmtrx
c***begin prologue     kmtrx
c***date written       940213   (yymmdd)
c***revision date               (yymmdd)
c***keywords           k-matrix, match
c***author             schneider, barry (nsf)
c***source             
c***purpose            construct k-matrix from r-matrix and
c***                   asymptotic solutions
c***description        
c***references       
c
c***routines called
c***end prologue       kmat
      subroutine kmtrx(rmat,kmat,coef,kmatc,tmat,xsect,y0,dy0,y1,
     1                 dy1,j,jp,y,yp,wron,scr,ipvt,l,ec,energy,
     2                 rstart,nc,nopen,prnt)
      implicit integer (a-z)
      common /io/ inp, iout
      real*8 rmat, kmat, y0, dy0, y1, dy1, ec, energy, rstart, tmp
      real*8 arg, j, jp, y, yp, wron, scr, coef, ekc, sqrekc, xsect
      complex*16 eye, kmatc, tmat
      character*80 title
      logical prnt
      dimension rmat(nc,nc), kmat(nc,nc), ec(nc), y0(nc), dy0(nc)
      dimension y1(nc), dy1(nc), l(nc), coef(nc,nc)
      dimension j(1,0:*), jp(1,0:*), y(1,0:*), yp(1,0:*)
      dimension wron(*), scr(*), ipvt(nc), kmatc(nc,nc)
      dimension tmat(nc,nc), xsect(nc,nc)
      data eye/(0.d0,1.d0)/
      do 10 i=1,nc
         tmp=l(i)*30.d0
         ltop=l(i)+sqrt(tmp)
         ltop=max(ltop,l(i))
         ekc=2.d0*(energy-ec(i))
         if (ekc.ge.0.d0) then
             ekc=sqrt(ekc)
             sqrekc=1.d0/sqrt(ekc)
             arg=ekc*rstart
             call rcbes(arg,j,jp,y,yp,wron,scr,1,l(i),ltop,
     1                  'derivatives',.false.)
             y0(i)=sqrekc*j(1,l(i))
             dy0(i)=sqrekc*ekc*jp(1,l(i))
             y1(i)=-sqrekc*y(1,l(i))     
             dy1(i)=-sqrekc*ekc*yp(1,l(i))
             if (prnt) then
                 write(iout,1) i, ekc, l(i)
                 write(iout,2) y0(i), dy0(i), y1(i), dy1(i)
             endif
         else
             ekc=-ekc
             ekc=sqrt(ekc)
             arg=ekc*rstart
             sqrekc=1.d0/sqrt(ekc)
             call modbes(arg,j,jp,y,yp,wron,scr,1,l(i),ltop,
     1                   'derivatives',.true.,.false.)
             y0(i)=0.d0
             dy0(i)=0.d0
             y1(i)=sqrekc*j(1,l(i))     
             dy1(i)=sqrekc*ekc*jp(1,l(i))
             if (prnt) then
                 write(iout,1) i, ekc, l(i)
                 write(iout,2) y0(i), dy0(i), y1(i), dy1(i)
             endif
         endif             
   10 continue
      call rzero(coef,nc*nc)
      do 20 i=1,nc
         coef(i,i)=coef(i,i)-y1(i)
         do 30 k=1,nc
            coef(k,i)=coef(k,i)+rmat(k,i)*dy1(i)
   30    continue
   20 continue
      nopen=0
      call rzero(kmat,nc*nc)
      do 40 i=1,nc                   
         ekc=2.d0*(energy-ec(i))
         if (ekc.ge.0.d0) then
             nopen=nopen+1
             kmat(i,i)=y0(i)
             do 50 k=1,nc
                kmat(k,i)=kmat(k,i)-rmat(k,i)*dy0(i)
   50        continue
         endif               
   40 continue
      call sgefa(coef,nc,nc,ipvt,info)
      do 60 i=1,nopen
         call sgesl(coef,nc,nc,ipvt,kmat(1,i),0)
   60 continue
      if (prnt) then         
          title='k-matrix'
          call prntrm(title,kmat,nc,nopen,nc,nc,iout)
      endif
      call czero(kmatc,nc*nc)
      do 70 i=1,nopen
         kmatc(i,i)=(1.d0,0.d0)
         do 80 k=1,nopen
            kmatc(i,k)=kmatc(i,k)-eye*kmat(i,k)
            tmat(i,k)=kmat(i,k)
   80    continue
   70 continue
      call cgefa(kmatc,nc,nopen,ipvt,info)
      do 90 i=1,nopen
         call cgesl(kmatc,nc,nopen,ipvt,tmat(1,i),0)
   90 continue
      do 100 i=1,nopen
         ekc=2.d0*(energy-ec(i))
         ekc=1.d0/ekc
         do 110 k=1,nopen
            xsect(i,k)=4.d0*tmat(i,k)*conjg(tmat(i,k))*ekc
  110    continue
  100 continue
      if (prnt) then         
          title='cross sections in pi a0**2'
          call prntrm(title,xsect,nopen,nopen,nc,nc,iout)
      endif                                        
      return
    1 format (/,1x,'channel = ',i3,1x,'k or kappa = ',e15.8,1x,
     1             'l = ',i2)   
    2 format (/,1x,'y0 = ',e15.8,5x,'dy0 = ',e15.8,/,1x,
     1             'y1 = ',e15.8,5x,'dy1 = ',e15.8)
      end
