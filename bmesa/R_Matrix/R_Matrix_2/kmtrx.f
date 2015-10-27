*deck kmtrx.f
c***begin prologue     kmtrx
c***date written       940213   (yymmdd)
c***revision date               (yymmdd)
c***keywords           k-matrix
c***author             schneider, barry (nsf)
c***source             
c***purpose            
c***description        calculate the k-matrix from the r-matrix
c***                   and the linearly independent solutions  
c***references       
c
c***routines called
c***end prologue       kmtrx
      subroutine kmtrx(rmat,sn,dsn,cn,dcn,echn,typ,energy,
     1                 kmat,rc,work,ipvt,nchan,nopen,lwork,prn)
      implicit integer (a-z)
      common /io/ inp, iout
      real*8 rmat, sn, dsn, cn, dcn, energy, echn, kmat, rc, work
      character*4 typ
      character*16 fptoc
      character*80 title
      logical prn
      dimension rmat(nchan,nchan)
      dimension sn(nchan), dsn(nchan), cn(nchan), dcn(nchan)
      dimension echn(nchan), typ(nchan)
      dimension kmat(nchan,nopen), rc(nchan,nchan)
      dimension work(*), ipvt(nchan)
      call rzero(kmat,nchan*nopen)
      call rzero(rc,nchan*nchan)
c
c     the matrix is filled so that the open channels come first.
c
      nclosed=nchan-nopen
      open=0
      closed=nopen
      do 10 chn=1,nchan
         if(typ(chn).eq.'o') then
            open=open+1
            call vscale(rc(1,open),rmat(1,chn),dcn(chn),nchan)
            call vneg(rc(1,open),rc(1,open),nchan)
            rc(open,open)=cn(chn)+rc(open,open)
            call vscale(kmat(1,open),rmat(1,chn),dsn(chn),nchan)
            kmat(open,open)=kmat(open,open)-sn(chn)
         else
            closed=closed+1
            call vscale(rc(1,closed),rmat(1,chn),dcn(chn),nchan)
            call vneg(rc(1,closed),rc(1,closed),nchan)
            rc(closed,closed)=cn(chn)+rc(closed,closed)
         endif
 10   continue
      if(prn) then
         title='coefficient matrix'   
         call prntrm(title,rc,nchan,nchan,nchan,nchan,iout)
         title='right hand side'   
         call prntrm(title,kmat,nchan,nopen,nchan,nopen,iout)
      endif
c
c     solve for K-matrix
c
      call dsysv('u',nchan,nopen,rc,nchan,ipvt,kmat,nchan,
     1               work,lwork,info)
      if(info.ne.0) then
         call lnkerr('error in solving for k-matrix')
      endif
      if(prn) then
         title='K(o,o)-matrix:r energy = '//fptoc(energy)
         call prntrm(title,kmat,nopen,nopen,nchan,nopen,iout)
         if(nclosed.ne.0) then
            title='K(on,c)-matrix: energy = '//fptoc(energy)
            call prntrm(title,kmat(1,nopen+1),nopen,nclosed,
     1                  nchan,nopen,iout)
            title='K(c,o)-matrix: energy = '//fptoc(energy)
            call prntrm(title,kmat(nopen+1,1),nclosed,nopen,
     1                  nchan,nopen,iout)
            title='K(c,c)-matrix: energy = '//fptoc(energy)
            call prntrm(title,kmat(nopen+1,nopen+1),nclosed,nclosed,
     1                  nchan,nopen,iout)
         endif 
      endif
      return
      end
