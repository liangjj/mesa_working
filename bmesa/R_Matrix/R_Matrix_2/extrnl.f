*deck extrnl.f
c***begin prologue     extrnl
c***date written       940213   (yymmdd)
c***revision date               (yymmdd)
c***keywords           external solutions
c***author             schneider, barry (nsf)
c***source             
c***purpose            
c***description        calculate the linearly independent external solutions  
c***references       
c
c***routines called
c***end prologue       extrnl
      subroutine extrnl(sn,dsn,cn,dcn,eigc,echn,typ,rbox,energy,
     1                  nc,nopen,prn)
      implicit integer (a-z)
      common /io/ inp, iout
      real*8 sn, dsn, cn, dcn, eigc, echn, rbox, energy, ek
      character*4 typ
      character*80 title
      logical prn
      dimension sn(nc), dsn(nc), cn(nc), dcn(nc)
      dimension eigc(nc), echn(nc), typ(nc)
      call rzero(sn,nc)
      call rzero(dsn,nc)
      call rzero(cn,nc)
      call rzero(dcn,nc)
      nopen=0
      do 10 chn=1,nc
         echn(chn)=energy-eigc(chn)
         if(echn(chn).gt.0.d0) then
            typ(chn)='o'
            nopen=nopen+1
            ek=sqrt(2.d0*echn(chn))
            sn(chn)=sin(ek*rbox)
            cn(chn)=cos(ek*rbox)
            dsn(chn)=ek*cn(chn)
            dcn(chn)=-ek*sn(chn)
         else
            typ(chn)='c'
            ek=sqrt(-2.d0*echn(chn))
            sn(chn)=exp(-ek*rbox)
            dsn(chn)=sn(chn)-ek*rbox*sn(chn)
            sn(chn)=rbox*sn(chn)
         endif
 10   continue
      if(prn) then
         title='regular external solution'   
         call prntrm(title,sn,nc,1,nc,1,iout)
         title='derivative of regular external solution'   
         call prntrm(title,dsn,nc,1,nc,1,iout)
         title='irregular external solution'   
         call prntrm(title,cn,nc,1,nc,1,iout)
         title='derivative of irregular external solution'   
         call prntrm(title,dcn,nc,1,nc,1,iout)
      endif
      return
      end
