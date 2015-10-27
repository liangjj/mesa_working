*deck genwpt.f
c***begin prologue     genwpt
c***date written       951229   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c***keywords           
c***author             schneider, barry (nsf)
c***source             
c***purpose            lobatto points and weights.
c***                   
c***references         
c
c***routines called    
c***end prologue       genwpt
      subroutine genwpt(ppt,xl,xr,ngrid,npt,qdtyp)
      implicit integer (a-z)
      real*8 pt, sc, xl, xr
      logical nfix
      character*(*) qdtyp
      character*80 title
      dimension nfix(2)
      data nfix/.true.,.true./
      common/io/inp, iout
      pointer (ppt,pt(1))
      pointer (psc,sc(1))
      mxpt=0
      do 10 grd=1,ngrid
         mxpt=max(mxpt,npt(grd))
 10   continue   
      call memory(wpadti(1+mxpt*mxpt),psc,nscr,'scratch',0)
      begin=1
      do 20 grd=1,ngrid
         q=begin
         wt=q+npt(grd)
         begin=wt+npt(grd)
         call getqpt(pt(q),pt(wt),xl,xr,qdtyp,'before',
     1               sc,nfix,npt(grd),npt(grd),1,.false.)
 20   continue   
      call memory(-nscr,psc,idum,'scratch',idum)
      return
      end       


