*deck @(#)nsvlam.f	
c***begin prologue     nsvlam
c***date written       920402   (yymmdd)
c***revision date      
c***keywords           nsvlam, link 6050
c***authors            Schneider,B (NSF)
c***                   
c***source             m6050
c***purpose            vlamda special function for optical potential
c***                   construction for two electron systems for non
c***                   s waves.
c***references       
c
c***routines called    
c***end prologue       nsvlam
      subroutine nsvlam (vlamda,pl,grid,clm,scr,gam,del,fact,lval,
     1                   mval,ntrms,angmom,nlam,npt,type,spin,
     2                   pvlam,pvlm)
      implicit real*8 (a-h,o-z)
      character*8 spin
      character*80 title
      complex*16 vlamda, gam, del, clm, scr, gamp1, delp1, pmul, ctwo
      integer type
      logical pvlam, pvlm, ppass 
      dimension vlamda(npt,nlam), clm(ntrms,nlam), lval(ntrms)
      dimension mval(ntrms), grid(4,npt), scr(npt), fact(0:100)
      dimension gam(ntrms), del(ntrms), pl(npt)
      common /io/ inp, iout
      data pi /3.14159265358979323846d+00/
      data ctwo / (2.d0,0.d0) /
c
      ifac = 1
      if (spin.eq.'triplet') then
          ifac = -1
      endif
      call czero(scr,npt)
c          ppass=pvlm
      do 10 i=1,ntrms
         gamp1=gam(i)+1.d0
         delp1=del(i)+1.d0
         pmul=dcmplx(1.d0,0.d0)
         call elm (scr,lval(i),mval(i),gamp1,del(i),pmul,fact,grid,
     1             npt,ppass)
         pmul=-4.d0*fact(lval(i)+2) / (gamp1)**( lval(i)+3)
         call elm (scr,0,mval(i),ctwo,del(i),pmul,fact,grid,
     1             npt,ppass)
         pmul=ifac
         call elm(scr,mval(i),lval(i),delp1,gam(i),pmul,fact,grid,
     1            npt,ppass)
         pmul=-ifac*4.d0*fact(mval(i)+2) /(delp1)**(mval(i)+3)
         call elm(scr,0,lval(i),ctwo,gam(i),pmul,fact,grid,
     1            npt,ppass)
         do 20 lam=1,nlam
            do 30 j=1,npt
               vlamda(j,lam) = vlamda(j,lam) + clm(i,lam) * scr(j)
   30       continue
   20    continue
   10 continue
c     normaliztion for phi angular part for m=0
      rmul=1.d0/sqrt(2.d0*pi)
      call smul(vlamda,vlamda,rmul,2*nlam*npt)
      do 40 lam=1,nlam
         do 50 j=1,npt
            vlamda(j,lam) = vlamda(j,lam)*pl(j)
   50    continue
   40 continue      
      if (pvlam) then
          title='vlamda'
          call prntcmn(title,vlamda,npt,nlam,npt,nlam,iout,'e')
      endif
      return
      end














