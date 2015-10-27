*deck @(#)svlam.f	
c***begin prologue     svlam
c***date written       920402   (yymmdd)
c***revision date      
c***keywords           svlam, link 6050
c***authors            Schneider,B (NSF)
c***                   
c***source             m6050
c***purpose            vlamda special function for optical potential
c***                   construction for two electron systems for s
c***                   waves. 
c***references       
c
c***routines called    
c***end prologue       svlam
      subroutine svlam (vlamda,grid,clm,scr,gam,del,z,fact,lval,mval,
     1                 ntrms,nlam,npt,type,spin,pvlam,pvlm)
      implicit real *8 (a-h,o-z)
      character *8 spin
      character *80 title
      complex *16 vlamda, gam, del, clm, scr, facl, facm
      complex *16 cz
      integer type
      logical pvlam, pvlm, ppass 
      dimension vlamda(npt,nlam,5), clm(ntrms,nlam), lval(ntrms)
      dimension mval(ntrms), grid(4,npt), scr(npt,2), fact(0:100)
      dimension gam(ntrms), del(ntrms)
      common /io/ inp, iout
c
      a = 2.d0 * sqrt ( 8.d0*z*z*z )
      b = 4.d0 *z*z*z
      cz = cmplx(z,0.d0)
c----------------------------------------------------------------------c
c                 factor needed to convert potential to Hartrees       c
c                 is done by multiplying vlamda by .5                  c
c----------------------------------------------------------------------c
      pre1 = .5d0 * a
      pre2 = .5d0 * a * b
      pre3 = .5d0 * a * b * b
c
      ifac = 1
      if (spin.eq.'triplet') then
          ifac = -1
      endif
      if ( type.eq.1 ) then
c          ppass=pvlm
           do 10 i=1,ntrms
              call vlm (scr(1,1),lval(i),mval(i),gam(i),del(i),z,
     1                  fact,grid,npt,ppass)
              call vlm (scr(1,2),mval(i),lval(i),del(i),gam(i),z,
     1                  fact,grid,npt,ppass)
              do 20 lam=1,nlam
                 do 30 j=1,npt
                    vlamda(j,lam,2) = vlamda(j,lam,2) +
     1                                pre1 * clm(i,lam) * 
     2                                ( scr(j,1) + ifac* scr(j,2) )
   30            continue
   20         continue
   10      continue
           if (pvlam) then
               title='first vlamda'
               call prntcmn(title,vlamda(1,1,2),npt,nlam,npt,
     1                      nlam,iout,'e')
           endif
      elseif (type.eq.2) then
           ppass=.false.
           do 40 i=1,ntrms
              call vlm(scr(1,2),0,mval(i),cz,del(i),z,fact,grid,
     1                 npt,ppass)
              call vlm(scr(1,1),0,lval(i),cz,gam(i),z,fact,grid,
     1                 npt,ppass)
              facl = fact(lval(i)+2) / ( gam(i) + z )**( lval(i)+3 ) 
              facm = fact(mval(i)+2) / ( del(i) + z )**( mval(i)+3 ) 
              do 50 lam=1,nlam
                 do 60 j=1,npt
                    vlamda(j,lam,3) = vlamda(j,lam,3) + 
     1                                pre2 * clm(i,lam) * 
     2                                ( facl * scr(j,2) + 
     3                                 ifac * facm * scr(j,1) )
   60            continue
   50         continue
   40      continue
           if (pvlam) then
               title='second vlamda'
               call prntcmn(title,vlamda(1,1,3),npt,nlam,npt,
     1                      nlam,iout,'e')
           endif
      elseif (type.eq.3) then
           ppass=.false.
           do 70 i=1,ntrms
              call vlm(scr(1,1),lval(i),0,gam(i),cz,z,fact,grid,
     1                 npt,ppass)
              call vlm(scr(1,2),mval(i),0,del(i),cz,z,fact,grid,
     1                 npt,ppass)
              facl = fact(lval(i)+2) / ( gam(i) + z )**( lval(i)+3 ) 
              facm = fact(mval(i)+2) / ( del(i) + z )**( mval(i)+3 ) 
              do 80 lam=1,nlam 
                 do 90 j=1,npt
                    vlamda(j,lam,4) = vlamda(j,lam,4) + 
     1                                pre2 * clm(i,lam) * 
     2                                ( facm * scr(j,1) +
     3                                  ifac * facl * scr(j,2) )
   90            continue
   80         continue
   70      continue
           if (pvlam) then
               title='third vlamda'
               call prntcmn(title,vlamda(1,1,4),npt,nlam,npt,
     1                      nlam,iout,'e')
           endif
      elseif (type.eq.4) then
           ppass=.false. 
           if (spin.eq.'singlet') then
               call vlm(scr(1,1),0,0,cz,cz,z,fact,grid,npt,ppass)
               do 100 i=1,ntrms
                  facl = fact(lval(i)+2) / ( gam(i) + z )**( lval(i)+3 ) 
                  facm = fact(mval(i)+2) / ( del(i) + z )**( mval(i)+3 ) 
                  do 200 lam=1,nlam
                     do 300 j=1,npt
                        vlamda(j,lam,5) = vlamda(j,lam,5) + 
     1                                    2.d0 * pre3 * clm(i,lam) *
     2                                    facl *facm *scr(j,1)
  300                continue
  200             continue
  100          continue
           endif
           if (pvlam) then
               title='fourth vlamda'
               call prntcmn(title,vlamda(1,1,5),npt,nlam,npt,
     1                      nlam,iout,'e')
           endif
      endif
      return
      end














