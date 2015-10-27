*deck @(#)swlam.f	
c***begin prologue     swlam
c***date written       920402   (yymmdd)
c***revision date      
c***keywords           swlam, link 6050
c***authors            Schneider,B (NSF)
c***                   
c***source             m6050
c***purpose            wlamda special function for optical potential
c***                   construction for two electron systems for s
c***                   waves.
c***references       
c
c***routines called    
c***end prologue       swlam
      subroutine swlam (wlamda,grid,clm,scr,gam,del,z,fact,lval,mval,
     1                  ntrms,nlam,npt,spin,pwlam,pvlm)
      implicit real *8 (a-h,o-z)
      character *8 spin
      character *80 title
      complex *16 wlamda, total
      complex *16  ilm, gam, del, clm, gampz, delpz, term, ctwoz
      logical pwlam, pvlm, ppass 
      dimension wlamda(npt,nlam), clm(ntrms,nlam), lval(ntrms)
      dimension scr(npt), mval(ntrms), grid(4,npt), fact(0:100), term(7)
      dimension gam(ntrms), del(ntrms)
      common /io/ inp, iout
c
      do 10 i=1,npt
         rval=sqrt( grid(1,i)*grid(1,i) + grid(2,i)*grid(2,i) +
     1              grid(3,i)*grid(3,i) )         
         scr(i) =  rval * exp (-z*rval)
   10 continue
c     
      a = 8.d0 * z * z * z * sqrt( 2.d0*z*z*z )
      b = 4.d0 * z * z * z
      ctwoz = dcmplx(2.d0*z,0.d0)
c----------------------------------------------------------------------c
c                 factor needed to convert potential to Hartrees       c
c                 is done by multiplying wlamda by .5                  c
c----------------------------------------------------------------------c
      pre = .5d0 * a
c
      ifac = 1
      if (spin.eq.'triplet') then
          ifac = -1
      endif
c          ppass=pvlm
      call czero(wlamda,npt)
      do 20 i=1,ntrms
         l2=lval(i)+2
         l3=lval(i)+3
         m2=mval(i)+2
         m3=mval(i)+3
         gampz = gam(i) + z
         delpz = del(i) + z
         term(1) = ilm(lval(i),mval(i),gampz,delpz,fact,ppass)
         term(2) = ifac*ilm(mval(i),lval(i),delpz,gampz,fact,ppass)
         term(3) = -b*ilm(0,mval(i),ctwoz,delpz,fact,ppass)
     1                     *fact(l2)/gampz**l3
         term(4) = -b*ifac*ilm(0,lval(i),ctwoz,gampz,fact,ppass)
     1                     *fact(m2)/delpz**m3   
         term(5) = -b*ilm(lval(i),0,gampz,ctwoz,fact,ppass)
     1                     *fact(m2)/delpz**m3
         term(6) = -b*ifac*ilm(mval(i),0,delpz,ctwoz,fact,ppass)
     1                     *fact(l2)/gampz**l3    
         total=term(1)+term(2)+term(3)+term(4)+term(5)+term(6)
         term(7) = 0.d0
         if (spin.eq.'singlet') then
             term(7) = 4.d0* b * b * ilm(0,0,ctwoz,ctwoz,fact,ppass)
             term(7) = fact(l2) * fact(m2) *term(7) /
     1                      ( ( gampz**l3 ) * ( delpz**m3 ) )      
         endif
         total = pre * ( total + term(7) )
         do 30 lam=1,nlam
            do 40 j=1,npt
               wlamda(j,lam) = wlamda(j,lam) +
     1                          clm(i,lam) * total * scr(j)
   40       continue
   30    continue
   20 continue     
      if (pwlam) then
          title='wlamda'
          call prntcmn(title,wlamda,npt,nlam,npt,
     1                 nlam,iout,'e')
      endif
      return
      end














