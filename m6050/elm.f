*deck @(#)elm.f	
c***begin prologue     elm
c***date written       920402   (yymmdd)
c***revision date      
c***keywords           elm, link 6050
c***authors            Schneider,B (NSF)
c***                   
c***source             m6050
c***purpose            special functions for optical potential
c***                   construction for two electron systems
c***references       
c
c***routines called    
c***end prologue       elm
      subroutine elm (ef,l,m,gam,del,pmul,fact,grid,npt,prnt)
      implicit real *8 (a-h,o-z)
      common /io/ inp, iout
      logical prnt
      complex*16 ef, gam, del, pre, expfac, term1, sum, gamfac
      complex*16 start, pmul
      character *80 title
      character *2 itoc
      dimension ef(npt), grid(4,npt), fact(0:100)
c
      l1=l+1
      l2=l+2
      l3=l+3
      start=dcmplx(0.d0,0.d0)
      pre= fact(l2) / ( gam )**(l3)
      do 10 i=1,npt
         rval=sqrt( grid(1,i)*grid(1,i) + grid(2,i)*grid(2,i) +
     1              grid(3,i)*grid(3,i) )
         expfac=exp(-gam*rval)
         term1= ( 1.d0 -expfac) * pre
         sum=start
         rfac=rval**l1
         gamfac=gam*gam
         do 20 j=1,l1
            sum = sum + j * ( fact(l1) / fact(l2-j) )
     1                    * ( rfac / gamfac )
            rfac=rfac/rval
            gamfac=gamfac*gam
   20    continue
         sum = sum * expfac
         ef(i) = ef(i) + pmul * 2.d0 * ( term1 - sum ) * exp(-rval*del)
     1                       * rval**m
   10 continue
      if (prnt) then
          title='elm-'//itoc(l)
          call prntcmn(title,ef,npt,1,npt,1,iout,'e')
      endif
      return
      end






