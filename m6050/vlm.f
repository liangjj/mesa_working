*deck @(#)vlm.f	
c***begin prologue     vlm
c***date written       920402   (yymmdd)
c***revision date      
c***keywords           vlm, link 6050
c***authors            Schneider,B (NSF)
c***                   
c***source             m6050
c***purpose            special functions for optical potential
c***                   construction for two electron systems
c***references       
c
c***routines called    
c***end prologue       vlm
      subroutine vlm (vf,l,m,gam,del,z,fact,grid,npt,prnt)
      implicit real *8 (a-h,o-z)
      common /io/ inp, iout
      logical prnt
      complex *16 vf, gam, del, exfac, pre1, pre2, sum
      character *80 title
      character *2 itoc
      dimension vf(npt), grid(4,npt), fact(0:100)
c
      exfac=gam+z
      pre1=fact(l+2) / ( gam + z)**(l+3)
      if (l.eq.0) then
          do 10 i=1,npt
             rval=sqrt( grid(1,i)*grid(1,i) + grid(2,i)*grid(2,i) +
     1                  grid(3,i)*grid(3,i) )
             vf(i) = pre1 + rval / ( (gam+z)*(gam+z) ) 
             vf(i) = vf(i) * exp(-exfac*rval)  
             vf(i) = ( pre1 - vf(i) ) * ( rval**m ) * exp( -del*rval)
   10     continue
      else
          do 20 i=1,npt
             rval=sqrt( grid(1,i)*grid(1,i) + grid(2,i)*grid(2,i) +
     1                  grid(3,i)*grid(3,i) )
             pre2 = rval**(l+1) / ( (gam+z)*(gam+z) ) 
             sum = pre2 / fact(l+1) 
             do 30 j=2, l+1
                pre2 = pre2 / ( rval * ( gam + z ) ) 
                sum= sum + j * pre2 / fact(l+2-j)
   30        continue
             sum = fact(l+1) * sum
             vf(i) = pre1 - ( pre1 + sum ) * exp (-exfac*rval)
             vf(i) = vf(i) * exp(-del*rval) * ( rval**m )
   20     continue
      endif
      if (prnt) then
          title='vlm-'//itoc(l)//'-'//itoc(m)
          call prntcmn(title,vf,npt,1,npt,1,iout,'e')
      endif
      return
      end






