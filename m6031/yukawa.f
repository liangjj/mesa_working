*deck yukawa 
c***begin prologue     yukawa 
c***date written       921801   (yymmdd)
c***revision date               (yymmdd)
c***keywords           m6005, link 6005, model potential
c***author             schneider, barry (nsf)
c***source             m6005
c***purpose            complex yukawa model potential 
c***references         none
c
c***routines called    none
c***end prologue       yukawa 
      subroutine yukawa(v,vcple,grid,nchan,npnts,ntri)
      implicit integer (a-z)
      parameter (ncmx=9)
      real *8 grid, rsq, r, v, vcple
      dimension v(npnts,ntri), grid(4,npnts), vcple(ncmx,ncmx)
      jk=0
      do 30 j=1,nchan
         do 20 k=1,j
            jk=jk+1	
            do 10 i=1,npnts
	       rsq=grid(1,i)*grid(1,i)+grid(2,i)*grid(2,i)+
     1     	                       grid(3,i)*grid(3,i)
	       r=sqrt(rsq)
	       v(i,jk)=-exp(vcple(j,k)*r)/r
   10       continue
   20    continue
   30 continue
      return
      end
