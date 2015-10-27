*deck well 
c***begin prologue     well 
c***date written       921701   (yymmdd)
c***revision date               (yymmdd)
c***keywords           m6005, link 6005, model potential
c***author             schneider, barry (nsf)
c***source             m6005
c***purpose            complex square well model potential 
c***references         none
c
c***routines called    none
c***end prologue       well 
      subroutine well(v,vcple,grid,size,nchan,npnts,ntri)
      implicit integer (a-z)
      parameter (ncmx=9)
      real *8 grid, size, tstsiz
      real *8 v, vcple
      dimension v(npnts,ntri), grid(4,npnts), vcple(ncmx,ncmx)
      do 10 i=1,npnts
         tstsiz=grid(1,i)*grid(1,i) + grid(2,i)*grid(2,i) +
     1                     grid(3,i)*grid(3,i)
         tstsiz=sqrt(tstsiz)
         if (tstsiz.le.size) then
	     jk=0
             do 20 j=1,nchan
	        do 30 k=1,j
                   jk=jk+1	
	           v(i,jk)=vcple(j,k)
   30           continue
   20        continue
         else
             jk=0
             do 40 j=1,nchan
	        do 50 k=1,j
                   jk=jk+1	
                   v(i,jk)=0.d+00
   50           continue
   40        continue
         endif
   10 continue
      return
      end






















