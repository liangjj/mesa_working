*deck zproj.f
c***begin prologue     zproj
c***date written       960723   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c***keywords           time development
c***author             schneider, barry (nsf)
c***source             
c***purpose            tabulate gaussian wavepacket 
c***                   on grid
c***references         
c
c***routines called    
c***end prologue       zproj
      subroutine zproj(q,wt,p,sigma,x_0,beta,z,n,typke,coord,prn)
      implicit integer (a-z)
      real*8 q, wt, p, sigma, x_0, beta
      complex*16 z, eye
      character*(*) typke, coord
      character*80 title
      logical prn
      dimension p(n,n), q(n), wt(n), z(n)
      data eye / (0.d0,1.d0) /
      common/io/inp, iout
      if(coord.eq.'cartesian') then
         do 10 i=1,n
            z(i) = exp( - ( q(i) - x_0 ) * ( q(i) - x_0 ) 
     1                                   / ( 2.d0*sigma*sigma) 
     2                  + eye*beta*q(i) )
 10      continue
      else
         do 20 i=1,n
            z(i) = q(i)*exp( - ( q(i) - x_0 ) * ( q(i) - x_0 )
     1                                   / ( 2.d0*sigma*sigma) 
     2                       + eye*beta*q(i) )
 20      continue
      endif
      if(typke.eq.'dvr'.or.typke.eq.'packed') then
         do 30 i=1,n
            z(i) = z(i) * wt(i) * p(i,i)
 30      continue
      endif   
      if(prn) then
         title='initial wavepacket'
         call prntcm(title,z,n,1,n,1,iout)
      endif
      return
      end       

