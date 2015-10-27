*deck conrmat.f
c***begin prologue     conrmat
c***date written       951229   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c***keywords           
c***author             schneider, barry (nsf)
c***source             
c***purpose            one dimensional r-matrix
c***                   
c***references         
c
c***routines called    
c***end prologue       conrmat
      subroutine conrmat(rmat,drmat,eigval,gama,gamb,energy,
     1                   na,nb,nc,n,prn,der)
      implicit integer (a-z)
      real*8 rmat, drmat, eigval, gama, gamb, energy
      logical prn, der
      character*80 title
      character*16 fptoc
      dimension rmat(nc,*), drmat(nc,*), eigval(n)
      dimension gama(nc,n), gamb(nc,n), prn(2)
      common/io/inp, iout
      do 10 chna=1,na
         do 20 chnb=1,nb 
            do 30 q=1,n
               rmat(chna,chnb) = rmat(chna,chnb) 
     1                            + 
     2                       gama(chna,q)*gamb(chnb,q)
     3                           /(eigval(q)-energy)
 30         continue
            rmat(chna,chnb)=.5d0*rmat(chna,chnb)
 20      continue   
 10   continue
      if(der) then
         do 40 chna=1,na
            do 50 chnb=1,nb 
               do 60 q=1,n
                  drmat(chna,chnb) = drmat(chna,chnb) 
     1                               - 
     2                          gama(chna,q)*gamb(chnb,q)
     3                             /( (eigval(q)-energy)
     4                                       *           
     5                                (eigval(q)-energy) )
 60            continue
               drmat(chna,chnb)=.5d0*drmat(chna,chnb)
 50         continue   
 40      continue
      endif
      if(prn(2)) then
         title='R(c|c'//')'
	 write(iout,1) energy
         call prntrm(title,rmat,na,nb,nc,nc,iout)
      endif
      return
 1    format(/,2x,'energy = ',d15.8)
      end       



