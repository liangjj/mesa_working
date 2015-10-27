*deck trevod.f
c***begin prologue     trevod
c***date written       971214   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c***keywords           polynomials
c***author             schneider, barry (nsf)
c***source             
c***purpose            special routine to transform matrix to good parity
c***                   functions.
c***                                                          
c***references         
c
c***routines called    
c***end prologue       trevod
      subroutine trevod(t,neven,nodd,n)
      implicit integer (a-z)
      real*8 t, fac
      character*80 title
      dimension t(n,n)
      common/io/inp, iout
      call rzero(t,n*n)
      write(iout,1)
      fac=1.d0/sqrt(2.d0)
c
c     determine if n is even or odd.
c
      ntst=mod(n,2)
      if(ntst.eq.0) then
         neven=n/2
         nodd=neven
         nback=n
         count=neven
         do 10 i=1,neven
            count=count+1
            t(i,i)=fac
            t(nback,i)=fac
            t(i,count)=fac
            t(nback,count)=-fac
            nback = nback - 1   
 10      continue
      else
         nodd=n/2
         neven=nodd+1
         nback = n
         count=neven
         do 20 i=1,nodd
            count=count+1
            t(i,i)=fac
            t(nback,i)=fac
            t(i,count)=fac
            t(nback,count)=-fac
            nback = nback - 1   
 20      continue
         t(nback,nback)=1.d0
      endif   
      write(iout,2) neven, nodd
      return
 1    format(/,5x,'transforming hamiltonian to functions of'
     1            ' good parity')
 2    format(/,5x,'number of even parity functions = ',i4,/,5x,
     1            'number of odd parity functions  = ',i4)
      end       


