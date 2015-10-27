*deck plotit.f
c***begin prologue     plotit
c***date written       960723   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c***keywords           plot, time-dependent, wavefunction
c***author             schneider, barry (nsf)
c***source             
c***purpose            plot wavefunction.
c***                   
c***references         
c
c***routines called    
c***end prologue       plotit
      subroutine plotit(psi,q1,q2,q3,type,t,dim,nd,n,nt,nc,pun)
      implicit integer (a-z)
      real*8 psi, q1, q2, q3, t
      real*8 psimax, psimin
      character*(*) type
      dimension nd(dim)
      dimension psi(n,nc)
      dimension q1(nd(1)), q2(nd(2)), q3(nd(3))
      dimension t(nt)
      common/io/inp, iout	 
      locte=0
      psimax=0.d0
      psimin=0.d0
      if(dim.eq.0) then
         do 10 i=1,nt,5
c            write(pun,1) t(i), psi(i,1)
            write(pun,1) psi(i,1)
 10      continue   
      elseif(dim.eq.1) then
         write(iout,*) pun
         do 20 i=1,nt
            write(pun,2) t(i)
            do 30 j=1,nd(1)
               locte=locte+1
               psimax=max(psimax,psi(locte,1))
               psimin=min(psimin,psi(locte,1))
               write(pun,1) q1(j), psi(locte,1)
 30         continue   
 20      continue   
      elseif(dim.eq.2) then
         do 40 i=1,nt
            write(pun,2) t(i)
            do 50 j=1,nd(1)
               do 60 k=1,nd(2)
                  locte=locte+1
                  psimax=max(psimax,psi(locte,1))
                  psimin=min(psimin,psi(locte,1))
                  write(pun,3) q1(j), q2(k), psi(locte,1)
 60            continue   
 50         continue
 40      continue   
      elseif(dim.eq.3) then
         do 70 i=1,nt
            write(pun,2) t(i)
            do 80 j=1,nd(1)
               do 90 k=1,nd(2)
                  do 100 l=1,nd(3)
                     locte=locte+1
                     psimax=max(psimax,psi(locte,1))
                     psimin=min(psimin,psi(locte,1))
                     write(pun,4) q1(j), q2(k), q3(l), psi(locte,1)
 100              continue   
 90            continue   
 80         continue
 70      continue   
      else
         call lnkerr('error in dimensions')
      endif
      write(iout,5) psimax, psimin
      return
 1    format(e15.8,1x,e15.8)
 2    format(e15.8)
 3    format(e15.8,1x,e15.8,1x,e15.8)
 4    format(e15.8,1x,e15.8,1x,e15.8,1x,e15.8 )
 5    format(/,1x,'psimax = ',e15.8,1x,'psimin = ',e15.8)
      end       
