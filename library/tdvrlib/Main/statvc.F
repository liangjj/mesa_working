*deck statvc.f
c***begin prologue     statvc
c***date written       960723   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c***keywords           plot, time-dependent, wavefunction
c***author             schneider, barry (nsf)
c***source             
c***purpose            plot or massage the complex time-dependent 
c***                   for analysis wavefunction.
c***                   
c***references         
c
c***routines called    
c***end prologue       statvc
      subroutine statvc(psi,psi0,q1,q2,q3,t,t0,pt,energy,dim,
     1                  nd,n,m,nt,nc,pnch)
      implicit integer (a-z)
      real*8 psi, psi0, q1, q2, q3, t, t0, pt
      real*8 energy
      real*8 repsi, impsi, psimax, psimin, fac
      complex*16 exact, approx, eye, cfac
      character*3 itoc
      logical pnch
      dimension nd(dim)
      dimension psi(n,nc,2), psi0(m,nc,2)
      dimension q1(nd(1)), q2(nd(2)), q3(nd(3)), t(nt), fac(2)
      dimension pt(nt,nt), cfac(3)
      common/io/inp, iout
      common/punch/pun(3)
      data eye/(0.d0,1.d0)/
      cfac(2)=exp(-eye*energy*t0)
      cfac(3)=1.d0/cfac(2)
      write(iout,1) t0
      if(dim.eq.1) then
         ij=0
         do 10 i=1,nt
            cfac(1) =  exp(-eye*energy*t(i)) 
            write(iout,2) t(i)
            do 20 j=1,nd(1)
               ij=ij+1
               fac(1) = psi0(j,1,1)
               fac(2) = psi0(j,1,2)
               exact = ( fac(1) + eye*fac(2) )*cfac(1)*cfac(3)
               psi(ij,1,1) = psi(ij,1,1)*pt(i,i) + fac(1)
               psi(ij,1,2) = psi(ij,1,2)*pt(i,i) + fac(2)
               approx = psi(ij,1,1) + eye*psi(ij,1,2) 
c               psi(ij,1,1)=real(exact)
c               psi(ij,1,2)=imag(exact)
               write(iout,3) q1(j)
               write(iout,4) approx, exact
 20         continue
 10      continue 
      elseif(dim.eq.2) then
         ijk=0
         do 30 i=1,nt
            cfac(1) =  exp(-eye*energy*t(i)) 
            write(iout,2) t(i)
            jk=0
            do 40 j=1,nd(1)
               do 50 k=1,nd(2)
                  ijk=ijk+1
                  jk=jk+1
                  fac(1) = psi0(jk,1,1)
                  fac(2) = psi0(jk,1,2)
                  exact = ( fac(1) + eye*fac(2) )*cfac(1)*cfac(3)
                  psi(ijk,1,1) = psi(ijk,1,1)*pt(i,i) + fac(1)
                  psi(ijk,1,2) = psi(ijk,1,2)*pt(i,i) + fac(2)
                  approx = psi(ijk,1,1) + eye*psi(ijk,1,2) 
                  psi(ijk,1,1)=real(exact)
                  psi(ijk,1,2)=imag(exact)
                  write(iout,5) q1(j), q2(k)
                  write(iout,4) approx, exact
 50            continue   
 40         continue
 30      continue  
      elseif(dim.eq.3) then
         ijkl=0
         do 60 i=1,nt
            cfac(1) =  exp(-eye*energy*t(i)) 
            write(iout,2) t(i)
            jkl=0
            do 70 j=1,nd(1)
               do 80 k=1,nd(2)
                  do 90 l=1,nd(3)
                     ijkl=ijkl+1
                     jkl=jkl+1
                     fac(1) = psi0(jkl,1,1)
                     fac(2) = psi0(jkl,1,2)
                     exact = ( fac(1) + eye*fac(2) )*cfac(1)*cfac(3)
                     psi(ijkl,1,1) = psi(ijkl,1,1)*pt(i,i) + fac(1)
                     psi(ijkl,1,2) = psi(ijkl,1,2)*pt(i,i) + fac(2)
                     approx = psi(ijkl,1,1) + eye*psi(ijkl,1,2) 
                     psi(ijkl,1,1)=real(exact)
                     psi(ijkl,1,2)=imag(exact)
                     write(iout,6) q1(j), q2(k), q3(l)
                     write(iout,4) approx, exact
 90               continue   
 80            continue   
 70         continue
 60      continue    
      else
         call lnkerr('error in dimensions')
      endif
c      if(pnch) then
c         call plotit(psi(1,1,1),q1,q2,q3,'real',t,dim,nd,n,nt,nc,pun(1))
c         call plotit(psi(1,1,2),q1,q2,q3,'imaginary',t,dim,nd,n,nt,
c     1               nc,pun(2))
c         do 100 i=1,n
c            psi(i,1,1) = psi(i,1,1)*psi(i,1,1) + psi(i,1,2)*psi(i,1,2)
c 100     continue   
c         call plotit(psi(1,1,1),q1,q2,q3,'absolute',t,dim,nd,n,nt,
c     1               nc,pun(3))
c      endif
      return
 1    format(/,1x,'comparison between approximate and exact solution ',
     1       /,1x,'for right hand side an eigenstate',/,1x,
     2            'initial time = ',e15.8)
 2    format(/,1x,'time        = ',e15.8)
 3    format(/,1x,'x = ',e15.8)
 4    format(/,10x,'approximate = ',e15.8,1x,e15.8,/,10x, 
     1             'exact       = ',e15.8,1x,e15.8)                   
 5    format(/,1x,'x = ',e15.8,1x,'y = ',e15.8)
 6    format(/,1x,'x = ',e15.8,1x,'y = ',e15.8,1x,'z = ',e15.8)
      end       





