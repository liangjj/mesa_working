*deck anlyze.f
c***begin prologue     anlyze
c***date written       960723   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c***keywords           plot, time-dependent, wavefunction
c***author             schneider, barry (nsf)
c***source             
c***purpose            analyze and/or plot the complex time-dependent 
c***                   wavefunction.
c***                   
c***references         
c
c***routines called    
c***end prologue       anlyze
      subroutine anlyze(psi,psi0,q1,q2,q3,t,p1,p2,p3,pt,u1,u2,u3,
     1                  t0,dim,nd,n,m,nt,nc,proj,pnch)
      implicit integer (a-z)
      real*8 psi, psi0, q1, q2, q3, t, p1, p2, p3, pt, t0
      real*8 gvec, psisq
      logical proj, pnch
      dimension nd(dim)
      dimension psi(n,nc,2), psi0(m,nc,2)
      dimension q1(nd(1)), q2(nd(2)), q3(nd(3))
      dimension t(nt), p1(nd(1),nd(1)), p2(nd(2),nd(2))
      dimension p3(nd(3),nd(3)), pt(nt,nt)
      dimension u1(nd(1),nd(1)), u2(nd(2),nd(2)), u3(nd(3),nd(3))
      common/io/inp, iout
      common/punch/pun
      pointer (pgvec,gvec(1))
      write(iout,1) t0
      if(proj) then
         write(iout,2)
         if(dim.eq.1) then
            need=wptoin(nd(1)) 
            call getmem(need,pgvec,ngot,'gvec',0)
            call vcoef1(psi,u1,gvec,t,nd(1),nt,nc)
            call getmem(-ngot,pgvec,idum,'gvec',idum)
         elseif(dim.eq.2) then
            need=wptoin(nd(2)*nd(1)*2) 
            call getmem(need,pgvec,ngot,'gvec',0)
            call vcoef2(psi,u1,u2,gvec,t,nd(1),nd(2),nt,nc)
            call getmem(-ngot,pgvec,idum,'gvec',idum)
         elseif(dim.eq.3) then
            need=wptoin(nd(3)*nd(2)*nd(1)) 
            call getmem(need,pgvec,ngot,'gvec',0)
            call vcoef2(psi,u1,u2,u3,gvec,t,nd(1),nd(2),nd(3),nt,nc)
            call getmem(-ngot,pgvec,idum,'gvec',idum)
         else
            call lnkerr('error in dimension')
         endif 
      endif
      if(dim.eq.1) then
         do 10 ic=1,nc
            write(iout,3) ic
            ij=0
            do 20 i=1,nt
               write(iout,4) t(i)
               do 30 j=1,nd(1)
                  ij=ij+1
                  psi(ij,ic,1) = psi(ij,ic,1)*pt(i,i) + psi0(j,ic,1)
                  psi(ij,ic,2) = psi(ij,ic,2)*pt(i,i) + psi0(j,ic,2)
                  psisq = psi(ij,ic,1)*psi(ij,ic,1) 
     1                                + 
     2                    psi(ij,ic,2)*psi(ij,ic,2)
                  write(iout,5) q1(j), psi(ij,ic,1), psi(ij,ic,2), psisq
 30            continue   
 20         continue
 10      continue 
      elseif(dim.eq.2) then
         do 40 ic=1,nc
            write(iout,3) ic
            ijk=0
            do 50 i=1,nt
               write(iout,4) t(i)
               jk=0
               do 60 j=1,nd(1)
                  do 70 k=1,nd(2)
                     jk=jk+1
                     ijk=ijk+1
                     psi(ijk,ic,1) = psi(ijk,ic,1)*pt(i,i) 
     1                                            + 
     2                               psi0(jk,ic,1)
                     psi(ijk,ic,2) = psi(ijk,ic,2)*pt(i,i) 
     1                                            + 
     2                               psi0(jk,ic,2)
                     psisq = psi(ijk,ic,1)*psi(ijk,ic,1) 
     1                                   + 
     2                       psi(ijk,ic,2)*psi(ijk,ic,2)
                     write(iout,6) q1(j), q2(k), psi(ijk,ic,1), 
     1                                           psi(ijk,ic,2), psisq
 70               continue
 60            continue   
 50         continue
 40      continue 
      elseif(dim.eq.3) then
         do 80 ic=1,nc
            write(iout,3) ic
            ijkl=0
            do 90 i=1,nt
               write(iout,4) t(i)
               jkl=0
               do 100 j=1,nd(1)
                  do 110 k=1,nd(2)
                     do 120 l=1,nd(3)
                        jkl=jkl+1
                        ijkl=ijkl+1
                        psi(ijkl,ic,1) = psi(ijkl,ic,1)*pt(i,i) 
     1                                               + 
     2                                   psi0(jkl,ic,1)
                        psi(ijkl,ic,2) = psi(ijkl,ic,2)*pt(i,i) 
     1                                               + 
     2                                   psi0(jkl,ic,2)
                        psisq = psi(ijkl,ic,1)*psi(ijkl,ic,1) 
     1                                      + 
     2                          psi(ijkl,ic,2)*psi(ijkl,ic,2)
                        write(iout,7) q1(j), q2(k), q3(l), 
     1                                psi(ijkl,ic,1), psi(ij,ic,2), 
     2                                psisq
 120                 continue
 110              continue   
 100           continue
 90         continue   
 80      continue
      else
         call lnkerr('error in dimensions')
      endif
      if(pnch) then
         do 200 ic=1,nc
            write(iout,8) ic
            call plotit(psi(1,ic,1),q1,q2,q3,'real',t,dim,nd,n,nt,nc)
            call plotit(psi(1,ic,2),q1,q2,q3,'imaginary',t,dim,nd,n,
     1                  nt,nc)
            do 300 i=1,n
               psi(i,ic,1) = psi(i,ic,1)*psi(i,ic,1) 
     1                                   + 
     2                       psi(i,ic,2)*psi(i,ic,2)
 300        continue   
            call plotit(psi(1,ic,1),q1,q2,q3,'abs',t,dim,nd,n,nt,nc)
 200     continue   
      endif
      return
 1    format(/,1x,'analysis of time dependent wavefunction',/,1x,
     1            'initial time = ',e15.8)
 2    format(/,1x,'eigenstate analysis as a function of time')
 3    format(/,1x,'tabulate and plot wavefunction for channel = ',i3)
 4    format(/,10x,'time = ',e15.8)
 5    format(/,1x,' x = ',e15.8,
     1      /,1x,'real psi = ',e15.8,1x,'imag psi = ',e15.8,
     2      /,1x,'|psi|**2 = ',e15.8)
 6    format(/,1x,' x = ',e15.8,1x,'y = ',e15.8,
     1      /,1x,'real psi = ',e15.8,1x,'imag psi = ',e15.8,
     2      /,1x,'|psi|**2 = ',e15.8)
 7    format(/,1x,' x = ',e15.8,1x,'y = ',e15.8,1x,'z = ',e15.8,
     1      /,1x,'real psi = ',e15.8,1x,'imag psi = ',e15.8,
     2      /,1x,'|psi|**2 = ',e15.8)
 8    format(/,1x,'plotting psi for channel = ',i3)
      end       








