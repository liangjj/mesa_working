*deck lschr.f
c***begin prologue     lschr
c***date written       951229   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c***keywords           3-d schroedinger equation
c***author             schneider, barry (nsf)
c***source             trap3d
c***purpose            calculate hamiltonian matrix explicitly 
c***                   and diagonalize on request.
c***                   
c***                   
c***references         
c
c***routines called    
c***end prologue       lschr
      subroutine lschr(h01,ch01,h02,ch02,h03,ch03,v,ham,cham,eig,
     1                 ceig,vec,vecl,vecr,cvecl,cvecr,
     2                 work,lwork,rwork,ind,n,n1,n2,n3,dim,withv,
     3                 op,prnt,mattyp)
      implicit integer (a-z)
      real*8 h01, h02, h03, v, ham, eig, vec
      real*8 vecl, vecr, rwork
      complex*16 ch01, ch02, ch03, cham, ceig, cvecl, cvecr, work
      character*80 title
      character*(*) mattyp
      logical withv, op, prnt
      dimension h01(n1,n1), h02(n2,n2), h03(n3,n3)
      dimension ch01(n1,n1), ch02(n2,n2), ch03(n3,n3)
      dimension v(n), ham(n,n), eig(n), ind(n,*)
      dimension cham(n,n), ceig(n), vecl(n,n), vecr(n,n) 
      dimension cvecl(n,n), cvecr(n,n)
      dimension work(*), rwork(*)
      common/io/inp, iout
c
      if(mattyp.eq.'complex'.or.mattyp.eq.'real-unsymmetric') then
         call czero(cham,n*n)
         if(withv) then            
            do 10 i=1,n
               cham(i,i) = v(i) 
 10         continue
         endif            
         if(dim.eq.1) then
            title='1d'
            do 20 i=1,n
               do 30 j=1,n
                  cham(i,j) = cham(i,j) + ch01(i,j)
 30            continue
 20         continue
         elseif(dim.eq.2) then
            title='2d'
            cnti=0
            do 40 i=1,n
               cnti=cnti+1
               i1=ind(cnti,1)
               j1=ind(cnti,2)
               cntj=0
               do 50 j=1,n
                  cntj=cntj+1
                  i2=ind(cntj,1)
                  j2=ind(cntj,2)
                  di1i2=0
                  dj1j2=0
                  if(i1.eq.i2) then       
                     di1i2=1
                  endif 
                  if(j1.eq.j2) then
                     dj1j2=1
                  endif
                  cham(cnti,cntj) = cham(cnti,cntj) + ch01(i1,i2)*dj1j2
     1                                              + ch02(j1,j2)*di1i2
 50            continue
 40         continue   
         elseif(dim.eq.3) then   
            title='3d'
            cnti=0
            do 60 i=1,n
               cnti=cnti+1
               i1=ind(cnti,1)
               j1=ind(cnti,2)
               k1=ind(cnti,3)
               cntj=0
               do 70 j=1,n
                  cntj=cntj+1
                  i2=ind(cntj,1)
                  j2=ind(cntj,2)
                  k2=ind(cntj,3)
                  di1i2=0
                  dj1j2=0
                  dk1k2=0
                  if(i1.eq.i2) then
                     di1i2=1
                  endif
                  if(j1.eq.j2) then
                     dj1j2=1
                  endif
                  if(k1.eq.k2) then
                     dk1k2=1
                  endif
                  cham(cnti,cntj) = cham(cnti,cntj) + 
     1                              ch01(i1,i2)*dj1j2*dk1k2 +
     2                              ch02(j1,j2)*di1i2*dk1k2 +
     3                              ch03(k1,k2)*di1i2*dj1j2
 70            continue
 60         continue   
         endif
         if(prnt) then
            title='hamiltonian-'//title(1:2)
            call prntcm(title,cham,n,n,n,n,iout)
         endif                        
         if(op) then
            write(iout,1) title
            call rdiag(ham,cham,eig,ceig,vec,cvecl,cvecr,
     1                 work,lwork,rwork,n,n,mattyp)
         endif
      else               
         call rzero(ham,n*n)
         if(withv) then            
            do 100 i=1,n
               ham(i,i) = v(i) 
 100        continue
         endif            
         if(dim.eq.1) then
            title='1d'
            do 200 i=1,n
               do 300 j=1,i
                  ham(i,j) = ham(i,j) + h01(i,j)
 300           continue
 200        continue   
         elseif(dim.eq.2) then
            title='2d'
            cnti=0
            do 400 i=1,n
               cnti=cnti+1
               i1=ind(cnti,1)
               j1=ind(cnti,2)
               cntj=0
               do 500 j=1,i
                  cntj=cntj+1
                  i2=ind(cntj,1)
                  j2=ind(cntj,2)
                  di1i2=0
                  dj1j2=0
                  if(i1.eq.i2) then       
                     di1i2=1
                  endif 
                  if(j1.eq.j2) then
                     dj1j2=1
                  endif
                  ham(cnti,cntj) = ham(cnti,cntj) + h01(i1,i2)*dj1j2
     1                                            + h02(j1,j2)*di1i2
 500            continue
 400     continue   
         elseif(dim.eq.3) then   
            title='3d'
            cnti=0
            do 600 i=1,n
               cnti=cnti+1
               i1=ind(cnti,1)
               j1=ind(cnti,2)
               k1=ind(cnti,3)
               cntj=0
               do 700 j=1,i
                  cntj=cntj+1
                  i2=ind(cntj,1)
                  j2=ind(cntj,2)
                  k2=ind(cntj,3)
                  di1i2=0
                  dj1j2=0
                  dk1k2=0
                  if(i1.eq.i2) then
                     di1i2=1
                  endif
                  if(j1.eq.j2) then
                     dj1j2=1
                  endif
                  if(k1.eq.k2) then
                     dk1k2=1
                  endif
                  ham(cnti,cntj) = ham(cnti,cntj) + 
     1                                     h01(i1,i2)*dj1j2*dk1k2 +
     2                                     h02(j1,j2)*di1i2*dk1k2 +
     3                                     h03(k1,k2)*di1i2*dj1j2
 700           continue
 600        continue   
         endif
         do 800 i=1,n
            do 900 j=1,i
               ham(j,i) = ham(i,j)
 900        continue
 800     continue
         if(prnt) then
            title='hamiltonian-'//title(1:2) 
            call prntrm(title,ham,n,n,n,n,iout)
         endif            
         if(op) then
            write(iout,1) title
            call rdiag(ham,cham,eig,ceig,vec,cvecl,cvecr,
     1                 work,lwork,rwork,n,n,mattyp)
         endif
      endif
      return         
    1 format(/,5x,'diagonalize time-independent ',a3,' dimensional '
     1            'schroedinger equation')
      end       
