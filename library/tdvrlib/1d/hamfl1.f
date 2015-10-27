*deck hamfl1.f
c***begin prologue     hamfl1
c***date written       951229   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c***keywords           two-dim
c***author             schneider, barry (nsf)
c***source             3-dim
c***purpose            form all or a sub-block of the two dimensional 
c***                   hamiltonian.
c***                   
c***                   
c                            Structure of Subroutine
c
c                      [ H11  H12 ]    H11 = - h1 - v = H22
c                      [          ]    H12 = - d/dt
c                      [ H21  H22 ]    H21 =   d/dt
c                                                                  
c***references         
c
c***routines called    
c***end prologue       hamfl1
      subroutine hamfl1(h11,h12,h21,h22,h1,ht,v,ind,n1,nt,
     1                  n,nc,m,prn)
      implicit integer (a-z)
      real*8 h11, h12, h21, h22, h1, ht, v
      character*80 title 
      logical prn
      dimension h11(m,*), h12(m,*), h22(m,*), h21(m,*)
      dimension ind(n1,nt,nc), h1(n1,n1), ht(nt,nt), v(n,nc,nc)
      common/io/inp, iout
c      title='space matrix'
c      call prntrm(title,h1,n1,n1,n1,n1,iout)
c      title='time matrix'
c      call prntrm(title,ht,nt,nt,nt,nt,iout)
      do 10 ci=1,nc
         do 20 ti=1,nt
            do 30 xi=1,n1
               xtci=ind(xi,ti,ci)
               do 40 xj=1,n1
                  xtcj=ind(xj,ti,ci)
                  h11(xtci,xtcj) = h11(xtci,xtcj) - h1(xi,xj)
                  h22(xtci,xtcj) = h22(xtci,xtcj) - h1(xi,xj)
 40            continue   
 30         continue   
 20      continue
 10   continue
      do 50 ci=1,nc   
         do 60 ti=1,nt
            do 70 xi=1,n1
               xtci=ind(xi,ti,ci)
               do 80 tj=1,nt
                  xtcj=ind(xi,tj,ci)
                  h12(xtci,xtcj) = h12(xtci,xtcj) - ht(ti,tj)
                  h21(xtci,xtcj) = h21(xtci,xtcj) + ht(ti,tj)
 80            continue
 70         continue   
 60      continue
 50   continue
      do 100 ci=1,nc
         do 200 cj=1,nc
            vcnt=0
            do 300 ti=1,nt
               do 400 xi=1,n1
                  xtci=ind(xi,ti,ci)
                  xtcj=ind(xi,ti,cj)
                  vcnt=vcnt+1
                  h11(xtci,xtcj) = h11(xtci,xtcj) 
     1                                     - 
     2                             v(vcnt,ci,cj)   
                  h22(xtci,xtcj) = h22(xtci,xtcj) 
     1                                     - 
     2                             v(vcnt,ci,cj)   
 400           continue   
 300        continue   
 200     continue
 100  continue   
      if(prn) then
         title='H11'
         call prntrm(title,h11,n,n,m,m,iout)
         title='H22'
         call prntrm(title,h22,n,n,m,m,iout)
         title='H12'
         call prntrm(title,h12,n,n,m,m,iout)
         title='H21'
         call prntrm(title,h21,n,n,m,m,iout)
      endif
      return
      end       

