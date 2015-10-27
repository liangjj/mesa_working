*deck hamfl2.f
c***begin prologue     hamfl2
c***date written       951229   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c***keywords           two-dim
c***author             schneider, barry (nsf)
c***source             3-dim
c***purpose            form three dimensional hamiltonian.
c***                   
c***                   
c                            Structure of Subroutine
c
c                      [ H11  H12 ]    H11 = - h1 - h2 - v = H22
c                      [          ]    H12 = - d/dt
c                      [ H21  H22 ]    H21 =   d/dt
c                                                                  
c***references         
c
c***routines called    
c***end prologue       hamfl2
      subroutine hamfl2(h11,h12,h21,h22,h1,h2,ht,v,ind,nd,nt,
     1                  n,nc,m,prn)
      implicit integer (a-z)
      real*8 h11, h12, h21, h22, h1, h2, ht, v
      logical prn
      dimension nd(2)
      dimension h11(m,*), h12(m,*), h22(m,*), h21(m,*)
      dimension ind(nd(2),nd(1),nt,nc)
      dimension h1(nd(1),nd(1)), h2(nd(2),nd(2))
      dimension ht(nt,nt), v(n,nc,nc)
      common/io/inp, iout
      do 10 ci=1,nc
         do 20 ti=1,nt
            do 30 xi=1,nd(1)
               do 40 yi=1,nd(2)
                  yxtci=ind(yi,xi,ti,ci)
                  do 50 yj=1,nd(2)
                     yxtcj=ind(yj,xi,ti,ci)
                     h11(yxtci,yxtcj) = h11(yxtci,yxtcj) 
     1                                     - 
     2                                  h2(yi,yj)
                     h22(yxtci,yxtcj) = h22(yxtci,yxtcj) 
     1                                     - 
     2                                  h2(yi,yj)
 50               continue   
                  do 60 xj=1,nd(1)
                     yxtcj=ind(yi,xj,ti,ci)   
                     h11(yxtci,yxtcj) = h11(yxtci,yxtcj) 
     1                                      - 
     2                                  h1(xi,xj)
                     h22(yxtci,yxtcj) = h22(yxtci,yxtcj) 
     1                                      - 
     2                                  h1(xi,xj)
 60               continue   
 40            continue   
 30         continue
 20      continue   
 10   continue   
      do 100 ci=1,nc
         do 110 ti=1,nt
            do 120 xi=1,nd(1)
               do 130 yi=1,nd(2)
                  yxtci=ind(yi,xi,ti,ci)
                  do 140 tj=1,nt
                     yxtcj=ind(yi,xi,tj,ci) 
                     h12(yxtci,yxtcj) = h12(yxtci,yxtcj) 
     1                                      - 
     2                                  ht(ti,tj)
                     h21(yxtci,yxtcj) = h21(yxtci,yxtcj) 
     1                                      + 
     2                                  ht(ti,tj)
 140              continue   
 130           continue   
 120        continue
 110     continue
 100  continue
      do 200 ci=1,nc
         do 210 cj=1,nc
            vcnt=0
            do 220 ti=1,nt
               do 230 xi=1,nd(1)
                  do 240 yi=1,nd(2)
                     vcnt=vcnt+1
                     yxtci=ind(yi,xi,ti,ci)
                     yxtcj=ind(yi,xi,ti,cj)
                     h11(yxtci,yxtcj) = h11(yxtci,yxtcj) 
     1                                      - 
     2                                  v(vcnt,ci,cj)
                     h22(yxtci,yxtcj) = h22(yxtci,yxtcj) 
     1                                      - 
     2                                  v(vcnt,ci,cj)   
 240              continue
 230           continue   
 220        continue   
 210     continue
 200  continue   
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


