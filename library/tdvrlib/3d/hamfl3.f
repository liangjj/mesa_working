*deck hamfl3.f
c***begin prologue     hamfl3
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
c                      [ H11  H12 ]    H11 = - hx - vxt = H22
c                      [          ]    H12 = - d/dt
c                      [ H21  H22 ]    H21 =   d/dt
c                                                                  
c***references         
c
c***routines called    
c***end prologue       hamfl3
      subroutine hamfl3(h11,h12,h21,h22,h1,h2,h3,ht,v,ind,
     1                  nd,nt,n,nc,m,prn)
      implicit integer (a-z)
      real*8 h11, h12, h21, h22, h1, h2, h3, ht, v
      logical prn
      dimension nd(3)
      dimension h11(m,*), h12(m,*), h22(m,*), h21(m,*)
      dimension ind(nd(3),nd(2),nd(1),nt,nc)
      dimension h1(nd(1),nd(1)), h2(nd(2),nd(2)), h3(nd(3),nd(3))
      dimension ht(nt,nt), v(n,nc,nc)
      common/io/inp, iout
      do 10 ci=1,nc
         do 20 ti=1,nt
            do 30 xi=1,nd(1)
               do 40 yi=1,nd(2)
                  do 50 zi=1,nd(3)
                     zyxtci=ind(zi,yi,xi,ti,ci)
                     do 60 zj=1,nd(3)
                        zyxtcj=ind(zj,yi,xi,ti,ci)
                        h11(zyxtci,zyxtcj) = h11(zyxtci,zyxtcj) 
     1                                            - 
     2                                       h3(zi,zj)
                        h22(zyxtci,zyxtcj) = h22(zyxtci,zyxtcj) 
     1                                            - 
     2                                       h3(zi,zj)
 60                  continue   
                     do 70 yj=1,nd(2)
                        zyxtcj=ind(zi,yj,xi,ti,ci)
                        h11(zyxtci,zyxtcj) = h11(zyxtci,zyxtcj) 
     1                                            - 
     2                                       h2(yi,yj)
                        h22(zyxtci,zyxtcj) = h22(zyxtci,zyxtcj) 
     1                                            - 
     2                                       h2(yi,yj)
 70                  continue
                     do 80 xj=1,nd(1)
                        zyxtcj=ind(zi,yi,xj,ti,ci)
                        h11(zyxtci,zyxtcj) = h11(zyxtci,zyxtcj) 
     1                                            - 
     2                                       h1(xi,xj)
                        h22(zyxtci,zyxtcj) = h22(zyxtci,zyxtcj) 
     1                                            - 
     2                                       h1(xi,xj)   
 80                  continue   
 50               continue
 40            continue   
 30         continue   
 20      continue
 10   continue   
      do 100 ci=1,nc
         do 110 ti=1,nt
            do 120 xi=1,nd(1)
               do 130 yi=1,nd(2)
                  do 140 zi=1,nd(3)
                     zyxtci=ind(zi,yi,xi,ti,ci)
                     do 150 tj=1,nt
                        zyxtcj=ind(zi,yi,xi,tj,ci)
                        h12(zyxtci,zyxtcj) = h12(zyxtci,zyxtcj) 
     1                                            - 
     2                                       ht(ti,tj)
                        h21(zyxtci,zyxtcj) = h21(zyxtci,zyxtcj) 
     1                                            - 
     2                                       ht(ti,tj)
 150                 continue
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
                     do 250 zi=1,nd(3)
                        zyxtci=ind(zi,yi,xi,ti,ci)
                        zyxtcj=ind(zi,yi,xi,ti,cj)
                        h11(zyxtci,zyxtcj) = h11(zyxtci,zyxtcj) 
     1                                               - 
     2                                       v(vcnt,ci,cj)   
                        h22(zyxtci,zyxtcj) = h22(zyxtci,zyxtcj) 
     1                                               - 
     2                                       v(vcnt,ci,cj)   
 250                 continue
 240              continue   
 230           continue   
 220        continue   
 210     continue
 200  continue   
      return
      end       



