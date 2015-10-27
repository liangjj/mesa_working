*deck @(#)nrmylm.f	1.2  10/27/94
c***begin prologue     nrmylm
c***date written       921223   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c***keywords           test norm of spherical harmonics
c***author             schneider, barry (nsf)
c***source             
c***references         none
c
c***routines called
c***end prologue       nrmylm
      subroutine nrmylm (plm,phifn,ylm,wtthe,wtphi,wtang,ovlp,nthet,
     1                   nphi,nang,lval,mval,totplm,totylm,
     2                   totphm,nonsep)
      implicit integer (a-z)
      real *8 plm, phifn, ylm, wtthe, wtphi, wtang
      real*8 ovlp
      logical nonsep
      dimension plm(nthet,totplm), phifn(nphi,totphm), ylm(nang,totylm)
      dimension wtthe(nthet,2), wtphi(nphi,2), wtang(nang,2), ovlp(*)
      common /io/ inp, iout
c----------------------------------------------------------------------c     
c             what we do depends on whether the angular integration    c
c             is separable or non-separable in (theta,phi)             c
c----------------------------------------------------------------------c
      if (nonsep) then
c      
c              calculate square root of weights
c
          call vsqrt(wtang(1,2),wtang(1,2),nang)
c          
c              scale ylm by square root of weights to vectorize
c 
c                           integrals
          call vmmul(wtang(1,2),ylm,ylm,nang,totylm)          
          locylm=1
          do 10 m=0,mval
             nm=2
             if (m.eq.0) then
                 nm=1
             endif        
             locyln=1
             do 20 n=0,m
                nn=2
                if (n.eq.0) then
                    nn=1
                endif
                if (m.eq.n) then
                    call ylmint(ylm(1,locylm),ylm(1,locylm),ovlp,
     1                          lval,m,m,nm,nm,nang)
                else
                    call ylmint(ylm(1,locylm),ylm(1,locyln),ovlp,
     1                          lval,m,n,nm,nn,nang)
                endif
                locyln=locyln+(lval-n+1)*nn
   20        continue
             locylm=locylm+(lval-m+1)*nm                                
   10     continue
c   
c         take inverse of square root of weights and rescale ylm
c         returning them to their previous values. copy back original
c         weights into wtang
c  
          call vinv(wtang(1,2),wtang(1,2),nang)
          call vmmul(wtang(1,2),ylm,ylm,nang,totylm)
          call copy(wtang(1,1),wtang(1,2),nang)
      else
c      
c         same procedure as above for getting square root of weights
c         and scaling but done for each angle separately
c
          call vsqrt(wtthe(1,2),wtthe(1,2),nthet)
          call vsqrt(wtphi(1,2),wtphi(1,2),nphi)
          call vmmul(wtthe(1,2),plm,plm,nthet,totplm)
          call vmmul(wtphi(1,2),phifn,phifn,nphi,totphm)
          locplm=1
          locphm=1
          do 30 m=0,mval
             nm=2
             if (m.eq.0) then
                 nm=1
             endif   
             call legint(plm(1,locplm),phifn(1,locphm),ovlp,nthet,nphi,
     1                   lval,m,nm)
             locplm=locplm+(lval-m+1)
             locphm=locphm+nm
   30     continue
c         same procedure as above for rescaling
          call vinv(wtthe(1,2),wtthe(1,2),nthet)
          call vinv(wtphi(1,2),wtphi(1,2),nphi)
          call vmmul(wtthe(1,2),plm,plm,nthet,totplm)
          call vmmul(wtphi(1,2),phifn,phifn,nphi,totphm)
          call copy(wtthe(1,1),wtthe(1,2),nthet)
          call copy(wtphi(1,1),wtphi(1,2),nphi)
      endif    
      return
c
      end


