*deck blmat.f
c***begin prologue     blmat
c***date written       951229   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c***keywords           
c***author             schneider, barry (nsf)
c***source             
c***purpose            lobatto blmat matrix elements
c***                   
c***description        computes the blmat matrix elements
c                            
c                      M = - 1./( 2 * mass ) [ fa(xleft) * fb'(xleft) 
c                          - fa(xright)  * fb'(xright) ]
c***references         
c
c***routines called    
c***end prologue       blmat
      subroutine blmat(bl,fa,dfb,pt,mass,coord,parity,n,prn,sec)
      implicit integer (a-z)
      real*8 bl, fa, dfb, pt, mass, scale
      character*(*) coord, parity
      logical prn
      character*80 title
      character*3 itoc
      dimension bl(n,n), fa(n,n), dfb(n,n), pt(n) 
      common/io/inp, iout
      call rzero(bl,n*n)
      scale = - .5d0/mass
      if(coord.ne.'rho') then
         do 10 i=1,n
            bl(n,i) = bl(n,i) - fa(n,n)*dfb(n,i)
            bl(1,i) = bl(1,i) + fa(1,1)*dfb(1,i)
 10      continue   
      else
         if(parity.eq.'even') then
            do 20 i=1,n
               bl(n,i) = bl(n,i) - 2.d0*pt(n)*pt(n)*fa(n,n)*dfb(n,i)
               bl(1,i) = bl(1,i) + 2.d0*pt(1)*pt(1)*fa(1,1)*dfb(1,i)
 20         continue
         elseif(parity.eq.'odd'.or.parity.eq.'none') then 
            do 30 i=1,n
               bl(n,i) = bl(n,i) - pt(n)*fa(n,n)*dfb(n,i)
               bl(1,i) = bl(1,i) + pt(1)*fa(1,1)*dfb(1,i)
 30         continue
         endif
      endif
      call vscale(bl,bl,scale,n*n)   
      if(prn) then
         title='unnormalized Bloch matrix for sector = '//itoc(sec)
         call prntrm(title,bl,n,n,n,n,iout)
      endif
      return
      end       



