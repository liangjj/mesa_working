*deck rstrct
c***begin prologue     rstrct
c***date written       941208   (yymmdd)
c***revision date               (yymmdd)
c***keywords           m1200, link 1200, multigrid
c***author             schneider, b. i.(nsf)
c***source             m1200
c***purpose            restrict fine grid to coarse grid solution
c***references        
c
c***routines called    iosys, util and mdutil
c***end prologue       rstrct
      subroutine rstrct(uc,uf,nc,nf,type)
c
      implicit integer (a-z)
      real*8 uc, uf
      character*80 title
      character*(*) type
      dimension uc(nc,nc), uf(nf,nf)
      common/io/inp, iout
      if (type.eq.'full-weighting') then
          do 10 jc=2,nc-1
             jf=2*jc-1
             do 20 ic=2,nc-1
                if=2*ic-1
                uc(ic,jc)=.25d0*uf(if,jf) +
     1                    .125d0*(  uf(if+1,jf) + uf(if-1,jf) 
     2                                          + uf(if,jf+1)
     3                                          + uf(if,jf-1) )   +
     4                    .0625d0*( uf(if+1,jf+1) + uf(if+1,jf-1) +
     5                              uf(if-1,jf+1) + uf(if-1,jf-1) )    
 20          continue   
 10       continue
       elseif (type.eq.'half-weighting') then
          do 30 jc=2,nc-1
             jf=2*jc-1
             do 40 ic=2,nc-1
                if=2*ic-1
                uc(ic,jc)=.5d0*uf(if,jf) +
     1                    .125d0*( uf(if+1,jf) + uf(if-1,jf) 
     2                                         + uf(if,jf+1)
     2                                         + uf(if,jf-1) )
 40          continue   
 30       continue
      elseif(type.eq.'injection') then
          do 50 jc=2,nc-1
             jf=2*jc-1
             do 60 ic=2,nc-1
                if=2*ic-1
                uc(ic,jc)=uf(if,jf)
 60          continue   
 50       continue
      else
          call lnkerr('error in restriction operator type')
      endif
      nfx=2*nc-1                              
      do 70 ic=1,nc
c        lower x-axis edge      
         uc(ic,1)=uf(2*ic-1,1)
c        upper x-axis edge         
         uc(ic,nc)=uf(2*ic-1,nfx)
 70   continue
      do 80 jc=1,nc
c        left y-axis edge      
         uc(1,jc)=uf(1,2*jc-1)
c        right y-axis edge         
         uc(nc,jc)=uf(nfx,2*jc-1)
 80   continue
c      title='coarse solution after injection'
c      call prntrm(title,uc,nc,nc,nc,nc,iout)
      return
      end





