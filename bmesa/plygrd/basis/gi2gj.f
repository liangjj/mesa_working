*deck gi2gj.f
c***begin prologue     gi2gj
c***date written       951229   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c***keywords           polynomials
c***author             schneider, barry (nsf)
c***source             
c***purpose            interpolate a vector from a coarse to a fine grid
c***                   
c***                                        
c***description        the initial vector, vi, contains the expansion 
c***                   coefficients of the three dimensional object in 
c***                   terms of dvr basis set (grid) i.  the final vector
c***                   vj is the original vector interpolated onto 
c***                   basis set (grid) j.  basis set j may be either a
c***                   finer or a coarser basis
c***references         
c
c***routines called    
c***end prologue       gi2gj
      subroutine gi2gj(vi,vj,p1ji,p2ji,p3ji,p4ji,n1j,n2j,n3j,n4j,
     1                 n1i,n2i,n3i,n4i,nj,ni,nvc,dim)
      implicit integer (a-z)
      real*8 vi, vj, p1ji, p2ji, p3ji, p4ji
      real*8  sc
      dimension vi(ni,2), vj(nj,2)
      dimension p1ji(n1j,n1i), p2ji(n2j,n2i)
      dimension p3ji(n3j,n3i), p4ji(n4j,n4i)
      common/io/inp, iout
      pointer(p,scr(1))
      if(dim.eq.1) then
         call i2j1d(p1ji,vi,vj,n1j,n1i,nvc)
      elseif(dim.eq.2) then
         need=wptoin(n1i*n2j*2)
         call memory(need,p,nscr,'gi2gj',0)
         call i2j2d(p1ji,p2ji,vi,vj,scr,n1j,n2j,n1i,n2i,nvc)
         call memory(-nscr,p,idum,'gi2gj',idum)
      elseif(dim.eq.3) then
         need=wptoin(2*n1i*n2j*n3j)
         call memory(need,p,nscr,'gi2gj',0)
         call i2j3d(p1ji,p2ji,p3ji,vi,vj,vj,scr,n1j,n2j,n3j,
     1              n1i,n2i,n3i,nvc,type)
         call memory(-nscr,p,idum,'gi2gj',idum)
      elseif(dim.eq.4) then
         need=wptoin(n4j*n3j*n2j*n1i*2)
         call memory(need,p,nscr,'gi2gj',0)
         call i2j4d(p1ji,p2ji,p3ji,p4ji,vi,vj,scr,vf,scr,n1j,n2j,n3j,
     1              n4j,n1i,n2i,n3i,n4i,nvc)
         call memory(-nscr,p,idum,'gi2gj',idum)
      else
         call lnkerr('error in dimension')
      endif
      return
      end









