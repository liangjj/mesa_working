*deck @(#)fixrep.f	4.1  7/7/93
 
      subroutine fixrep(t)
      implicit real*8(a-h,o-z)
c
c     the symmetry operations have been passed to the orbital
c     symmetry routines as if the molecule were re-oriented by the
c     symmetry program.  in certain cases, however, after this has
c     been done, the symmetry program changes its mind, and decides
c     not to re-orient the molecule after all.  when this happens,
c     this routine should be called to transform the matrices that
c     were passed to the orbital assignment routines.  the only
c     calling argument, "t", is the 3x3 matrix that would have been
c     used to re-orient the molecule.
c     r. a. whiteside - august 1979
c     m. j. frisch - april 1983
c
      dimension t(3,3),c(9)
      common/repcom/nsymop,nreps,lblrep(32),chrtbl(10,16),symops(9,10)
c
c     call rtrace(6hfixrep,1)
      if(nsymop.eq.0) return
      do 100 iop=1,nsymop
         call mul3x3(symops(1,iop),t,c)
         do 100 i=1,9
 100        symops(i,iop)=c(i)
            do 200 i=2,3
               im1=i-1
               do 200 j=1,im1
                  temp=t(i,j)
                  t(i,j)=t(j,i)
 200              t(j,i)=temp
                  do 300 iop=1,nsymop
                     call mul3x3(t,symops(1,iop),c)
                     do 300 i=1,9
 300                    symops(i,iop)=c(i)
                        return
                        end
