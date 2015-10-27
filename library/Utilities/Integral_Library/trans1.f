*deck @(#)trans1.f	5.1  11/6/94
      subroutine trans1(orig,final,ia,ib,fa,fb,a,b,t1,len1,lenblk,
     $                  imin,imax,jmin,jmax,ncart)
c***begin prologue     trans1
c***date written       850601  (yymmdd)
c***revision date      yymmdd  (yymmdd)
c***keywords           one-electron, integrals, trnasformation
c***author             saxe, paul (lanl)
c***source
c***purpose                                                     t
c                      vectorized one-index transformation: v'=c vc
c***description
c                      call trans1(orig,final,ia,ib,fa,fb,a,b,t1,len1,lenblk,
c                                  imin,imax,jmin,jmax,ncart)
c
c***references
c***routines called    ebtc(math),ebc(math)
c***end prologue       trans1
      implicit integer (a-z)
c
      real*8 orig(ia,ib,lenblk),final(fa,fb,lenblk)
      real*8 a(ia,fa,imin:imax),b(ib,fb,jmin:jmax),t1(len1)
      integer ncart(0:*)
c
c     ----- start timing -----
c
c
      ij=0
      do 4 imom=imin,imax
         do 3 i=1,ncart(imom)
            do 2 jmom=jmin,jmax
               do 1 j=1,ncart(jmom)
                  ij=ij+1
                  call ebtc(t1,a(1,1,imom),orig(1,1,ij),fa,ia,ib)
                  call ebc(final(1,1,ij),t1,b(1,1,jmom),fa,ib,fb)
    1          continue
    2       continue
    3    continue
    4 continue
c
c     ----- stop timing -----
c
c
      return
      end
