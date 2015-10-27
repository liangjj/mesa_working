*deck ssseq
      subroutine ssseq(maxap3,natoms,nop,maxop,iatflg,nperm)
      implicit real*8(a-h,o-z)
c
c     flag (set iatflg to 1) the atoms which are in symmetric
c     subspaces equivalent to the current one.
c
      dimension iatflg(1), nperm(maxap3,maxop)
c
c     call rtrace(6hssseq ,1)
      if (nop .eq. 1) return
      do 100 iat=1,natoms
         if (iatflg(iat) .ne. 2) goto 100
         do 50 iop=2,nop
            do 50 jat=1,natoms
               if (nperm(jat,iop) .ne. iat) goto 50
               if (iatflg(jat) .eq. 0) iatflg(jat) = 1
 50         continue
 100     continue
         return
         end
