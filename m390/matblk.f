*deck @(#)matblk.f	1.1  8/1/91
      subroutine matblk(mat,blksiz,nblks,iout)
c
c***begin prologue     matblk
c***date written       871125   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c
c***keywords           blocked matrix printing
c***author             saxe, paul (lanl)
c***source             @(#)matblk.f	1.1   8/1/91
c
c***purpose            to print a blocked matrix, block by block.
c
c***description        
c
c***references         
c
c***routines called    (none)
c
c***end prologue       matblk
c
      implicit integer (a-z)
c
      real*8 mat(*)
      integer blksiz(nblks)
      integer iout
c
c     ----- print each block if it exists -----
c
      pt=1
      do 2 blk=1,nblks
         size=blksiz(blk)
         if (size.gt.0) then
            write (iout,1) blk
 1          format (/,t10,'block ',i3)
            call matout(mat(pt),size,size,size,size,iout)
            pt=pt+size**2
         end if
 2    continue
c
c
      return
      end
