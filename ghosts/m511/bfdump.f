      subroutine bfdump(nat,nprim,ncont,ntypes,c,ex,cont,ptprim,noprim,
     $     nocont,ptcont,iout)
      integer nat,nprim,ncont,ntypes,iout
      real*8 c(3,nat),ex(nprim),cont(ncont)
      integer ptprim(nat,ntypes),noprim(nat,ntypes),ptcont(nat,ntypes), 
     $     nocont(nat,ntypes)
      write(iout,*)"Basis set info:"
      write(iout,*)
      write(iout,*)' atomic coordinates'
      do 6661 i=1,nat
         write(iout,*)'   atom i',c(1,i),c(2,i),c(3,i)
 6661 continue 
      write(iout,*)' exponents:'
      do 6662 i=1,nprim
         write(iout,*)'   ',i,ex(i)
 6662 continue 
      write(iout,*)' contraction coefficients:'
      do 6664 i=1,ncont
         write(iout,*)'   ',i,cont(i)
 6664 continue 
      write(iout,*)' pointers:'
      write(iout,*)'   ptprim,noprim,nocont,ptcont:'
      do 6665 i=1,nat
         write(iout,*)'   atom ',i
         do 6666 j=1,ntypes
            write(iout,*)'     ',j,ptprim(i,j),noprim(i,j),
     $           nocont(i,j),ptcont(i,j)
 6666    continue 
 6665 continue 
      return
      end
