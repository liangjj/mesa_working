*deck %W%  %G%
      subroutine mcprti (nsym,nbf,nob,ntab,iic)
c
c***begin prologue
c***date written       871022   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c
c***keywords
c***author             lengsfield, byron (brl)
c***source             %W%   %G%
c
c***purpose
c
c***description
c
c***references
c
c***routines called    (none)
c
c***end prologue
c
      implicit real*8 (a-h,o-z)
cc
cmp   extended dummy iic
cc
      dimension nbf(2),nob(2),ntab(2),iic(2)
c
      common /io/ inp,iout
c
c
c     lll=ntab(nsym+1)
c     write(iout,81) lll
c     write(iout,82) (iic(mmm),mmm=1,lll)
c  81 format('  mcprti  lll ',i4,/,' icc ')
c  82 format(8(1x,i3))
c
      write (iout,83)
 83   format('0orbital mixings')
c
 10   ne=0
      do 40 l=1,nsym
         if (nbf(l).eq.0) go to 40
         if (nob(l).eq.0) go to 40
         ia=ntab(l)+1
         na=nob(l)
         nb=nbf(l)
c
         write (iout,200) l
         do 30 n = 1, na, 10
            nc=n + 9
            if (nc.gt.na) nc=na
            nd=nc-n+1
            write (iout,202) (i,i=n,nc)
            ib=ia+nd*nb-1
            write (iout,206)
            do 20 k=1,nb
               write (iout,208) k,(iic(i),i=ia,ib,nb)
 20         ia=ia+1
 30      ia=ib+1
         ne=ne+nob(l)
 40   continue
      return
c
 200  format ('0symmetry',i3)
 202  format ('0',9x,10i6)
 206  format (1x)
 208  format (1x,i9,8(1x,i5))
c
      end
