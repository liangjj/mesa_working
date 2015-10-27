      subroutine mcprtv (nsym,nbf,nob,ntab,c)
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
c     implicit real*8 (a-h,o-z)
cc
cmp   extended dummy c
cc
      dimension nbf(2),nob(2),ntab(2),c(2)
c
 10   ne=0
      do 40 l=1,nsym
         if (nbf(l).eq.0) go to 40
         if (nob(l).eq.0) go to 40
         ia=ntab(l)
         na=nob(l)
         nb=nbf(l)
c
cc    write (iout,200) l
         do 30 n=1,na,8
            nc=n+7
            if (nc.gt.na) nc=na
            nd=nc-n+1
cc    write (iout,202) (i,i=n,nc)
            ib=ia+nd*nb-1
cc    write (iout,206)
            do 20 k=1,nb
cc    write (iout,208) k,(c(i),i=ia,ib,nb)
 20         ia=ia+1
 30      ia=ib+1
         ne=ne+nob(l)
 40   continue
      return
c
c 200 format ('0sym=',i3)
c 202 format ('0',9x,8i15)
c 204 format ('0',9x,8f15.8)
 206  format (1x)
 208  format (1x,i9,8f15.8)
c
      end
