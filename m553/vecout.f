*deck @(#)vecout.f	5.1  11/6/94
      subroutine vecout(vec,nbf,nob)
c
c***begin prologue
c***date written       871022   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c
c***keywords
c***author             lengsfield, byron (brl)
c***source             @(#)vecout.f	5.1   11/6/94
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
      real*8 vec(nbf,2)
c
      common /io/ inp,iout
c
c
      write(iout,2)
 2    format(//,'  vector output ')
      do 20 i=1,nob,6
         lim=min(i+6-1,nob)
         write(iout,1)(k,k=i,lim)
 1       format(/,6(4x,i3,4x))
         write(iout,3)
 3       format(/,6(1x,'---------',1x))
         do 10 j=1,nbf
            write(iout,50)(vec(j,m),m=i,lim)
 10      continue
 20   continue
c
 50   format(6(1x,f10.6))
c
      return
      end
