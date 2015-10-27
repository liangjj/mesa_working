*deck @(#)aintch.f	5.1  11/6/94
      subroutine aintch(am,nv,i,i0,ii,li)
      implicit none
c
      integer nv,i,i0,ii,li
      real*8 am(1)
c
      integer int,j,ili,li0,ij,lij,jli
      real*8 work
c     
      li=int(am(li))
      if (li.ne.i) then
   15    ili=i0+li
         li0=li*nv-(li*(li-1))/2
         do 11 j=i,nv
            ij=i0+j
            work=am(ij)
            if(li.le.j) then
               lij=li0+j
               am(ij)=am(lij)
               am(lij)=work
            else
               jli=j*nv+li-(j*(j-1))/2
               am(ij)=am(jli)
               am(jli)=work
            endif
   11    continue
         work=am(ii)
         am(ii)=am(ili)
         am(ili)=work
      endif
c
c
      return
      end
