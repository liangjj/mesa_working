*deck expand.f
c***begin prologue     expand
c***date written       951229   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c***keywords           polynomials
c***author             schneider, barry (nsf)
c***source             
c***purpose            
c***                   
c***                                                          
c***references         
c
c***routines called    
c***end prologue       expand
      subroutine expand(y,z,iloc,type,ngrid,q,wt,npt,maxd,dim)
      implicit integer (a-z)
      real*8 y, z, tmp
      character*(*) type
      character*80 title
      dimension y(*), z(*)
      dimension iloc(dim,dim), q(dim), wt(dim), npt(dim)
      common/io/inp, iout
      pointer(p,tmp(1))
      fine=1
      corse=fine+maxd
      coef=corse+maxd
      need=wpadti(coef+maxd)
      call memory(need,p,ngot,'expand',0)
      write(iout,*) type
c
c     just do the simpleminded interpolation at the coarse grid
c     using the fine grid functions.
c
      do 10 i=1,ngrid
         pff=iloc(i,i)
         call cfine(tmp(fine),y(pff),z(q(i)),z(wt(i)),type,npt(i))
         do 20 j=1,ngrid
c            if(npt(i).ge.npt(j)) then       
               pcf=iloc(j,i)
               write(iout,1) npt(i), npt(j)
               call fin2cor(tmp(fine),tmp(corse),y(pcf),z(q(j)),type,
     1                      npt(i),npt(j))     
                pcc=iloc(j,j)
               call cor2cor(tmp(corse),tmp(coef),y(pcc),npt(j))  
c            endif
 20      continue
 10   continue   
      call memory(-ngot,p,idum,'expand',idum)
      return
 1    format(/,5x,'interpolating from grid = ',i4,' points',/,5x,
     1            '               to',/,5x
     2            '                   grid = ',i4,' points')          
      end
