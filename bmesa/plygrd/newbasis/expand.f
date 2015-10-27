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
      subroutine expand(z,p,pn,q,wt,qm,wtm,phrse,titphr,
     1                  npt,nmax,type,maxgrd,maxd,
     2                  modfed,prnt)
      implicit integer (a-z)
      real*8 z, tmp, dum
      character*(*) phrse, titphr, type
      character*80 title
      logical modfed, prnt
      dimension z(*)
      dimension p(maxgrd,maxgrd), pn(maxgrd,maxgrd), q(*), wt(*)
      dimension qm(*), wtm(*), npt(*), nmax(*), titphr(2)
      common/io/inp, iout
      pointer(p1,tmp(1))
      call iosys('read integer "no. subgrids for '
     1            //phrse//'" from bec',1,ngrid,0,' ')
      call iosys('read integer "no. points for '
     1            //phrse//'" from bec',ngrid,npt,0,' ')
      call iosys('read integer "mod. no. points for '
     1            //phrse//'" from bec',ngrid,nmax,0,' ')
      call iosys('read character "'//titphr(2)//'" from bec',
     1            -1,0,0,title)
      call iosys('read integer '//title//' from bec without '//
     1           'rewinding',ngrid,q,0,' ')
      call iosys('read integer '//title//' from bec without '//
     1           'rewinding',ngrid,wt,0,' ')
      call iosys('read integer '//title//' from bec without '//
     1           'rewinding',ngrid,qm,0,' ')
      call iosys('read integer '//title//' from bec without '//
     1           'rewinding',ngrid,wtm,0,' ')
      call iosys('read integer '//title//' from bec without '//
     1           'rewinding',maxgrd*maxgrd,p,0,' ')
      offset = 2*maxgrd*maxgrd
      call iosys('read integer '//title//' from bec without '//
     1           'rewinding',maxgrd*maxgrd,pn,offset,' ')
      call iosys('rewind '//title//' on bec read-and-write',
     1            0,0,0,' ')
      call iosys('read character "'//titphr(1)//'" from bec',0,0,
     1            0,title)
      do 10 i=1,ngrid
         call iosys('read real '//title//' from bec',npt(i),
     1               z(q(i)),q(i)-1,' ')
         call iosys('read real '//title//' from bec',npt(i),
     1               z(wt(i)),wt(i)-1,' ')
         call iosys('read real '//title//' from bec',npt(i),
     1               z(qm(i)),qm(i)-1,' ')
         call iosys('read real '//title//' from bec',npt(i),
     1               z(wtm(i)),wtm(i)-1,' ')
 10   continue              
      fine=1
      corse=fine+maxd
      coef=corse+maxd
      need=wpadti(coef+maxd)
      call memory(need,p1,ngot,'expand',0)
      write(iout,*) type
c
c     just do the simpleminded interpolation at the coarse grid
c     using the fine grid functions.
c
      do 20 i=1,ngrid
         if(.not.modfed) then
            call iosys('read real '//title//' from bec',npt(i)*npt(i),
     1                  z(p(i,i)),p(i,i)-1,' ')
            call cfine(tmp(fine),z(p(i,i)),z(q(i)),z(wt(i)),type,
     1                 npt(i),prnt)
         else
            call iosys('read real '//title//' from bec',nmax(i)*nmax(i),
     1                  z(pn(i,i)),pn(i,i)-1,' ')
            call cfine(tmp(fine),z(pn(i,i)),z(qm(i)),z(wtm(i)),
     1                 type,nmax(i),prnt)
         endif
         do 30 j=1,ngrid
            if(.not.modfed) then
               call iosys('read real '//title//' from bec',
     1                     npt(i)*npt(j),z(p(j,i)),p(j,i)-1,' ')
               call iosys('read real '//title//' from bec',
     1                     npt(j)*npt(j),z(p(j,j)),p(j,j)-1,' ')
               write(iout,1) npt(i), npt(j)
               call gi2gj(tmp(fine),tmp(corse),z(p(j,i)),dum,dum,dum,
     1                    npt(j),ndum,ndum,ndum,npt(i),ndum,ndum,ndum,
     2                    npt(j),npt(i),1,1)
               call cmpare(tmp(corse),tmp(coef),z(p(j,j)),dum,dum,
     1                     dum,z(q(j)),dum,dum,dum,npt(j),ndum,ndum,
     2                     ndum,1,npt(j),1,type,prnt)
            else
               call iosys('read real '//title//' from bec',
     1                     nmax(i)*nmax(j),z(pn(j,i)),pn(j,i)-1,' ')
               call iosys('read real '//title//' from bec',
     1                     nmax(j)*nmax(j),z(pn(j,j)),pn(j,j)-1,' ')
               write(iout,1) nmax(i), nmax(j)
               call gi2gj(tmp(fine),tmp(corse),z(pn(j,i)),dum,dum,dum,
     1                    nmax(j),ndum,ndum,ndum,nmax(i),ndum,ndum,ndum,
     3                    nmax(j),nmax(i),1,1)
               call cmpare(tmp(corse),tmp(coef),z(pn(j,j)),dum,dum,
     1                     dum,z(qm(j)),dum,dum,dum,nmax(j),ndum,ndum,
     2                     ndum,1,nmax(j),1,type,prnt)
            endif
 30      continue
 20   continue   
      call memory(-ngot,p1,idum,'expand',idum)
      return
 1    format(/,5x,'interpolating from grid = ',i4,' points',/,5x,
     1            '               to',/,5x
     2            '                   grid = ',i4,' points')          
      end
