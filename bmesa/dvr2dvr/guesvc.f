*deck guesvc.f
c***begin prologue     guesvc
c***date written       960723   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c***keywords           basis, transformation
c***author             schneider, barry (nsf)
c***source             
c***purpose            transformation of vectors.
c***                   
c***references         
c
c***routines called    
c***end prologue       guesvc
      subroutine guesvc(veci,tp2pn,rhs,veco,smat,ns,nrhs,ntrial,nb,
     1                  prnt,drctv,grid)
      implicit integer (a-z)
      real*8 veci, tp2pn, veco, rhs, smat
      dimension veci(ns,*), tp2pn(nb,nb), veco(nb,*), rhs(nb,*), smat(*)
      character*2 itoc
      character*80 title
      character*(*) drctv
      logical prnt
      common/io/inp, iout
      if(grid.eq.1) then
         call copy(rhs,veco,nb*nrhs)
         ntrial=nrhs
      else
         if(drctv.eq.'all-vectors') then
            call iosys('read integer "number of vectors for grid '
     1                  //itoc(grid-1)//'" from ham',1,ntrial,0,' ')     
            call iosys('read real "vectors for grid '//itoc(grid-1)
     1                 //'" from ham',ns*ntrial,veci,0,' ')
            if(prnt) then
               title='overlap matrix for basis vectors for '//
     1               'grid = '//itoc(grid-1)         
               call ebtcxx(smat,veci,veci,ntrial,ns,ntrial,ntrial,ns,ns)
               call prntrm(title,smat,ntrial,ntrial,ntrial,ntrial,iout) 
            endif
            call iosys('read real "transformation matrix for grid '
     1                 //itoc(grid)//'" from ham',nb*nb,tp2pn,0,' ')        
            call ebtcxx(veco,tp2pn,veci,nb,ns,ntrial,nb,nb,ns) 
            if(prnt) then
               call ebtcxx(smat,veco,veco,ntrial,nb,ntrial,ntrial,nb,nb)
               title='pseudo overlap matrix for basis vectors for '//
     1               'grid = '//itoc(grid-1)
               call prntrm(title,smat,ntrial,ntrial,ntrial,ntrial,iout)          
               title='basis vectors in new basis'
               call prntrm(title,veco,nb,ntrial,nb,ntrial,iout)
            endif
         elseif(drctv.eq.'solutions') then
            ntrial=nrhs
            call iosys('read real "solutions for grid '//itoc(grid-1)
     1                 //'" from ham',ns*ntrial,veci,0,' ')
            call iosys('read real "transformation matrix for grid '
     1                 //itoc(grid)//'" from ham',nb*nb,tp2pn,0,' ')        
            call ebtcxx(veco,tp2pn,veci,nb,ns,ntrial,nb,nb,ns) 
            ntrial=nrhs
            title='basis vectors in new basis'
            call prntrm(title,veco,nb,ntrial,nb,ntrial,iout)
         endif                             
      endif
      return
      end       
