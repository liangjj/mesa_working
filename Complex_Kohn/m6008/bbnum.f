*deck @(#)bbnum.f	1.1 9/8/91
c***begin prologue     bbnum
c***date written       890524   (yymmdd)
c***revision date               (yymmdd)
c***keywords           kohn integrals
c***author             schneider barry (lanl)
c***source             m6008
c***purpose            read bound-bound one body integrals and
c***                   form final bound-bound matrix
c***references       
c
c***routines called
c***end prologue       bbnum  
      subroutine bbnum(hamin,hamnum,echan,zeroc,nchan,nmo,matbb,nbtot,
     1                 orblst,finlst,dimmo,dimc,prnt)
      implicit integer (a-z)
      real *8 hamin, hamnum, echan, rowv
      character *80 title
      character *8 colt, rowt
      character *4 itoc
      logical prnt, zeroc
      common /io/ inp, iout
      dimension hamin(nmo,nmo), hamnum(matbb,matbb)
      dimension nbtot(dimc), finlst(dimmo,dimc), orblst(dimmo,dimc)
      dimension echan(dimc), zeroc(dimc)
      call rzero(hamnum,matbb*matbb)
      do 3 i=1,nchan
         if (.not.zeroc(i)) then
             do 4 j=1,i
                if (.not.zeroc(j)) then
                     call iosys ('read real "mo kohn hamiltonian'//
     1                            itoc(i)//itoc(j)//'" from kohndt',
     2                            nmo*nmo,hamin,0,' ')
                     if ( i.eq.j) then
                         do 5 i1=1,nbtot(i)
                            icnt=finlst(i1,i)
                            ipnt=orblst(i1,i)
                            do 6 j1=1,nbtot(j)
                               jcnt=finlst(j1,j)
                               jpnt=orblst(j1,j)
                               hamnum(icnt,jcnt)=hamnum(icnt,jcnt)
     1                                           +hamin(ipnt,jpnt)
    6                       continue
    5                    continue
                     else
                         do 7 i1=1,nbtot(i)
                            icnt=finlst(i1,i)
                            ipnt=orblst(i1,i)
                            do 8 j1=1,nbtot(j)
                               jcnt=finlst(j1,j)
                               jpnt=orblst(j1,j)
                               hamnum(icnt,jcnt)=hamnum(icnt,jcnt)
     1                                           +hamin(ipnt,jpnt)
                               hamnum(jcnt,icnt)=hamnum(icnt,jcnt)
    8                       continue
    7                    continue
                     endif
                endif
    4        continue
          endif
    3 continue
c----------------------------------------------------------------------c
c             put echan(i) on diagonal it is needed later              c
c----------------------------------------------------------------------c
      do 1 i=1,nchan
         do 2 i1=1,nbtot(i)
              count=finlst(i1,i)
              hamnum(count,count)=hamnum(count,count)+echan(i)
    2    continue
    1 continue    
      if (prnt) then
          title='bound-bound numerator matrix'
          rowv=-99.
          colv=-99
          call mprir(hamnum,rowv,colv,matbb,matbb,matbb,matbb,
     1               title,rowt,colt,iout)
      endif
      return
      end
