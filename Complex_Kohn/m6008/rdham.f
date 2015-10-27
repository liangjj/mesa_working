*deck @(#)rdham.f	1.1 9/8/91
c***begin prologue     rdham
c***date written       880423   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c***keywords           rdham, link 1106, kohn variational
c***author             schneider, barry (lanl), rescigno, tom(llnl)  
c***source             m1106
c***purpose            read non-local optical potential
c***description        read in non-local optical potential matrices or
c***                   exchange matrices formed in m9000
c***              
c
c***references         schneider and rescigno, physical review
c
c***routines called    iosys, util and mdutil
c***end prologue       rdham
      subroutine rdham(hamin,hambb,haminc,hambbc,orblst,nchan,nbscat,
     1                 npvec,matbv,nl2,nsmall,energy,dimmo,dimc,
     2                 prnt,statex,opt)
      implicit integer(a-z)
      logical prnt, statex
      character *24 ftitl
      character *8 colt, rowt
      character *80 title
      character *16 fptoc
      character *3 opt
      real*8 energy, hambb, hamin, rowv
      complex *16 haminc, hambbc
      dimension hamin(npvec,npvec), hambb(matbv,matbv), nbscat(dimc)
      dimension orblst(dimmo,dimc), haminc(npvec,npvec)
      dimension hambbc(matbv,matbv)
      data title/'bound-bound optical potential'/
      common /io/ inp, iout
c----------------------------------------------------------------------c
c          read in optical potential over p-space eigenvectors         c
c          these correspond to a target eigenstate multiplied          c
c                       by a virtual orbital                           c
c----------------------------------------------------------------------c
      ftitl='v(pp)-'//fptoc(energy)
      if (statex) then
          ftitl='h(pp)-'//fptoc(energy)
      endif
      if (opt.ne.'yes') then
          call iosys ('read real '//ftitl//' from kohndt',npvec*npvec,
     1                 hamin,0,' ')
      else
          call iosys ('read real '//ftitl//' from kohndt',
     1                 2*npvec*npvec,haminc,0,' ')
      endif
c----------------------------------------------------------------------c
c       construct the proper channel optical potential by placing the  c
c              matrix elements in the proper slots in hambb            c
c----------------------------------------------------------------------c
      if (opt.ne.'yes') then
          icnt=0
          do 10 i=1,nchan
             do 20 i1=1,nbscat(i)
                ipnt=nl2*(i-1)+orblst(i1,i)-nsmall
                icnt=icnt+1
                jcnt=0
                do 30 j=1,i
                   do 40 j1=1,nbscat(j)
                      jcnt=jcnt+1   
                      jpnt=nl2*(j-1)+orblst(j1,j)-nsmall
                      hambb(icnt,jcnt)=hamin(ipnt,jpnt)
                      hambb(jcnt,icnt)=hambb(icnt,jcnt)
   40              continue
                   jadd=jadd+nbscat(j)
   30           continue
   20        continue
             iadd=iadd+nbscat(i)
   10     continue
          if (prnt) then
              title='bound-bound denominator matrix'
              rowv=-99.d0
              colv=-99
              call mprir (hambb,rowv,colv,matbv,matbv,matbv,matbv,title,
     1                    rowt,colt,iout)
          endif
      else
          icnt=0
          do 50 i=1,nchan
             do 60 i1=1,nbscat(i)
                ipnt=nl2*(i-1)+orblst(i1,i)-nsmall
                icnt=icnt+1
                jcnt=0
                do 70 j=1,i
                   do 80 j1=1,nbscat(j)
                      jcnt=jcnt+1   
                      jpnt=nl2*(j-1)+orblst(j1,j)-nsmall
                      hambbc(icnt,jcnt)=haminc(ipnt,jpnt)
                      hambbc(jcnt,icnt)=hambbc(icnt,jcnt)
   80              continue
                   jadd=jadd+nbscat(j)
   70           continue
   60        continue
             iadd=iadd+nbscat(i)
   50     continue
          if (prnt) then
              title='bound-bound denominator matrix'
              call prntcm (title,hambbc,matbv,matbv,matbv,matbv,iout)
          endif
      endif
      return
      end
