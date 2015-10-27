*deck david.f
c***begin prologue     david
c***date written       951229   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c***keywords           eigenvalues, eigenvectors, davidson
c***author             schneider, barry (nsf)
c***source             
c***purpose            driver for davidson algorithm.
c***references         
c
c***routines called    
c***end prologue       david
      subroutine david(hambuf,eigbuf,vec,hvec,root,dvdmat,dvdvec,cvn,
     1                 scr,thresh,cnverg,ops,n,nroots,nattim,mxvc,iter,
     2                 prtflg)      
      implicit integer (a-z)
      real*8 hambuf, eigbuf
      real*8 vec, hvec, root, dvdvec, dvdmat, cvn, scr
      real*8 cnverg, thresh, fpkey, rjunk
      character*(*) ops
      character*20 status
      character*8 prtflg
      dimension eigbuf(n), hambuf(n,n), dvdmat(*), dvdvec(*), cvn(n)
      dimension vec(n,mxvc), hvec(n,mxvc), scr(n,nroots), root(nroots)
      common/io/inp, iout
      call rzero(scr,n*nroots)
      if (prtflg.eq.'minimum') then
         status='noprint'
      else
         status='print'
      end if 
      cntg=intkey(ops,'davidson=number-of-guess-vectors',nroots,' ')     
      gsize=intkey(ops,'davidson=size-of-guess-matrix',cntg,' ')
      call rzero(vec,n*mxvc)
      do 10 i=1,gsize
         do 20 j=1,gsize
            vec(i,j)=hambuf(i,j)
 20      continue
 10   continue
      call guess(vec,eigbuf,hvec,n,gsize,nroots)   
      do 30 i=1,n
         eigbuf(i)=hambuf(i,i)
         hambuf(i,i)=0.d0
 30   continue         
c
c      write(iout,*) ' initializing davdag '
      call davdag('initialize',status,eigbuf,thresh,dvdvec,scr,n,iter,
     1             nroots,iout,nattim,0.d0,cnverg,dvdvec,cvn)
      nvc=cntg   
      call honv(vec,hvec,hambuf,n,nvc)
c***********************************************************************
c                 c. check matrix size     
c      write(iout,*) 'check matrix size'
      call davdag('check matrix size',status,eigbuf,vec,hvec,scr,n,iter,
     1             nvc,junk,junk,rjunk,rjunk,rjunk,rjunk)
c***********************************************************************
c                 b. solve the small eigenvalue equation and test
c                    for convergence. for unconverged roots add
c                    new vectors using perturbation theory. the number
c                    of these vectors is returned as nresid and they
c                    reside in the scr array upon return.
        do while (status.ne.'converged')                    
c           write(iout,*) 'solving small matrix '
           call davdag('solve',status,eigbuf,vec,hvec,scr,n,iter,nresid,
     1                  junk,junk,dvdmat,root,dvdvec,cvn)
           if(status.ne.'converged') then
              test=nvc+nresid
              if(test.gt.mxvc) then
                 write(iout,4)
              endif
c              write(iout,*) 'passing in trial vectors'
              call davdag('new trials',status,eigbuf,vec,hvec,
     1                     scr,n,iter,added,junk,junk,rjunk,root,
     2                     rjunk,rjunk)
c              write(iout,*) status
              if(status.eq.'none left') then
                 write(iout,5)
                 call davdag('finish',status,rjunk,rjunk,rjunk,rjunk,
     1                       junk,junk,junk,junk,junk,rjunk,rjunk,
     2                       rjunk,rjunk)
                 return
              elseif(status.eq.'done') then
c                     write(iout,*) 'multiplying by hamiltonian'
                     call honv(vec(1,nvc+1),hvec(1,nvc+1),hambuf,n,
     1                         added)         
                     nvc=nvc+added              
c                     write(iout,*) 'check matrix size'         
                     call davdag('check matrix size',status,eigbuf,vec,
     1                            hvec,scr,n,iter,added,junk,junk,rjunk,
     2                            rjunk,rjunk,rjunk)
              endif
           endif
      enddo
      call davdag('finish',status,rjunk,rjunk,rjunk,rjunk,junk,     
     1             junk,junk,junk,junk,rjunk,rjunk,rjunk,rjunk)
      write(iout,6)
      write(iout,7)
      do 40 i=1,nroots
         write(iout,8) i, eigbuf(i), root(i)
 40   continue
      call iosys('write real "davidson eigenvalues" to lamdat',nroots,
     1            root,0,' ')
      call iosys('write real "davidson eigenvectors" to lamdat',
     1            n*nroots,vec,0,' ')                       
      return
 1    format(/,20x,'iterative diagonalization procedure',/,1x,
     1             'number of zeroth order blocks = ',i3)
 2    format(/,1x,'diagonalizing block = ',i2,1x,'size = 'i4)
 3    format(/,1x,'block = ',i2,' not prediagonalized')
 4    format(/,1x,'warning: the next set of vectors added may exceed',
     1       /,1x,'         the maximum')
 5    format(/,1x,'no more trial vectors can be generated;quit')
 6    format(//,25x,'eigenvalue summary information')     
 7    format(/,5x,'     root     ',5x,' guess energy ',7x,
     1            ' final energy ')
 8    format(9x,i3,10x,e15.8,6x,e15.8)     
      end       
