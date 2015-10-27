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
     1                 scr,work,ipvt,thresh,cnverg,ops,n,nroots,nattim,
     2                 mxvc,iter,prgues,prvec,prtflg,prres)      
      implicit integer (a-z)
      real*8 hambuf, eigbuf
      real*8 vec, hvec, root, dvdvec, dvdmat, cvn, scr, work
      real*8 cnverg, thresh, fpkey, rjunk
      logical logkey, schmdt, prgues, prvec, restrt, prres
      character*(*) ops
      character*6 file
      character*80 title
      character*20 status, resslv, chrkey
      character*8 prtflg
      dimension eigbuf(n), hambuf(n,n), dvdmat(*), dvdvec(*), cvn(n)
      dimension vec(n,mxvc), hvec(n,mxvc), scr(n,nroots), root(nroots)
      dimension work(*), ipvt(n)
      common/io/inp, iout
      call rzero(scr,n*nroots)
      if (prtflg.eq.'minimum') then
         status='noprint'
      else
         status='print'
      end if 
      file='lamdat'
      restrt=logkey(ops,'davidson=restart',.false.,' ')
      resslv=chrkey(ops,'davidson=solver','diagonal',' ')
      gsiter=1
      if (resslv.eq.'gauss-seidel') then
          gsiter=intkey(ops,'davidson=gauss-seidel-iterations',10,' ')
      endif
      schmdt=logkey(ops,'davidson=two-schmidt-orthogonalizations',
     1              .false.,' ')
      call rzero(vec,n*mxvc)
      write(iout,1) nroots, nattim, mxvc, iter, thresh, cnverg,
     1              resslv, gsiter
      if(.not.restrt) then
          gsize=intkey(ops,'davidson=size-of-guess-matrix',nroots,' ')
          cntg=intkey(ops,'davidson=number-of-guess-vectors',nroots,' ')     
          cntg=min(cntg,gsize)
          write(iout,2) gsize, cntg
          call rzero(work,nroots)
          call guess(hambuf,eigbuf,vec,work,hvec,n,gsize,nroots,prgues)
          call iosys('write real "guess eigenvalues" to lamdat',nroots,
     1                work,0,' ')
      else
          call iosys('read integer "number of davidson iterates" '//
     1               'from lamdat',1,gsize,0,' ')
          cntg=intkey(ops,'davidson=number-of-guess-vectors',
     1                gsize,' ')
          write(iout,3) gsize, cntg
          call oldin(vec,work,n,gsize,nroots,prgues)
c          cntg=gsize
      endif
c      call trnham(hambuf,work,vec,hvec,n,gsize)
c
      call davdag('initialize',status,file,rjunk,thresh,rjunk,
     1             rjunk,n,iter,nroots,iout,nattim,rjunk,0.d0,
     2             cnverg,rjunk,rjunk,schmdt,prres)
      newtr=cntg   
      if (newtr.le.iter) then
          call honv(vec,hvec,hambuf,n,newtr)
      else
          call lnkerr('problem with first step in davdag')
      endif
c***********************************************************************
c
c***********************************************************************
c                 b. solve the small eigenvalue equation and test
c                    for convergence. for unconverged roots add
c                    new vectors using perturbation theory. the number
c                    of these vectors is returned as nresid and they
c                    reside in the scr array upon return.
        nvc=newtr
        do while (status.ne.'all roots converged')                    
           call davdag('check matrix size',status,file,rjunk,rjunk,
     1                  rjunk,rjunk,n,iter,newtr,junk,junk,junk,rjunk,
     2                  rjunk,rjunk,rjunk,schmdt,prres)
           if(status.ne.'ok') then
              write(iout,4)
              nvc=nvc-newtr
              call davdag('cleanup',status,file,rjunk,vec,hvec,
     1                     rjunk,n,iter,nvc,junk,junk,junk,rjunk,root,
     2                     rjunk,rjunk,schmdt,prres)
              call davdag('finish',status,file,rjunk,rjunk,rjunk,
     1                     rjunk,junk,junk,junk,junk,junk,junk,rjunk,
     2                     rjunk,rjunk,rjunk,schmdt,prres)
              return
           endif
           call davdag('solve',status,file,eigbuf,vec,hvec,scr,n,iter,
     1                  ntoadd,numcnv,nvc,junk,dvdmat,root,dvdvec,
     2                  cvn,schmdt,prres)
           if(status.eq.'all roots converged') then
              call davdag('cleanup',status,file,rjunk,vec,hvec,
     1                     rjunk,n,iter,nvc,junk,junk,junk,rjunk,root,
     2                     rjunk,rjunk,schmdt,prres)
              call davdag('finish',status,file,rjunk,rjunk,rjunk,
     1                     rjunk,junk,junk,junk,junk,junk,junk,rjunk,
     2                     rjunk,rjunk,rjunk,schmdt,prres)
           elseif(status.eq.'continue') then
                  status=resslv 
                  call davdag('new trials',status,file,eigbuf,vec,
     1                         hambuf,scr,n,iter,newtr,numcnv,gsiter,
     2                         ipvt,dvdmat,root,dvdvec,work,schmdt,
     3                         prres)
                  if(status.eq.'none left') then
                     write(iout,5) 
                     call davdag('cleanup',status,file,rjunk,vec,
     1                            hvec,rjunk,n,iter,nvc,junk,junk,
     2                            junk,rjunk,root,rjunk,rjunk,
     3                            schmdt,prres)
                     call davdag('finish',status,file,rjunk,rjunk,
     1                            rjunk,rjunk,junk,junk,junk,junk,
     2                            junk,junk,rjunk,rjunk,rjunk,rjunk,
     3                            schmdt,prres)
                     return                        
                  elseif(status.eq.'done') then
                     call honv(vec(1,nvc+1),hvec(1,nvc+1),hambuf,n,
     1                         newtr)         
                     nvc=nvc+newtr
                  endif
           endif
      enddo
      call iosys('read real "guess eigenvalues" from lamdat',nroots,
     1            work,0,' ')
      write(iout,6)
      write(iout,7)
      do 60 i=1,nroots
         write(iout,8) i, work(i), root(i)
 60   continue
      if(prvec) then
         title='final eigenvectors'
         call prntrm(title,vec,n,nroots,n,n,iout)
      endif
      call copy(vec,hambuf,n*nroots)
      call copy(root,eigbuf,nroots)          
      write(iout,9) nvc
      return
 1    format(/,5x,'davidson iterative diagonalization',/,5x,
     1            'number of roots                          = ',i3,/,5x,
     2            'number calculated at a time              = ',i3,/,5x,
     3            'maximum number of vectors stored in core = ',i3,/,5x,
     4            'maximum number of iterations             = ',i3,/,5x,
     5            'threshold criterion for vectors          = ',e15.8,
     6                                                           /,5x,
     7            'convergence criterion                    = ',e15.8,
     8                                                          /,5x,
     9            'preconditioner/smoother                  = ',a20,
     x                                                        /5x,
     x            'number of smoothing iterations           = ',i3)
 2    format(/,5x,'initial vectors from guess',/,5x,
     1             'size of guess matrix    = ',i3,/,5x,
     2             'number of guess vectors = ',i3)
 3    format(/,5x,'initial vectors from davidson restart file',/,5x,
     1             'size of guess matrix    = ',i3,/,5x,
     2             'number of guess vectors = ',i3)
 4    format(/,5x,'size of vector space exceeded. cleanup and quit')
 5    format(/,1x,'no more trial vectors can be generated;quit')
 6    format(//,25x,'eigenvalue summary information')     
 7    format(/,5x,'     root     ',5x,' guess energy ',7x,
     1            ' final energy ')
 8    format(9x,i3,10x,e15.8,6x,e15.8)     
 9    format(/,1x,'final number of davidson iterates = ',i5)
      end       
