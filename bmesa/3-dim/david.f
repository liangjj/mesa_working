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
      subroutine david(hbuf,ibuf,diag,eig,vec,hvec,dvdmat,dvdvec,
     1                 cvn,resid,eig0,ipvt,thresh,cnverg,ops,n,nroots,
     2                 nattim,mxvc,iter,lenbuf,ntot,incore)      
      implicit integer (a-z)
      real*8 hbuf, diag, eig
      real*8 vec, hvec, dvdvec, dvdmat, cvn, eig0, resid
      real*8 cnverg, thresh, rjunk
      logical logkey, schmdt, prvec, prres, incore
      character*(*) ops
      character*6 file
      character*80 title
      character*20 status
      character*8 prtflg
      dimension hbuf(lenbuf), ibuf(2,lenbuf), diag(n)
      dimension eig(n), dvdmat(*), dvdvec(*), cvn(n)
      dimension vec(n,mxvc), hvec(n,mxvc), resid(n,mxvc)
      dimension eig0(n), ipvt(n)
      common/io/inp, iout
      call iosys ('read character "print flag" from rwf',-1,0,0,prtflg)      
      prvec=logkey(ops,'print=davidson=final-vectors',.false.,' ')
      prres=logkey(ops,'print=davidson=residuals',.false.,' ')
      schmdt=logkey(ops,'davidson=two-schmidt-orthogonalizations',
     1              .false.,' ')
      north=1
      if(schmdt) then
         north=2
      endif   
      if (prtflg.eq.'minimum') then
         status='noprint'
      else
         status='print'
      end if 
      file='lamdat'
      write(iout,1) 
      call copy(eig,eig0,nroots)
      call davdag('initialize',status,file,rjunk,thresh,rjunk,
     1             rjunk,n,iter,nroots,iout,nattim,rjunk,0.d0,
     2             cnverg,rjunk,rjunk,schmdt,prres)
      newtr=nroots   
      if (newtr.le.iter) then
          call honv(hbuf,ibuf,vec,hvec,n,newtr,lenbuf,ntot,
     1              incore,title)
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
c                    reside in the resid array upon return.
        nvc=newtr
        do while (status.ne.'all roots converged')                    
           call davdag('check matrix size',status,file,rjunk,rjunk,
     1                  rjunk,rjunk,n,iter,newtr,junk,junk,junk,rjunk,
     2                  rjunk,rjunk,rjunk,schmdt,prres)
           if(status.ne.'ok') then
              write(iout,2)
              nvc=nvc-newtr
              call davdag('cleanup',status,file,rjunk,vec,hvec,
     1                     rjunk,n,iter,nvc,junk,junk,junk,rjunk,eig,
     2                     rjunk,rjunk,schmdt,prres)
              call davdag('finish',status,file,rjunk,rjunk,rjunk,
     1                     rjunk,junk,junk,junk,junk,junk,junk,rjunk,
     2                     rjunk,rjunk,rjunk,schmdt,prres)
              return
           endif
           call davdag('solve',status,file,diag,vec,hvec,resid,n,iter,
     1                  ntoadd,numcnv,nvc,junk,dvdmat,eig,dvdvec,
     2                  cvn,schmdt,prres)
           if(status.eq.'all roots converged') then
              call davdag('cleanup',status,file,rjunk,vec,hvec,
     1                     rjunk,n,iter,nvc,junk,junk,junk,rjunk,eig,
     2                     rjunk,rjunk,schmdt,prres)
              call davdag('finish',status,file,rjunk,rjunk,rjunk,
     1                     rjunk,junk,junk,junk,junk,junk,junk,rjunk,
     2                     rjunk,rjunk,rjunk,schmdt,prres)
           elseif(status.eq.'continue') then
                  call davdag('new trials',status,file,diag,vec,
     1                         rjunk,resid,n,iter,newtr,numcnv,junk,
     2                         ipvt,dvdmat,eig,dvdvec,rjunk,schmdt,
     3                         prres)
                  if(status.eq.'none left') then
                     write(iout,3) 
                     call davdag('cleanup',status,file,rjunk,vec,
     1                            hvec,rjunk,n,iter,nvc,junk,junk,
     2                            junk,rjunk,eig,rjunk,rjunk,
     3                            schmdt,prres)
                     call davdag('finish',status,file,rjunk,rjunk,
     1                            rjunk,rjunk,junk,junk,junk,junk,
     2                            junk,junk,rjunk,rjunk,rjunk,rjunk,
     3                            schmdt,prres)
                     return                        
                  elseif(status.eq.'done') then
                     call honv(hbuf,ibuf,vec(1,nvc+1),hvec(1,nvc+1),n,
     1                         newtr,lenbuf,ntot,incore,title)
                     nvc=nvc+newtr
                  endif
           endif
      enddo
      write(iout,4)
      write(iout,5)
      do 60 i=1,nroots
         write(iout,6) i, eig0(i), eig(i)
 60   continue
      if(prvec) then
         title='final eigenvectors'
         call prntrm(title,vec,n,nroots,n,nroots,iout)
      endif
      write(iout,7) nvc
      return
 1    format(/,5x,'enter davidson diagonalization')
 2    format(/,5x,'size of vector space exceeded. cleanup and quit')
 3    format(/,1x,'no more trial vectors can be generated;quit')
 4    format(//,25x,'eigenvalue summary information')     
 5    format(/,5x,'     root     ',5x,' guess energy ',7x,
     1            ' final energy ')
 6    format(9x,i3,10x,e15.8,6x,e15.8)     
 7    format(/,1x,'final number of davidson iterates = ',i5)
      end       
