*deck drvlin.f
c***begin prologue     drvlin
c***date written       951229   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c***keywords           linear systems, variation-iteration
c***author             schneider, barry (nsf)
c***source             
c***purpose            driver large linear systems solver.
c***references         
c
c***routines called    
c***end prologue       drvlin
      subroutine drvlin(hbuf,ibuf,diag,vec,hvec,energy,rhs,dinv,
     1                 smat,svec,cvn,resid,space,ipvt,thresh,cnverg,
     2                 ops,n,nrhs,nattim,mxvc,iter,lenbuf,ntot,incore)      
      implicit integer (a-z)
      real*8 hbuf, diag
      real*8 vec, hvec, svec, smat, rhs, cvn, resid, space
      real*8 cnverg, thresh, rjunk, energy, dinv
      logical logkey, schmdt, prvec, prres, incore
      character*(*) ops
      character*6 file
      character*80 title
      character*32 status, type, chrkey
      character*8 prtflg
      dimension hbuf(lenbuf), ibuf(2,lenbuf), diag(n), rhs(n,nrhs)
      dimension smat(mxvc,mxvc,2), svec(mxvc,nrhs,2), cvn(n), dinv(n)
      dimension vec(n,mxvc), hvec(n,mxvc), resid(n,mxvc)
      dimension space(*), ipvt(n)
      common/io/inp, iout
      write(iout,1) 
      call iosys ('read character "print flag" from rwf',-1,0,0,prtflg)      
      prvec=logkey(ops,'print=linalg=final-vectors',.false.,' ')
      prres=logkey(ops,'print=linalg=residuals',.false.,' ')
      schmdt=logkey(ops,'linalg=two-schmidt-orthogonalizations',
     1              .false.,' ')
      type=chrkey(ops,'linslv=form','one-minus-matrix',' ')
      if(type.eq.'matrix') then
         do 10 i=1,n
            dinv(i)=1.d0
 10      continue      
      elseif(type.eq.'one-plus-matrix') then
         do 20 i=1,n
            dinv(i)=1.d0/diag(i)
 20      continue
      elseif(type.eq.'one-minus-matrix') then
         do 30 i=1,n
            dinv(i) = 1.d0/( energy - diag(i) )
 30      continue
      else
         call lnkerr('error in form of equations')
      endif         
      call vmmul(dinv,rhs,rhs,n,nrhs)   
      ntrial=nrhs
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
      debug=1
      nsol=nrhs
      call linslv('initialize',type,vec,rhs,thresh,cnverg,rjunk,
     1             rjunk,junk,mxvc,n,iout,nsol,nvec,debug,schmdt)
      maxtrp=mxvc/nsol+1
      do 40 trips=1,maxtrp
         nold=nvec
         test=nold+nsol
         if(test.le.mxvc) then
            call linslv('generate vectors',status,vec,rjunk,rjunk,rjunk,
     1                   svec(1,1,2),rjunk,rjunk,ntrial,n,mxvc,mxvc,nvec,
     2                   debug,schmdt)
            if(status.eq.'iteration limit exceeded') then
               write(iout,2) status, nvec
               return 
            endif
            if(ntrial.ne.0) then
               iold=iold+1
               call honv(hbuf,ibuf,vec(1,iold),hvec(1,iold),n,ntrial,
     1                   lenbuf,ntot,incore,title)
               if(type.eq.'matrix') then
                  call gvmmul(diag,vec(1,iold),hvec(1,iold),
     1                        hvec(1,iold),n,ntrial)   
               else
                  call vmmul(dinv,hvec(1,iold),hvec(1,iold),n,ntrial)            
               endif
               call linslv ('solve equations',status,vec,
     1                       hvec,smat(1,1,2),smat(1,1,1),svec(1,1,2),
     2                       svec(1,1,1),rhs,ipvt,n,iter,iter,nvec,debug,
     3                       schmdt)
               if(status.eq.'starting to diverge') then                
                  write(iout,3) nold
                  return
               endif
            else
               write(iout,4) nvec
               return
            endif
         else      
            call linslv('final best solution',status,vec,rjunk,
     1                   rjunk,rjunk,svec(1,1,2),rjunk,rhs,junk,
     2                   n,mxvc,mxvc,nvec,debug,schmdt)
               write(iout,2) status, nvec
               return
        endif 
   40 continue                                         
      return
 1    format(/,5x,'enter iterative linear system solver')
 2    format(/,5x,'message from linslv:',/,26x,a32,/,5x,
     1            'number of vectors used= ',i5)
 3    format(/,5x,'solution statring to diverge.  will quit',/,5x,
     1            'number of vectors used = ',i5)     
 4    format(/,5x,'no more vectors.  we are done.',/,10x,
     1            'number of vectors used = ',i5)     
      end       
