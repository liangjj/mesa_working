c \documentclass{article}
c \usepackage{graphicx}
c \setkeys{Gin}{width=\linewidth}
c \title{DVDDAT: Enter Data for Davidson DVR Code}
c \author{Barry I. Schneider}
c \date{}
c \def \<{\langle}
c \def \>{\rangle}
c \begin{document}
c \maketitle

*deck dvddat.f
c***begin prologue     dvddat
c***date written       010829   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c***keywords           
c***author             schneider, barry (nsf)
c***source             
c***purpose            data entry for dvr davidson routine.
c***                   
c***references         
c
c***routines called    
c***end prologue       dvddat
      subroutine dvddat(card,cpass,cnverg,thresh,precon,n,nroot,
     1                  nattim,ntrial,maxit,maxvec,nblck,prdvd)
      implicit integer (a-z)
      real*8 cnverg, thresh, fpkey
      character*(*) card, cpass, precon
      character*80 chrkey
      logical dollar, logkey, prdvd
      dimension prdvd(12)
      common/io/inp, iout
      if( dollar('$dvrdvd',card,cpass,inp) ) then 
         cnverg=fpkey(card,'convergence',1.d-08,' ')
         thresh=fpkey(card,'overlap-tolerance',1.d-08,' ')
         precon=chrkey(card,'preconditioner','none',' ')
         if(precon.ne.'none') then
            nblck=intkey(card,'maximum-size-of-preconditioning-block',
     1                         200,' ')                      
         endif
         nroot=intkey(card,'number-of-roots',n,' ')
         nroot=min(nroot,n)
         maxvec=intkey(card,'maximum-number-of-vectors',2*nroot,' ')
         maxvec=min(maxvec,n)
         maxit=intkey(card,'maximum-number-of-iterations',maxvec,' ')
         nattim=intkey(card,'number-of-roots-at-a-time',nroot,' ')
         ntrial=intkey(card,'number-of-trial-vectors',nroot,' ')
         ntrial = max(ntrial,nroot)
	 ntrial=min(ntrial,maxvec)
         prdvd(1)=logkey(card,'print=davidson=trials',.false.,' ')
         prdvd(2)=logkey(card,'print=davidson=orthogonalized-trials',
     1                   .false.,' ')
         prdvd(3)=logkey(card,'print=davidson=h-on-vectors',.false.,' ')
         prdvd(4)=logkey(card,'print=davidson=small-matrix',.false.,' ')
         prdvd(5)=logkey(card,'print=davidson=eigenvalues',.false.,' ')
         prdvd(6)=logkey(card,'print=davidson=iteration-information',
     1                   .false.,' ')
         prdvd(7)=logkey(card,'print=davidson=transformed-vectors',
     1                   .false.,' ')
         prdvd(8)=logkey(card,'print=davidson=transformed-h-on-vectors',
     1                   .false.,' ')
         prdvd(9)=logkey(card,'print=davidson=residuals',.false.,' ')
         prdvd(10)=logkey(card,'print=davidson=new-raw-vectors',
     1                    .false.,' ')
         prdvd(11)=logkey(card,'print=davidson=overlaps',.false.,' ')
         prdvd(12)=logkey(card,'print=davidson=all',.false.,' ')	 
         if(prdvd(12)) then
            do 10 i=1,11
               prdvd(i)=.true.
 10         continue
         endif
      endif
      write(iout,1) nroot, nattim, ntrial, 
     1              thresh, cnverg, maxit, maxvec, precon, nblck
      return
 1    format(/,15x,'iterative diagonalization information',/,/,5x,
     1             'number of roots                    = ',i3,/,5x,
     2             'number of roots at a time          = ',i3,/,5x
     3             'number of trials                   = ',i3,/,5x
     4             'overlap tolerance                  = ',e15.8,/,5x,
     5             'convergence criterion              = ',e15.8,/,5x,
     6             'maximum number of iterations       = ',i6,/,5x,
     7             'maximum number of vectors          = ',i6,/,5x,           
     8             'preconditioning                    = ',a24,/,5x,
     9             'block size                         = ',i5)           
      end




