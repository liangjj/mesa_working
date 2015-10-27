c \documentclass{article}
c \usepackage{graphicx}
c \setkeys{Gin}{width=\linewidth}
c \title{LINDAT: Enter Data For GMRES COde}
c \author{Barry I. Schneider}
c \date{}
c \def \<{\langle}
c \def \>{\rangle}
c \begin{document}
c \maketitle

*deck lindat.f
c***begin prologue     lindat
c***date written       010829   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c***keywords           
c***author             schneider, barry (nsf)
c***source             
c***purpose            data entry for iterative solver routine.
c***                   
c***references         
c
c***routines called    
c***end prologue       lindat
      subroutine lindat(card,cpass,cnverg,thresh,precon,n,
     1                  maxit,maxvec,nblck,prlin)
      implicit integer (a-z)
      real*8 cnverg, thresh, fpkey
      character*(*) card, cpass, precon
      character*80 chrkey
      logical dollar, logkey, prlin
      dimension prlin(*)
      common/io/inp, iout
      if( dollar('$gmres',card,cpass,inp) ) then 
         cnverg=fpkey(card,'convergence',1.d-08,' ')
         thresh=fpkey(card,'overlap-tolerance',1.d-08,' ')
         precon=chrkey(card,'preconditioner','none',' ')
         maxvec=intkey(card,'maximum-number-of-vectors',n,' ')
         maxit=intkey(card,'maximum-number-of-iterations',maxvec,' ')         
         nblck=0
         if(precon.ne.'none') then
            nblck=intkey(card,'maximum-size-of-preconditioning-block',
     1                         200,' ')                      
         endif
         prlin(1)=logkey(card,'print=iterative-solve=all',.false.,' ')
         prlin(2)=logkey(card,'print=iterative-solve=trials',
     1                  .false.,' ')
         prlin(3)=logkey(card,'print=iterative-solve=vectors',
     1                   .false.,' ')
         prlin(4)=logkey(card,'print=iterative-solve=h-on-vectors',
     1                   .false.,' ')
         prlin(5)=logkey(card,'print=iterative-solve=hamiltonian',
     1                   .false.,' ')
         prlin(6)=logkey(card,'print=iterative-solve=iteration-'//
     1                   'information',.false.,' ')
         prlin(7)=logkey(card,'print=iterative-solve=residuals',
     1                   .false.,' ')
         prlin(8)=logkey(card,'print=iterative-solve='//
     1                   'transformed-vectors',.false.,' ')
         prlin(9)=logkey(card,'print=iterative-solve='//
     1                   'transformed-h-on-vectors',
     2                   .false.,' ')               
         prlin(10)=logkey(card,'print=iterative-solve='//
     1                    'new-trial-vectors',.false.,' ')
         prlin(11)=logkey(card,'print=iterative-solve=overlaps',
     1                .false.,' ')
         if(prlin(1)) then
            do 10 i=2,11
               prlin(i)=.true.
 10         continue
         endif
      endif
      write(iout,1) thresh, cnverg, precon, maxit, maxvec, nblck
      return
 1    format(/,15x,'iterative linear system information',/,/,5x,
     1             'overlap tolerance                  = ',e15.8,/,5x,
     2             'convergence criterion              = ',e15.8,/,5x,
     3             'preconditioning                    = ',a24,/,5x,
     4             'maximum number of iterations       = ',i6,/,5x,
     5             'maximum number of vectors          = ',i6,/5x,           
     4             'block size                         = ',i5)           
      end


