*deck @(#)fptest.f	5.1  11/6/94
      subroutine fptest(nvar,archiv,convrg,abnrml,fpcycl,mxcycl,convf,
     $                  d1var,pool1,yold)
c***begin prologue     fptest.f
c***date written       yymmdd  
c***revision date      11/6/94      
c
c***keywords           
c***author             binkley, et al. (g82)
c***source             @(#)fptest.f	5.1   11/6/94
c***purpose            
c***description
c     
c    
c
c***references
c
c***routines called
c
c***end prologue       fptest.f
      implicit none
c     --- input variables -----
      integer nvar,fpcycl,mxcycl
      logical archiv,convrg,abnrml
      real*8 convf
c     --- input arrays (unmodified) ---
      real*8 d1var(nvar),pool1(nvar),yold(nvar)
c     --- input arrays (scratch) ---
c     --- output arrays ---
c     --- output variables ---
c     --- scratch arrays ---
c     --- local variables ---
      integer inp,iout
      real*8 frms,sdot
      character*3 result
c
      common/io/inp,iout
c
 1000 format(5x,'convergence tests:',
     $      /10x,16x,'item',8x,'value',4x,'threshold',2x,'converged?')
 1002 format(10x,7x,'maximum force',5x,f9.6,5x,f9.6,5x,a3)
 1004 format(10x,7x,'    rms force',5x,f9.6,5x,f9.6,5x,a3)
 1010 format(5x,'fletcher-powell optimization terminated.')
 1020 format(5x,'cycle number ',i3,' out of a maximum of ',i3,
     $      /5x,'all quantities printed in hartrees-bohr-radians.')
 1030 format(5x,'forces below ',e10.3,' after ',i3,' cycles.')
 1050 format(5x,'maximum number of steps exceeded...',
     $          ' forces not below threshold.')
 1060 format(5x,'archiving disabled.')
c
c     --- convergence checks.
      write(iout,1020) fpcycl,mxcycl
      write(iout,1000)
      archiv=.false.
      convrg=.false.
      frms=sqrt(sdot(nvar,d1var,1,d1var,1)/float(nvar))
      call convgd(frms,convf,result)
      write(iout,1004) frms,convf,result
      call iosys('write real "rms force" to rwf',1,frms,0,' ')
c
      if(frms.le.convf) then
         write(iout,1010)
         write(iout,1030) convf,fpcycl
         archiv=.true.
         convrg=.true.
      else if(fpcycl.ge.mxcycl) then
         write(iout,1010)
         write(iout,1050)
         write(iout,1060)
         abnrml=.true.
         archiv=.false.
         convrg=.false.
         call iosys('write integer archive to rwf',1,.false.,0,' ')
      endif
c
c
      return
      end
