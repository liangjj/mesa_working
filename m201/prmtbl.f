*deck @(#)prmtbl.f	5.1  11/6/94
      subroutine prmtbl(hding,vname,x,type,dx,nvar,lbl,nz,toang)
c***begin prologue     prmtbl.f
c***date written       yymmdd  
c***revision date      11/6/94      
c
c***keywords           
c***author             binkley, et al. (g82)
c***source             @(#)prmtbl.f	5.1   11/6/94
c***purpose            prints optimization information
c***description
c     
c    
c
c***references
c
c***routines called
c
c***end prologue       prmtbl.f
      implicit none
c     --- input variables -----
      integer hding,nvar,nz
      real*8 toang
c     --- input arrays (unmodified) ---
      integer type(nvar),lbl(nz)
      character*(*) vname(nvar)
      real*8 x(nvar),dx(nvar*(nvar+1)/2)
c     --- input arrays (scratch) ---
c     --- output arrays ---
c     --- output variables ---
c     --- scratch arrays ---
c     --- local variables ---
      integer inp,iout
      integer i,j,ii,nrep,cursor,decpl,cskipb
      character tname*16
      character*80 line
      logical degree,angs
      real*8 todeg,zero,one,f45,ten,hundrd,conver,value
c
      parameter (zero=0.0d+00,one=1.0d+00,f45=4.5d+01,ten=1.0d+01)
      parameter (hundrd=1.0d+02)
c
      common/io/inp,iout
c
 1010 format(5x,'initial parameters (angstroms and degrees)')
 1020 format(5x,'optimized parameters (angstroms and degrees)')
 1030 format(5x,'non-optimized parameters (angstroms and degrees)')
 1040 format(10x,4x,'name',7x,'value',5x,'derivative information'
     $                                   ,' (atomic units)')
 1050 format(1x,a80)
c
c     --- print heading.
      if(hding.eq.0) write(iout,1010)
      if(hding.eq.1) write(iout,1020)
      if(hding.eq.2) write(iout,1030)
      write(iout,1040)
c
c     --- loop over variables.
      todeg=f45/atan(one)
      do 10 i=1,nvar
c        --- append variable name.
         call crjust(vname(i),tname)
         line(1:)=' '//tname
c        --- convert from bohr/radian to angstroms/degree.
         conver=todeg
         degree=.true.
         angs=.false.
         if(nrep(i,lbl,nz).ne.0) then
            conver=toang
            angs=.true.
            degree=.false.
         endif
         value=x(i)*conver
c        --- tab, and insert value.
         cursor=20
         if(value.ge.zero) cursor=cursor+1
         if(value.lt.ten) cursor=cursor+1
         if(value.lt.hundrd) cursor=cursor+1
c
         if(angs) decpl=6
         if(degree) decpl=4
         call putfp(value,decpl,line,cursor)
c        --- append the derivative information.
         cursor=35
         j=abs(type(i))
         if(type(i).eq.-1) j=4
         if(j.eq.0) line(cursor:)='estimate d2e/dx2'
         if(j.eq.1) line(cursor:)='d2e/dx2='
         if(j.eq.2) line(cursor:)='update d2e/dxy'
         if(j.eq.3) line(cursor:)='calc d2e/dxdy, stepsize='
         if(j.eq.4) line(cursor:)='calc d2e/dx2 analytically'
         if(j.eq.5) line(cursor:)='d2e/dx2=identity'
         if(j.eq.97) line(cursor:)='d2e/dx2='
         if(j.eq.98) line(cursor:)='de/dx='
         if(j.eq.99) line(cursor:)='-de/dx='
c
c        --- find the end of the string.
         cursor=cskipb(line,' ')+1
c        --- append values if appropriate.
         if(j.ne.0.and.j.ne.2.and.j.ne.4.and.j.ne.5) then
            cursor=cursor+1
            ii=i*(i+1)/2
            if(j.eq.97.or.j.eq.98.or.j.eq.99) ii=i
            if(dx(ii).ge.zero) cursor=cursor+1
            if(abs(dx(ii)).lt.hundrd) cursor=cursor+1
            if(abs(dx(ii)).lt.ten) cursor=cursor+1
            call putfp(dx(ii),6,line,cursor)
         endif
c
c        --- write the line.
         write(iout,1050) line
   10 continue
c
c
      return
      end
