*deck conham.f
c***begin prologue     conham
c***date written       951229   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c***keywords           
c***author             schneider, barry (nsf)
c***source             
c***purpose            matrix elements of DVR/FEM basis.
c***                   
c***references         
c
c***routines called    
c***end prologue       conham
      subroutine conham(pmat,hmat,v,norm,n,npt,nreg,prn)
      implicit integer (a-z)
#ifdef DECPOINTER
      integer*8 pmat
#endif DECPOINTER
#ifdef SGIPOINTER
      integer*4 pmat
#endif SGIPOINTER
      real*8 mati, matj, hmat, v, norm
      real*8 dummy
      character*80 title
      logical prn, bridge
      dimension pmat(nreg), npt(nreg)
      dimension hmat(n,n), v(n), norm(n), prn(*)
      common/io/inp, iout
      pointer (pmati,mati(1))
      pointer (pmatj,matj(1))
      write(iout,1) n
      call rzero(hmat,n*n)
c
c     calculate the needed matrix elements
c
      row=1
      do 100 reg=1,nreg
c
c        locate the matrix elements for this region and the previous region
c
         pmati=pmat(reg)

c        locate the first word of the overlap, block and kinetic
c        energy matrices in this region.

         ovi=1
         bli=ovi+npt(reg)*npt(reg)
         kei=bli+npt(reg)*npt(reg)
         bridge=.true.
         start=2
         end=npt(reg)

c        since all regions except the first and last have a bridge
c        function, offset the row and column indices so they start
c        at the second function.  the first function is included as 
c        the bridge function from the previous interval.

         keii=kei + end + 1
         keij=kei + 1
         if(reg.ne.nreg) then
            pmatj=pmat(reg+1)
            ovj=1
            blj=ovj+npt(reg+1)*npt(reg+1)
            kej=blj+npt(reg+1)*npt(reg+1)
            kejj=kej
         endif
         if(reg.eq.1) then
            start=1
            keii=kei
         endif 
         if(reg.eq.nreg) then
            bridge=.false.
         endif 
         nfun=end-start+1
c
c        get the global starting and ending value
c
         last=row+nfun-1
         write(iout,2) reg, nfun, start, end, row, 
     1                                 last, bridge
c
c
c        if we consider only the lower triangle of the full matrix
c        and recognize that use of the bloch operator ensures
c        hermiticity, then a function in region i connects with
c        functions in the i - 1 region and the i region.  in the
c        i - 1 region it can only connect to a bridge function, if
c        it exists which has a piece in the i region.  
c        if the function in region i is a bridge function, it connects
c        to all lower functions in that region and with the first function
c        in the i + 1 region.
c

         if(reg.gt.1) then
c
c           do integrals between functions in region i and i - 1
c
c
c           ok we have a bridge function in region i - 1
c           it is a linear combination of F(n|i-1) + F(1|i).
c           we only need integrals from region i, suitably
c           weighted by the norms.
c            
c
            call makod(hmat(row,row-1),mati(keij),norm(row-1),
     #                 npt(reg),nfun,n,prn(1))
         endif
c
c        do integrals between functions all in region i 
c        we already have the normalizations
c
c        if we have a bridge function in region i
c        it is a linear combination of F(n|i) + F(1|i+1).
c        we only need integrals from region i and i + 1, suitably
c        weighted by the norms.
c
         call makd(hmat(row,row),mati(keii),matj(kejj),norm(row),
     #             bridge,npt(reg),nfun,n,prn(1))
         row=row+nfun
 100  continue
      do 200 i=1,n
         do 300 j=1,i
            hmat(j,i)=hmat(i,j)
 300     continue
 200  continue   
c
c     add the potential
c
      do 400 i=1,n
         hmat(i,i) = hmat(i,i) + v(i)
 400  continue   
      if(prn(2)) then
         title='full hamiltonian'
         call prntrm(title,hmat,n,n,n,n,iout)   
      endif
      return
 1    format(/,5x,'size of global basis set = ',i4)
 2    format(/,5x,'region = ',i3,1x,'number of functions = ',i3,/,5x,
     1            'starting function = ',i3,1x,'ending function = ',i3,
     2       /,5x,'global starting function = ',i3,1x,
     3            'global ending function = ',i3,
     4       /,5x,'bridge function = ',l1)
 4    format(10x,i4,12x,i4,1x,a38)
      end       



