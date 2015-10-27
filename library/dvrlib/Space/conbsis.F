c \documentclass{article}
c \usepackage{graphicx}
c \setkeys{Gin}{width=\linewidth}
c \title{Construct Final DVR Basis}
c \author{Barry I. Schneider}
c \date{}
c \def \<{\langle}
c \def \>{\rangle}
c \begin{document}
c \maketitle 

*deck conbsis.f
c***begin prologue     conbsis
c***date written       951229   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c***keywords           
c***author             schneider, barry (nsf)
c***source             
c***purpose            DVR/FEM basis.
c***                   
c***references         
c
c***routines called    
c***end prologue       conbsis
      subroutine conbsis(pgrid,x,xwt,f,df,ddf,norm,edge,n,npt,nreg,
     1                   prn)
      implicit integer (a-z)
#ifdef DECPOINTER
      integer*8 pgrid
#endif
#ifdef SGIPOINTER
      integer*4 pgrid
#endif
      real*8 grdi, grdj, x, xwt, f, df, ddf, flst, edge, norm
      character*80 title
      logical prn, bridge
      dimension pgrid(nreg), x(n), xwt(n), f(n,n), df(n,n), ddf(n,n)
      dimension norm(n), prn(*)
      dimension npt(nreg)
      dimension edge(nreg+1)
      common/io/inp, iout
      pointer (pgi,grdi(1))
      pointer (pgj,grdj(1))
      call rzero(f,n*n)
      call rzero(df,n*n)
      call rzero(ddf,n*n)

c
c     calculate the needed functions
c
      
      row=1
      do 10 reg=1,nreg
         
c
c        get the starting and ending function in this region and the 
c        pointers.
c

         pgi=pgrid(reg)
         qi=1
         wti=qi+npt(reg)
         pi=wti+npt(reg)
         dpi=pi+npt(reg)*npt(reg)
         ddpi=dpi+npt(reg)*npt(reg)
         if(reg.ne.nreg) then
            pgj=pgrid(reg+1)
            qj=1
            wtj=qj+npt(reg+1)
            pj=wtj+npt(reg+1)
            dpj=pj+npt(reg+1)*npt(reg+1)
            ddpj=dpj+npt(reg+1)*npt(reg+1)
         endif
         bridge=.true.
         start=2
         end=npt(reg)
         if(reg.eq.1) then
            start=1
         endif
         if(reg.eq.nreg) then
            bridge=.false.
         endif 
         nfun = end - start + 1
         last = row + nfun - 1
         write(iout,2) reg, nfun, start, end, row, last, bridge

c
c        calculate the functions and derivatives in this region
c        region i, including the bridge function if present.
c        if there is a bridge function, it contains pieces
c        from region i and region i + 1.
c

         call filfun(x(row),xwt(row),f(row,row),df(row,row),
     1               ddf(row,row),norm(row),grdi(qi),grdj(qj),
     2               grdi(wti),grdj(wtj),grdi(pi),grdj(pj),
     4               grdi(dpi),grdj(dpj),grdi(ddpi),grdj(ddpj),
     6               npt(reg),npt(reg+1),edge(reg),edge(reg+1),
     7               edge(reg+2),bridge,nfun,n,start)
         row = row + nfun
 10   continue   
      if(prn(1)) then
         title='global coordinates'
         call prntrm(title,x,n,1,n,1,iout)   
         title='global weights'
         call prntrm(title,xwt,n,1,n,1,iout)   
         title='global norm'
         call prntrm(title,norm,n,1,n,1,iout)   
      endif
      if(prn(2)) then
         title='global basis'
         call prntrm(title,f,n,n,n,n,iout)   
         title='1 derivative of global basis'
         call prntrm(title,df,n,n,n,n,iout)   
         title='2 derivative of global basis'
         call prntrm(title,ddf,n,n,n,n,iout)   
      endif
      return
 1    format(/,5x,'size of global basis set = ',i5)
 2    format(/,5x,'region = ',i3,1x,'number of functions = ',i5,/,5x,
     1            'starting function = ',i5,1x,'ending function = ',i5,
     2       /,5x,'global starting function = ',i5,1x,
     3            'global ending function = ',i5,
     4       /,5x,'bridge function = ',l1)
 4    format(10x,i4,12x,i4,1x,a38)
      end       




