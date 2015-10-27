c \documentclass{article}
c \usepackage{graphicx}
c \setkeys{Gin}{width=\linewidth}
c \title{DRVSPLT: Driver for Split Operator Solution of Schroedinger
c        Equation in DVR Representation}
c \author{Barry I. Schneider}
c \date{}
c \def \<{\langle}
c \def \>{\rangle}
c \begin{document}
c \maketitle

*deck drvsplt.f 
c***begin prologue     drvsplt
c***date written       000702   (yymmdd)
c***revision date               (yymmdd)
c***keywords           split-operator,  dvr
c***                   
c***author             schneider, b. i.(nsf)
c***source             drvsplt
c***purpose            driver for solution of time independent
c***                   schroedinger equation using split operator.
c
c***references       
c
c***routines called    iosys, util and mdutil
c***end prologue       drvsplt
      program drvsplt
c
      implicit integer (a-z)
      parameter ( maxreg=100 ) 
      character*4096 ops
      character*2 itoc
      character*80 cpass, chrkey, title
      character*24 coord, sym, typpot, key
      character*800 card
      logical logkey, prn, pack1
      real*8 edge, fpkey 
      integer*8 pham, px, py, pz
      real*8 hx, hy, hz
      dimension n(4), nphy(4), coord(4), key(4), typpot(4)
      dimension edge(maxreg+1)
      dimension prn(30)
      dimension pham(4), ngrds(4)
      dimension nwham(4)
      common/io/inp, iout      
      common/punch/pun
      pointer(px,hx(1)), (py,hy(1)), (pz,hz(1))
c
      call drum
      write(iout,*)
      call iosys ('read character options from rwf',-1,0,0,ops)
      sym=chrkey(ops,'symmetry','unsymmetric',' ')
      dim=intkey(ops,'number-of-space-dimensions',1,' ')
      pack1=logkey(ops,'pack',.false.,' ')
      if(sym.eq.'symmetric') then
         pack1=.false.
      endif
      write(iout,1) dim
c
c     read some basic grid information from file.  this is used for all
c     dimensions and regions. 
c
c
c     read in information for each dimension
c
      nphy(dim+1)=1
      do 10 i=1,dim
         coord(i)=chrkey(ops,'dimension-'//itoc(i),'x',' ')
         call basis(ops,cpass,card,i,pham(i),edge,ngrds(i),coord(i),
     1              key(i),typpot(i),n(i),nphy(i),nwham(i),
     2              .false.,prn)
         nphy(dim+1)=nphy(dim+1)*nphy(i) 
 10   continue
      px=pham(1) 
      py=pham(2)
      pz=pham(3)
      h1=1
      v1=h1+nphy(1)*nphy(1)
      u1=v1+nphy(1)
      s1=u1+nphy(1)*nphy(1)
      eig1=s1+2*nphy(1)
      srf1=eig1+nphy(1)
      add=srf1+2+2*nphy(1)*nphy(1)+3*nphy(1)
      q1=add
      wt1=q1+nphy(1)+1
      h2=1
      v2=h2+nphy(2)*nphy(2)
      u2=v2+nphy(2)
      s2=u2+nphy(2)*nphy(2)
      eig2=s2+2*nphy(2)
      srf2=eig2+nphy(2)
      add=srf2+2+2*nphy(2)*nphy(2)+3*nphy(2)
      q2=add
      wt2=q2+nphy(2)+1
      call chainx(0)               
      stop
 1    format(/,20x,'split operator dvr basis function code',//,20x,
     1             'number of spatial dimensions = ',i1)      
      end

