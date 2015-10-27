*deck drvh1.f
c***begin prologue     drvh1
c***date written       000710   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c***keywords           one-dim
c***author             schneider, barry (nsf)
c***source             3-dim
c***purpose            driver for one body hamiltonian matrix elements 
c***                   and indices for one dimensional dvr hamiltonian. 
c***                   
c***                   
c***references         
c
c***routines called    
c***end prologue       drvh1
      subroutine drvh1(pham,pv,phone,phamil,n,len,nonzro,ngot,pack,
     1                 drctv,prn)
      implicit integer (a-z)
      real*8 ham, v, one, hamil
      logical pack, prn
      character*80 title
      character*(*) drctv
      pointer (pham,ham(1))
      pointer (pv,v(1))
      pointer (phone,one(1))
      pointer (phone,ione(1))
      pointer (phamil,hamil(1))
      common/io/inp, iout
      if(pack) then
         len=n*n
         hbuf=1
         buf=wpadti(hbuf+len)
         diag=iadtwp(buf+2*len)
         need=wpadti(diag+n) 
         call getmem(need,phone,ngot,'one',0)
         call h1pac(ham,v,one(hbuf),ione(buf),one(diag),n,len,
     1              nonzro,prn)
         write(iout,1) nonzro
         if(prn) then
            call prpac1(one(hbuf),ione(buf),one(diag),n,len,nonzro)
         endif
      endif
      if(drctv.eq.'diagonalize') then
         htmp=1 
         eig=htmp+n*n
         work=eig+n
         need=wpadti(work+5*n)
         call getmem(need,phamil,junk,'hamil',0)
c
c        diagonalize the free-particle hamiltonian
c
         call copy(ham,hamil(htmp),n*n)
         call vec2di(ham,ham,v,'subtract',n)
         call dsyev('v','l',n,ham,n,hamil(eig),hamil(work),5*n,info)
         title='free particle eigenvalues'
         call prntrm(title,hamil(eig),n,1,n,1,iout)
         call copy(hamil(htmp),ham,n*n)
         call dsyev('v','l',n,hamil(htmp),n,hamil(eig),hamil(work),
     1              5*n,info)
         title='one particle eigenvalues'
         call prntrm(title,hamil(eig),n,1,n,1,iout)
      endif
      return
 1    format(/,1x,'number of non-zero matrix elements = ',i5)
      end       



