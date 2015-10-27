*deck poles
      subroutine poles(eig,hfbox,enrm,v0,rbox,convg,l,niter,ncon,exact)
      implicit integer (a-z)
      real*8 eig, hfbox, enrm, eguess, v0, rbox, convg, ecorr, dele
      real*8 rmte, drmte, rmta, drmta, beta, dbeta
      real*8 funct, dfunct, fac, dfac
      dimension eig(ncon), hfbox(ncon), enrm(ncon,2)
      logical exact
      common/io/inp, iout
      if (v0.gt.0.d0) then
          write(iout,*) 'must be a well not a barrier'
          return
      endif
      count=0
      do 10 ipol=1,ncon
         eguess=eig(ipol)+1.d-03
         if (eguess.lt.0.d0) then
             count=count+1 
             write(iout,1)
             do 20 i=1,niter
                beta=-2.d0*eguess
                beta=sqrt(beta)
                dbeta=-1.d0/beta
                fac=1.d0/beta
                dfac=-dbeta/(beta*beta)
                if (l.ne.0) then
                    fac=beta + 1.d0/(rbox*(1.d0+beta*rbox))
                    fac=1.d0/fac
                    dfac=-dbeta*fac*fac*(1.d0 -1.d0/((1.d0+beta*rbox)*
     1                                               (1.d0+beta*rbox)))
                endif
                call rmat(rmta,drmta,rmte,drmte,hfbox,eig,rbox,
     1                    eguess,v0,l,ncon,1,.false.)
                if (exact) then
                    funct=rmte + fac
                    dfunct=drmte + dfac 
                else
                    funct=rmta + fac
                    dfunct=drmta + dfac
                endif                     
                ecorr=eguess-funct/dfunct
                dele=abs(eguess-ecorr)
                write(iout,2) i, ecorr, dele
                if (dele.le.convg) then
                    enrm(count,1)=ecorr
                    enrm(count,2)=(beta*beta*drmate+.5d0/beta)*
     1                             exp(-2.d0*beta*rbox)
                    enrm(count,2)=1.d0/sqrt(enrm(count,2))
                    write(iout,3) i, ecorr, funct, enrm(count,2)
                    go to 10
                endif
                eguess=ecorr
   20        continue
             write(iout,4) i           
         endif   
   10 continue
    1 format(/,5x,'iteration',2x,'     root      ',5x,
     1            '   delta e   ',/)
    2 format(7x,i3,6x,e15.8,2x,e15.8)          
    3 format(/,5x,'eigenvalue converged after ',i3,' iterations',/5x,
     1            'value of converged eigenvalue = ',e15.8,/,5x,
     2            'value of function = ',e15.8,
     3         5x,'normalization = ',e15.8)   
    4 format(/,5x,'no convergence after ',i3,' iterations')
      return
      end
