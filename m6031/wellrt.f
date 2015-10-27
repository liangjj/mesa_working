*deck wellrt
      subroutine wellrt(e0,v0,rbox,convg,nroots,niter)
      implicit integer (a-z)
      real*8 e0, v0, rbox, gam, dgam, f, df, convg, ediff
      common /io/ inp, iout
      eold=e0
      enew=e0
      do 10 i=1,niter
         kappa=sqrt(2.d0*(abs(v0)-abs(eold)))
         fac1=tan(kappa*rbox)/kappa
         fac2=1.d0/sqrt(2.d0*abs(eold))
         dfac1=
         dfac2
                fac=1.d0
                if (eold.lt.0.d0) then
                    fac=-1.d0
                endif
                r=0.d0
                dr=0.d0 
                do 30 j=1,ncon
                   r=r+hfbox(i)*hfbox(j)/(eig(j)-eold)
                   dr=dr-hfbox(j)*hfbox(j)/((eig(j)-eold)*(eig(j)-eold))
   30           continue
                gam=sqrt(2.d0*eabs)
                dgam=fac/gam
                f=gam*r+1.d0
                df=r*dgam+dr*gam
                enew=eold-f/df
                ediff=abs(eold-enew)
                write(iout,3) i, enew, ediff
                if (ediff.le.convg) then
                    go to 40
                else
                    eold=enew
                endif
   20        continue
             write(iout,5) i
             go to 10
   40        write (iout,4) i, enew, f
         endif  
   10 continue
    1 format(/,5x,'bound eigenvalue for root = ',i3,2x,
     1            'starting at = ',e15.8)
    2 format(/,5x,'iteration',2x,'     root      ',5x,
     1            '   delta e   ',/)
    3 format(7x,i3,6x,e15.8,2x,e15.8)          
    4 format(/,5x,'eigenvalue converged after ',i3,' iterations',/5x,
     2            'value of converged eigenvalue = ',e15.8,/,5x,
     3            'value of function = ',e15.8)   
    5 format(/,5x,'no convergence after ',i3,' iterations')
      return
      end
