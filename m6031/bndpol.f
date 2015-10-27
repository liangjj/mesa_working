*deck bndpol
      subroutine bndpol(eig,hfbox,convg,nroots,ncon,niter)
      implicit integer (a-z)
      real*8 eig, hfbox, eold, enew, eabs, fac, r, dr, gam
      real*8 dgam, f, df, convg, ediff, fuzz
      dimension eig(ncon), hfbox(ncon)
      common /io/ inp, iout
      data fuzz/1.d-03/
      do 10 iene=1,nroots
         if (eig(iene).lt.0.d0) then
             write(iout,1) iene, eig(iene)
             write(iout,2)
             eold=eig(iene)+fuzz
             enew=eold
             do 20 i=1,niter
                eabs=abs(eold)
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
