       implicit double precision (a-h,o-z)
       dimension iprm(1:100000),ns2(0:100000),ns1(0:100000)
       read (1) n2, n1, n0
       write(2,*) n2, n1, n0
       read (1) (ns2(i),i=0,n2), (ns1(i),i=0,n1), ns, mhm, nhm
       write(2,*) (ns2(i),i=0,n2), (ns1(i),i=0,n1), ns, mhm, nhm
       write(*,*) 'ns2:',(ns2(i),i=0,n2)
       write(*,*) 'ns1:',(ns1(i),i=0,n1)
       write(*,*) 'ns,mhm,nhm=',ns,mhm,nhm
       read (1) (iprm(i),i=1,mhm)
       write(2,*) (iprm(i),i=1,mhm)
*       write(*,*) (iprm(i),i=1,mhm)
*
 10    read (1) ich, jch,  ishift, jshift,idim, jdim
       write(2,*) ich, jch,  ishift, jshift,idim, jdim
       write(*,*) ich, jch,  ishift, jshift,idim, jdim
       if (ich.le.0) goto 30
 20    read (1) i,j,S
       write(2,*) i,j,S
*       STOP
       if (i.gt.0) goto 20
       if (ich.gt.0) goto 10
*
       print *,'overlap done'
       print *,'i,j,S:',i,j,S
       print *,'ich, jch,  ishift, jshift,idim, jdim:',
     >   ich, jch,  ishift, jshift,idim, jdim
*       STOP
*
*       backspace(1)
 30    read (1) ich, jch, ishift, jshift, idim, jdim 
       write(2,*) ich, jch, ishift, jshift, idim, jdim 
       write(*,*) ich, jch, ishift, jshift, idim, jdim 
       if (ich.le.0) goto 50
 40    read(1) i,j,S
       write(2,*) i,j,S 
*       print *,'2nd: i,j,S:',i,j,S 
*       stop
       if(i.gt.0) go to 40
       if (ich.gt.0) goto 30
*
 50    print *,'hamiltonian done'
*
       stop
       end
