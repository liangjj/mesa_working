#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "f2c.h"

#define pi 3.14159265358979323846
#define NR_END 1
#define FREE_ARG char*
#define TOL 0.00000000001 /* rms error of the wavefunction */

#define PIM4 0.7511255444649425  /* 1/pi^(1/4) - normalization for H_n */

#define Nmin 1000  /* minimum number of atoms in the trap */
#define Nmax 1000  /* maximum number of atoms in the trap */
#define Nint 1000  /* interval number of atoms in loop */
#define tratio 0.001  /* the ratio of temperature to chemical potential */
#define a0 0.000000275  /* s-wave scattering length, 23Na */
/* #define a0 0.000000573 s-wave scattering length, Rb */
/* #define a0 0.00000031482837 s-wave scattering length, Cs */
/* #define d0 0.0002762 oscillator length of the trap along x, Cs */
#define d0 0.0005635 /* approximate length of Lene Hau's trap along x, 23Na */
/* #define d0 0.0005859 actual length of Lene Hau's trap along x, 23Na */
/* #define d0 0.0004774  axial oscillator length of Ketterle's trap, 23Na */
/* #define d0 0.00015761  (largest) oscillator length of Phillips' trap, 23Na */
/* #define alpha 1.414213562  omega_y/omega_x, Phillips' trap = sqrt(2) */
/* #define beta 2.0  omega_z/omega_x, Phillips' trap */
/* #define alpha 1.162162162  omega_y/omega_x, Lene' trap */
/* #define beta 1.162162162  omega_z/omega_x, Lene' trap */
/* #define alpha 12.9668  omega_y/omega_x, Ketterle's trap */
/* #define beta 12.9668  omega_z/omega_x, Ketterle's trap  */
#define alpha 1.0  /* omega_y/omega_x, isotropic case */
#define beta 1.0  /* omega_z/omega_x, isotropic case */
#define bi_max 1  /* number of times to ramp down the kinetic contribution */
#define kmax 1.0 /* maximum factor to multiply kinetic energy */
#define numgrid 2  /* number of grids in multigrid procedure */
#define Emin 0.0  /* minimum excitation energy relative to mu */
#define Emax 10.0  /* maximum excitation energy relative to mu */
#define px 0  /* parity state along x for the excitations */
#define py 0  /* parity state along x for the excitations */
#define pz 0  /* parity state along x for the excitations */

double **hermit(double *, long);
/* calculates the Hermite polynomials at a given list of points */
double radius(double *, double *, double *, long, long, long);
/* calculates the radius associated with and (x,y,z) values */
void print_psi(double *,double *,double *,double *,double, int, int, int, int,
	int, double);
/* prints out the wavefunction along the three Cartesian axes */
void honv2(int *,int *,int *,int *,double *, int *, int *,int *,double *,
	doublecomplex *,doublecomplex *,int *,int *);
/* multiplies the Bogoliubov matrix on a vector */
int *ivector(long,long);
double *fvector(long,long);
double *dvector(long,long);
double **fmatrix(long,long,long,long);
void free_ivector(int *, long, long);
void free_fvector(double *, long, long);
void free_dvector(double *, long, long);
void free_fmatrix(double **, long, long, long, long);
doublecomplex *cvector(long,long);
doublecomplex *zvector(long,long);
void free_cvector(doublecomplex *, long, long);
void free_zvector(doublecomplex *, long, long);

double **hx;
double **hy;
double **hz;
double **hxo;
double **hyo;
double **hzo;

FILE *tout;
FILE *tin;
FILE *fout;
FILE *fout2;
FILE *fin;
FILE *fin2;

main()
{
	int i,bi,j,k,kk,n,nn,ix,iy,iz,jx,jy,jz,kx,ky,kz,size,cols,dim,dum;
	int info,Np,Nx,Ny,Nz,done,maxsize,imax,grid,factor;
	int *ia, *ja, *ipiv,N1[numgrid+1],N2[numgrid+1],N3[numgrid+1];
	int maxi[numgrid+1], bi_low[numgrid+1], bi_high[numgrid+1];
	int nnt,Nxt,Nyt,Nzt;
	double mu,mu_old,U,norman,lambda,mu0,eta0,R0,En,test,Ubar,Tbar,gamma;
	double *xx, *xy, *xz, *weightx, *weighty, *weightz, *w, *w_good, *sa;
	double *xxo, *xyo, *xzo, *wxo, *wyo, *wzo, gam[numgrid+1];
	double dum1,dum2,dum3,dum4,rmax,rmaxo;
	double **p, **pa, *psi, *psi2, *mut, *out, **B, *C;
	double *rho_ab, *rho_yz, *psir, **phix, **phiy, **phiz, gammat;
	int Q,R,PINIT,MMAX,M;
	long entries,dim_alloc,dim_alloco,order,K_elem,V_elem;
	double EPS,D[20];
	char filename[100],waver[100];
/* parameters for the excitations */
	int e_index, sc_len, sb_len, maxit, n_eigs, order_old;
	int *ib, *jb, *ic, *jc, *id, *jd;
	double *sc, *sd, **rvec, **wer, *eigs_exact;
	double *e_normal, scale, toler;
	static logical drctv, nonsym;
	doublecomplex *vex, *eigs, *diag;
	/* parameters for the LAPACK routines */
	int jobl, jobe;
	int lda,ldvl,ldvr,lwork;
	double *sb, *work;
	doublecomplex *wr, *vr;

	printf("hello\n");

/* We have assumed that the Hamiltonian matrix is extremely sparse in the DVR
representation. In fact, for Nx=Ny=Nz=60, less than .1% of the matrix elements
are non-zero */

	fout=fopen("tempout","w");

/* parameters for the Lanczos routine */

  Q=20; /* number of vectors to keep in the diagonalisation */
  PINIT=6; /* initial block size */
  R=1; /* the number of eigenvalues and eigenvectors to calculate */
  MMAX=8000000; /* The maximum number of matrix/vector products */
  EPS=0.00000000000000001; /* relative precision for eigenvalue */

	gam[1]=0.1;  /* coarsening factor <= 1; the smaller, the more coarse */
	N1[1]=20;   /* number of basis functions along x */
	N2[1]=20;   /* number of basis functions along y */
	N3[1]=20;   /* number of basis functions along z */
	maxi[1]=300;  /* maximum number of iterations for self-consistency */
	bi_low[1]=bi_max;
	bi_high[1]=bi_max;
	gam[2]=0.1; 
	N1[2]=30;  
	N2[2]=30;  
	N3[2]=30; 
	maxi[2]=300; 
	bi_low[2]=bi_max;
	bi_high[2]=bi_max;
	/* gam[3]=0.2;
	N1[3]=80; 
	N2[3]=80;
	N3[3]=80;
	maxi[3]=500; 
	bi_low[3]=bi_max;
	bi_high[3]=bi_max;
	gam[4]=0.3; 
	N1[4]=80;  
	N2[4]=80; 
	N3[4]=80;
	maxi[4]=300; 
	bi_low[4]=bi_max;
	bi_high[4]=bi_max;
	gam[5]=1.0; 
	N1[5]=80;  
	N2[5]=80; 
	N3[5]=80;
	maxi[5]=100; 
	bi_low[5]=bi_max;
	bi_high[5]=bi_max; */

	for (nn=Nmin;nn<=Nmax;nn+=Nint) { /* loop over number of atoms */

		for (grid=1;grid<=numgrid;grid++) {
			
/* loop over kinetic contribution */
			for (bi=bi_low[grid];bi<=bi_high[grid];bi++) {

				size=1;
				done=0;  /* flag for printing out the wavefunction */

/* for the first iteration in the loop over bases, use the TF expression
for the initial guess. Subsequently, use the last self-consistent result. */

				if (grid == 1 && bi == bi_low[1]) { 
					Nx=N1[1];
					Ny=N2[1];
					Nz=N3[1];
					imax=maxi[1];
					maxsize=N1[1]*N2[1]*N3[1]/8;
					gamma=gam[1];
					printf("Nx,Ny,Nz = %d,%d,%d\n",Nx,Ny,Nz);
					printf("gamma = %lf\n",gamma);
					factor=(int)(1.0/gamma);
					/* entries=(maxsize/200)*(maxsize/factor); */
					entries=(maxsize)*(maxsize);
					/* dim_alloc=maxsize/factor; */
					dim_alloc=maxsize/1;
					printf("allocated entries = %ld\n",entries);
					printf("allocated dimension of H = %ld\n",dim_alloc);

					Np=nn;
					/* eta0=Np*a0/d0*pow(gamma,1.5)*pow(alpha*beta,1.0/6.0); */
					eta0=Np*a0/d0*pow(gamma,1.5);
					U=4.0*pi*eta0;  /* interaction strength */

					ia=ivector(0,entries);
					ja=ivector(0,entries);
					sa=fvector(0,entries);
					w=fvector(0,Q*dim_alloc);
					w_good=fvector(0,Q*dim_alloc);
					xx=fvector(1,Nx);
					xy=fvector(1,Ny);
					xz=fvector(1,Nz);
					weightx=fvector(1,Nx);
					weighty=fvector(1,Ny);
					weightz=fvector(1,Nz);
					hx=fmatrix(0,Nx,0,Nx);
					hy=fmatrix(0,Ny,0,Ny);
					hz=fmatrix(0,Nz,0,Nz);
					p=fmatrix(0,imax+1,0,dim_alloc);
					pa=fmatrix(0,imax+1,0,dim_alloc);
					psi=fvector(0,dim_alloc);
					psi2=fvector(0,dim_alloc);
					mut=fvector(0,imax+1);

/* compute the roots and weights of Nth order Hermite polynomial. */
					gauher(xx,weightx,Nx);
					gauher(xy,weighty,Ny);
					gauher(xz,weightz,Nz);
				
/* store the values of the normalized Hermite polynomials h_(1..N)[x(i)] at 
						the roots x[i]. */
					hx=hermit(xx,Nx);
					hy=hermit(xy,Ny);
					hz=hermit(xz,Nz);

					rmax=sqrt((xx[1]*xx[1]+xy[1]*xy[1]+xz[1]*xz[1])/3.0*gamma);
					printf("rmax = %lf\n",rmax);
					R0=pow(15.0*eta0*alpha*beta,0.2);
					/* mu0=13.0+0.5*R0*R0/sqrt(gamma); */
					mu0=5.0+0.5*R0*R0/sqrt(gamma);
					mu=mu0;
					printf("mu = %lf\n",mu);

					sprintf(waver,"psi_temp_%d.out",nn);
					fin=fopen(waver,"r");
					if (fin != NULL) {
						fscanf(fin,"%d  %d  %d  %d  %lf  %d\n",&nnt,&Nxt,&Nyt,&Nzt,&gammat,
							&dim);
						if (nnt != nn || Nxt != Nx || Nyt != Ny || Nzt != Nz
							|| gammat != gamma) {
							for (i=2;i<=numgrid;i++) {
								if (gammat == gam[i] && Nxt == N1[i] && Nyt == N2[i] 
									&& Nzt == N3[i] && nnt == nn) {
									grid=i;
									bi=bi_max;
									break;
								}
								if (i == numgrid) {
									printf("parameters inconsistent with file\n");
									exit(1);
								}
							}
						}
						for (i=0;i<dim;i++) {
							fscanf(fin,"%lf\n",&psi[i]);
							psi2[i]=psi[i]*psi[i];
						}
					}

					if (fin == NULL) {
						dim=0;
						norman=0.0;
						for (ix=1;ix<=Nx/2;ix++) {  /* initial conditions */
							for (iy=1;iy<=Ny/2;iy++) { 
								for (iz=1;iz<=Nz/2;iz++) {
									if (sqrt(xx[ix]*xx[ix]+xy[iy]*xy[iy]+xz[iz]*xz[iz]) > rmax)
										continue;
									/* psi2[dim]=0.0; */
									psi2[dim]=(mu0-0.5*(xx[ix]*xx[ix]+xy[iy]*xy[iy]*alpha
										+xz[iz]*xz[iz]*beta)/gamma)/U;
									if (psi2[dim] < 0.0) psi2[dim]=0.0;
									psi[dim]=sqrt(psi2[dim]);
									norman += psi2[dim];
									++dim;
								}
							}
						}
						printf("dim(Psi) = %d, norm = %lf\n",dim,norman);
						/* print_psi(psi,xx,xy,xz,rmax,done,Nx,Ny,Nz,maxsize,gamma); */
						/* for (i=0;i<dim;i++) psi2[dim] /= norman; */
					}
					fclose(fin);
				}

				if (grid > 1 && bi == bi_low[grid]) { 
					free_ivector(ia,0,entries);
					free_ivector(ja,0,entries);
					free_fvector(sa,0,entries);
					free_fmatrix(p,0,imax+1,0,dim_alloc);
					free_fmatrix(pa,0,imax+1,0,dim_alloc);
					free_fvector(psi,0,dim_alloc);
					free_fvector(psi2,0,dim_alloc);
					free_fvector(mut,0,imax+1);

					phix=fmatrix(0,N1[grid]/2,0,Nx/2);
					phiy=fmatrix(0,N2[grid]/2,0,Ny/2);
					phiz=fmatrix(0,N3[grid]/2,0,Nz/2);
					xxo=fvector(1,Nx);
					xyo=fvector(1,Ny);
					xzo=fvector(1,Nz);
					wxo=fvector(1,Nx);
					wyo=fvector(1,Ny);
					wzo=fvector(1,Nz);
					hxo=fmatrix(0,Nx,0,Nx);
					hyo=fmatrix(0,Ny,0,Ny);
					hzo=fmatrix(0,Nz,0,Nz);

					for (i=1;i<=Nx;i++) {
						xxo[i]=xx[i];
						wxo[i]=weightx[i];
						for (k=1;k<=Nx;k++) hxo[k][i-1]=hx[k][i-1];
					}
					for (i=1;i<=Ny;i++) {
						xyo[i]=xy[i];
						wyo[i]=weighty[i];
						for (k=1;k<=Ny;k++) hyo[k][i-1]=hy[k][i-1];
					}
					for (i=1;i<=Nz;i++) {
						xzo[i]=xz[i];
						wzo[i]=weightz[i];
						for (k=1;k<=Nz;k++) hzo[k][i-1]=hz[k][i-1];
					}
					rmaxo=rmax;
					dim_alloco=dim_alloc;

					free_fvector(xx,1,Nx);
					free_fvector(xy,1,Ny);
					free_fvector(xz,1,Nz);
					free_fvector(weightx,1,Nx);
					free_fvector(weighty,1,Ny);
					free_fvector(weightz,1,Nz);
					free_fmatrix(hx,0,Nx,0,Nx);
					free_fmatrix(hy,0,Ny,0,Ny);
					free_fmatrix(hz,0,Nz,0,Nz);

					Nx=N1[grid];
					Ny=N2[grid];
					Nz=N3[grid];
					imax=maxi[grid];
					maxsize=Nx*Ny*Nz/8;
					gamma=gam[grid];
					printf("Nx,Ny,Nz = %d,%d,%d\n",Nx,Ny,Nz);
					printf("gamma = %lf\n",gamma);

					factor=(int)(1.0/gamma);
					/* entries=(maxsize/200)*(maxsize/factor); */
					entries=(maxsize)*(maxsize);
					/* dim_alloc=maxsize/factor; */
					dim_alloc=maxsize/1;
					printf("allocated entries in H = %ld\n",entries);
					printf("allocated dimension of H = %ld\n",dim_alloc);

					Np=nn;
					eta0=Np*a0/d0*pow(gamma,1.5);
					U=4.0*pi*eta0;  /* interaction strength */

					ia=ivector(0,entries);
					ja=ivector(0,entries);
					sa=fvector(0,entries);
					p=fmatrix(0,imax+1,0,dim_alloc);
					pa=fmatrix(0,imax+1,0,dim_alloc);
					psi=fvector(0,dim_alloc);
					psi2=fvector(0,dim_alloc);
					mut=fvector(0,imax+1);
					wer=fmatrix(0,e_index,0,dim_alloc);
					e_normal=fvector(0,e_index);

					xx=fvector(1,Nx);
					xy=fvector(1,Ny);
					xz=fvector(1,Nz);
					weightx=fvector(1,Nx);
					weighty=fvector(1,Ny);
					weightz=fvector(1,Nz);
					hx=fmatrix(0,Nx,0,Nx);
					hy=fmatrix(0,Ny,0,Ny);
					hz=fmatrix(0,Nz,0,Nz);

/* compute the roots and weights of Nth order Hermite polynomial. */
					gauher(xx,weightx,Nx);
					gauher(xy,weighty,Ny);
					gauher(xz,weightz,Nz);

					for (i=1;i<=Nx;i++) xx[i] *= sqrt(gam[grid-1]/gamma);
					for (i=1;i<=Ny;i++) xy[i] *= sqrt(gam[grid-1]/gamma);
					for (i=1;i<=Nz;i++) xz[i] *= sqrt(gam[grid-1]/gamma);

/* store the values of the normalized Hermite polynomials h_(1..N)[x(i)] at 
						the NEW roots x[i]. */
					hx=hermit(xx,Nx);
					hy=hermit(xy,Ny);
					hz=hermit(xz,Nz);

					for (i=1;i<=Nx;i++) xx[i] /= sqrt(gam[grid-1]/gamma);
					for (i=1;i<=Ny;i++) xy[i] /= sqrt(gam[grid-1]/gamma);
					for (i=1;i<=Nz;i++) xz[i] /= sqrt(gam[grid-1]/gamma);

					rmax=sqrt((xx[1]*xx[1]+xy[1]*xy[1]+xz[1]*xz[1])/3.0*gamma);

/* First generate the Lagrange polynomials phi on the interpolated grid */

					for (i=0;i<Nx/2;i++) {
						for (j=0;j<N1[grid-1]/2;j++) {
							phix[i][j]=0.0;
							for (k=0;k<N1[grid-1];k+=2) {
								phix[i][j] += sqrt(2.0*wxo[j+1])*hx[k+1][i]*hxo[k+1][j];
							}
						}
					}
					for (i=0;i<Ny/2;i++) {
						for (j=0;j<N2[grid-1]/2;j++) {
							phiy[i][j]=0.0;
							for (k=0;k<N2[grid-1];k+=2) {
								phiy[i][j] += sqrt(2.0*wyo[j+1])*hy[k+1][i]*hyo[k+1][j];
							}
						}
					}
					for (i=0;i<Nz/2;i++) {
						for (j=0;j<N3[grid-1]/2;j++) {
							phiz[i][j]=0.0;
							for (k=0;k<N3[grid-1];k+=2) {
								phiz[i][j] += sqrt(2.0*wzo[j+1])*hz[k+1][i]*hzo[k+1][j];
							}
						}
					}

/* Second generate the interpolated wavefunction(s) on the new grid */
				
					psir=fvector(0,dim_alloc);

					dim=0;
					for (i=0;i<maxsize;i++) {
						ix=i/(Ny*Nz/4);
						iy=(i-ix*(Ny*Nz/4))/(Nz/2);
						iz=i-iy*(Nz/2)-ix*(Ny*Nz/4);
						if (radius(xx,xy,xz,ix,iy,iz) > rmax) continue;
						psi[dim]=0.0;
						psir[dim]=0.0;
						for (k=0;k<e_index;k++) wer[k][dim]=0.0;
						dum=0;
						for (j=0;j<N1[grid-1]*N2[grid-1]*N3[grid-1]/8;j++) {
							jx=j/(N2[grid-1]*N3[grid-1]/4);
							jy=(j-jx*(N2[grid-1]*N3[grid-1]/4))/(N3[grid-1]/2);
							jz=j-jy*(N3[grid-1]/2)-jx*(N2[grid-1]*N3[grid-1]/4);
							if (radius(xxo,xyo,xzo,jx,jy,jz) > rmaxo) continue;
							psi[dim] += pow(gam[grid-1]/gamma,0.75)*w[dum]
								*exp(-0.5*(xx[ix+1]*xx[ix+1]+xy[iy+1]*xy[iy+1]
								+xz[iz+1]*xz[iz+1])*gam[grid-1]/gamma)
								*phix[ix][jx]*phiy[iy][jy]*phiz[iz][jz];
							psir[dim] += pow(gam[grid-1]/gamma,0.75)*w[dum]
								*exp(-0.5*(xx[ix+1]*xx[ix+1]+xy[iy+1]*xy[iy+1]
								+xz[iz+1]*xz[iz+1])*gam[grid-1]/gamma/gamma)
								*phix[ix][jx]*phiy[iy][jy]*phiz[iz][jz];
							for (k=0;k<e_index;k++) {
								wer[k][dim] += pow(gam[grid-1]/gamma,0.75)*rvec[k][dum]
									*exp(0.5*(1.0-gam[grid-1]/gamma)*(xx[ix+1]*xx[ix+1]
									+xy[iy+1]*xy[iy+1]+xz[iz+1]*xz[iz+1]))
									*sqrt(8.0*weightx[ix+1]*weighty[iy+1]*weightz[iz+1])
									*phix[ix][jx]*phiy[iy][jy]*phiz[iz][jz];
							}
							++dum;
						}
						psi2[dim]=psi[dim]*psi[dim];
						++dim;
					}

					free_fmatrix(rvec,0,3*(Emax*Emax-Emin*Emin),0,order);

					vex=cvector(0,e_index*dim);
					for (j=0;j<e_index;j++) {
						for (i=0;i<dim;i++) {
							vex[j*dim+i].r=wer[j][i];
							vex[j*dim+i].i=0.0;
						}
					}

					free_fmatrix(phix,0,Nx/2,0,N1[grid-1]/2);
					free_fmatrix(phiy,0,Ny/2,0,N2[grid-1]/2);
					free_fmatrix(phiz,0,Nz/2,0,N3[grid-1]/2);
					free_fvector(xxo,1,N1[grid-1]);
					free_fvector(xyo,1,N2[grid-1]);
					free_fvector(xzo,1,N3[grid-1]);
					free_fvector(wxo,1,N1[grid-1]);
					free_fvector(wyo,1,N2[grid-1]);
					free_fvector(wzo,1,N3[grid-1]);
					free_fmatrix(hxo,0,N1[grid-1],0,N1[grid-1]);
					free_fmatrix(hyo,0,N2[grid-1],0,N2[grid-1]);
					free_fmatrix(hzo,0,N3[grid-1],0,N3[grid-1]);
					free_fvector(w,0,Q*dim_alloco);
					free_fvector(w_good,0,Q*dim_alloco);
					free_fmatrix(wer,0,e_index,0,dim_alloc);

					w=fvector(0,Q*dim_alloc);
					w_good=fvector(0,Q*dim_alloc);
					hx=hermit(xx,Nx);
					hy=hermit(xy,Ny);
					hz=hermit(xz,Nz);

					done=1;
					print_psi(psir,xx,xy,xz,rmax,done,Nx,Ny,Nz,maxsize,gamma);
					done=0;

					free_fvector(psir,0,dim_alloc);

					/* if (bi == bi_max && grid == numgrid) {
						sprintf(waver,"psi_temp_%5.3lf_%d.out",gamma,nn);
						fout2=fopen(waver,"w");
						fprintf(fout2,"%d  %d  %d  %d  %lf  %d\n",nn,Nx,Ny,Nz,gamma,dim);
						for (i=0;i<dim;i++) fprintf(fout2,"%30.25lf\n",psi[i]);
						fclose(fout2);
						break;
					} */
				}

/* Calculate the kinetic energy matrix, which only has to be done once. First,
check whether this has already been done. If not, do so and save the results
to a file. */

				sprintf(filename,"M_%5.3f_%d_%d_%d",gam[grid],Nx,Ny,Nz);
				tin=fopen(filename,"r");

/* test whether the file containing the elements of T already exists */

				if (tin == NULL) { /* if it doesn't, generate and save it */
					tout=fopen(filename,"w");
					dim=0;
					K_elem=0;
					for (i=0;i<maxsize;i++) { /* loop over positive coordinates only */
						ix=i/(Ny*Nz/4);
						iy=(i-ix*(Ny*Nz/4))/(Nz/2);
						iz=i-iy*(Nz/2)-ix*(Ny*Nz/4);
						if (radius(xx,xy,xz,ix,iy,iz) > rmax) continue;
						sa[K_elem]=0.5*gamma*(1.0+alpha+beta-(xx[ix+1]*xx[ix+1]
							+alpha*xy[iy+1]*xy[iy+1]+beta*xz[iz+1]*xz[iz+1]));
						ia[K_elem]=dim+1;
						ja[K_elem]=dim+1;
						++K_elem;
						dum=dim;
						for (j=i;j<maxsize;j++) {
							jx=j/(Ny*Nz/4);
							jy=(j-jx*(Ny*Nz/4))/(Nz/2);
							jz=j-jy*(Nz/2)-jx*(Ny*Nz/4);
							if (radius(xx,xy,xz,jx,jy,jz) > rmax) continue;
							if (iy == jy && iz == jz) {
								sa[K_elem]=0.0;
								for (k=0;k<Nx;k+=2) {
									sa[K_elem] += 2.0*sqrt(weightx[ix+1])*sqrt(weightx[jx+1])
										*gamma*(double)k*hx[k+1][ix]*hx[k+1][jx];
								}
								ia[K_elem]=dim+1;
								ja[K_elem]=dum+1;
								++K_elem;
							}
							if (ix == jx && iz == jz) {
								sa[K_elem]=0.0;
								for (k=0;k<Ny;k+=2) {
									sa[K_elem]+=2.0*alpha*sqrt(weighty[iy+1])*sqrt(weighty[jy+1])
										*gamma*(double)k*hy[k+1][iy]*hy[k+1][jy];
								}
								ia[K_elem]=dim+1;
								ja[K_elem]=dum+1;
								++K_elem;
							}
							if (ix == jx && iy == jy) {
								sa[K_elem]=0.0;
								for (k=0;k<Nz;k+=2) {
									sa[K_elem] += 2.0*beta*sqrt(weightz[iz+1])*sqrt(weightz[jz+1])
										*gamma*(double)k*hz[k+1][iz]*hz[k+1][jz];
								}
								ia[K_elem]=dim+1;
								ja[K_elem]=dum+1;
								++K_elem;
							}
							++dum;
						}
						++dim;
					}
					printf("Non-zero elements in T = %d\n",K_elem);
					fprintf(tout,"%d\n",K_elem);
					for (i=0;i<K_elem;i++) fprintf(tout,"%d  %d  %20.12lf\n",ia[i],ja[i],
						sa[i]);
				}
				fclose(tout);

				if (tin != NULL) { /* if it does, read it */
					fscanf(tin,"%d\n",&K_elem);
					for (i=0;i<K_elem;i++) fscanf(tin,"%d  %d  %lf\n",&ia[i],&ja[i],
						&sa[i]);
				}
				fclose(tin);

				for (i=0;i<K_elem;i++) {
					sa[i] *= kmax/(1.0+(kmax-1.0)*(double)bi/(double)bi_max);
				}

				V_elem=0; /* generate the (diagonal) potential energy matrix */
				dim=0;
				for (i=0;i<maxsize;i++) {
					ix=i/(Ny*Nz/4);
					iy=(i-ix*(Ny*Nz/4))/(Nz/2);
					iz=i-iy*(Nz/2)-ix*(Ny*Nz/4);
					if (radius(xx,xy,xz,ix,iy,iz) > rmax) continue;

/* neglect points where the wavefunction is identically zero, i.e. beyond
containing sphere */

					sa[K_elem+V_elem]=U*psi2[dim]+0.5*(xx[ix+1]*xx[ix+1]
						+alpha*xy[iy+1]*xy[iy+1]+beta*xz[iz+1]*xz[iz+1])/gamma;
					ia[K_elem+V_elem]=dim+1;
					ja[K_elem+V_elem]=dim+1;
					++V_elem;
					++dim;
				}
				printf("dimension = %d, non-zero elements in V = %d\n",dim,V_elem);

				order=dim;
				lda=dim;
				entries=K_elem+V_elem;
				printf("order = %d, non-zero elements = %d\n",order,entries);

				M=0; /* number of eigenvalues already computed */
				printf("loop 0\n");

				minval(&entries,ia,ja,sa,&order,&Q,&PINIT,&R,&MMAX,&EPS,&M,D,w,&info);

				printf("info = %d\n",info);
				for (i=0;i<R;i++) printf("D[%d]=%lf\n",i,D[i]);
				mut[0]=D[0];
				test=0.0;
				dim=0;
				for (j=0;j<maxsize;j++) {
					jx=j/(Ny*Nz/4);
					jy=(j-jx*(Ny*Nz/4))/(Nz/2);
					jz=j-jy*(Nz/2)-jx*(Ny*Nz/4);
					if (radius(xx,xy,xz,jx,jy,jz) > rmax) continue;
					psi[dim]=w[dim]/sqrt(8.0*weightx[jx+1]*weighty[jy+1]*weightz[jz+1])
						*exp(-0.5*(xx[jx+1]*xx[jx+1]+xy[jy+1]*xy[jy+1]+xz[jz+1]*xz[jz+1]));
					psi2[dim]=psi[dim]*psi[dim];
					p[0][dim]=w[dim]; /* store the resulting eigenvector */
					test += w[dim]*w[dim];
					++dim;
				}

				B=fmatrix(1,imax+1,1,imax+1); /* arrays for the DIIS routine */
				out=fvector(0,imax);

				for (n=1;n<=imax;n++) { /* start the self-consistency loop */

					M=0; /* number of eigenvalues already computed */
					printf("loop %d\n",n);

					V_elem=0; /* generate the (diagonal) potential energy matrix */
					dim=0;
					for (i=0;i<maxsize;i++) {
						ix=i/(Ny*Nz/4);
						iy=(i-ix*(Ny*Nz/4))/(Nz/2);
						iz=i-iy*(Nz/2)-ix*(Ny*Nz/4);
						if (radius(xx,xy,xz,ix,iy,iz) > rmax) continue;
						sa[K_elem+V_elem]=U*psi2[dim]+0.5*(xx[ix+1]*xx[ix+1]
							+alpha*xy[iy+1]*xy[iy+1]+beta*xz[iz+1]*xz[iz+1])/gamma;
						ia[K_elem+V_elem]=dim+1;
						ja[K_elem+V_elem]=dim+1;
						++V_elem;
						++dim;
					}

/* Begin the DIIS procedure */

					dim=0;
					for (i=0;i<maxsize;i++) {
						ix=i/(Ny*Nz/4);
						iy=(i-ix*(Ny*Nz/4))/(Nz/2);
						iz=i-iy*(Nz/2)-ix*(Ny*Nz/4);
						if (radius(xx,xy,xz,ix,iy,iz) > rmax) continue;
						pa[n-1][dim]=0.0; /* initialize alternate vector, H*w  */
						++dim;
					}
					for (j=0;j<K_elem+V_elem;j++) {
						pa[n-1][ia[j]-1] += sa[j]*p[n-1][ja[j]-1];
						if (ia[j] != ja[j]) pa[n-1][ja[j]-1] += sa[j]*p[n-1][ia[j]-1];
					}
					if (size < imax) ++size;

					for (k=1;k<=size-1;k++) { /* generate the matrix elements Bij */
						out[k-1]=0.0;
						dum1=dum2=dum3=dum4=0.0;
						for (i=0;i<dim;i++) {
							dum1 += pa[n-size+k][i]*pa[n-1][i];
							dum2 += p[n-size+k][i]*p[n-1][i];
							dum3 += pa[n-size+k][i]*p[n-1][i];
							dum4 += p[n-size+k][i]*pa[n-1][i];
						}
						B[size-1][k]=2.0*(dum1*dum2-dum3*dum4); /* trace of e(n)*e(k-1) */
						B[k][size-1]=B[size-1][k];
						B[size][k]=-1.0;
						B[k][size]=-1.0;
					}
					B[size][size]=0.0;
					out[size-1]=-1.0;
					cols=1;
					C=fvector(0,size*size);
					for (i=0;i<size;i++) {
						for (j=0;j<size;j++) {
							C[size*i+j]=B[i+1][j+1];
						}
					}
					ipiv=ivector(0,size);
					jobl=0;

					/* dgesv(&size,&cols,C,&size,ipiv,out,&size,&info); */
					sgefa(C,&size,&size,ipiv,&info);
					sgesl(C,&size,&size,ipiv,out,jobl);

					free_ivector(ipiv,0,size);
					free_fvector(C,0,size*size);
					/* for (i=0;i<size;i++) printf("%d  %lf\n",i+1,out[i]); */
					printf("grid %d: out[%d] = %le\n",grid,size,out[size-1]);
			
					dim=0;
					for (i=0;i<maxsize;i++) {
						ix=i/(Ny*Nz/4);
						iy=(i-ix*(Ny*Nz/4))/(Nz/2);
						iz=i-iy*(Nz/2)-ix*(Ny*Nz/4);
						if (radius(xx,xy,xz,ix,iy,iz) > rmax) continue;
						psi2[dim]=0.0;
						for (j=1;j<=size-1;j++) {
							psi2[dim] += out[j-1]*p[n-size+j][dim]*p[n-size+j][dim]
								/(8.0*weightx[ix+1]*weighty[iy+1]*weightz[iz+1])
								*exp(-(xx[ix+1]*xx[ix+1]+xy[iy+1]*xy[iy+1]+xz[iz+1]*xz[iz+1]));
						}
						++dim;
					}
					mu_old=mu;
					mu=0.0;
					for (j=1;j<=size-1;j++) mu += out[j-1]*mut[n-size+j];
					printf("bi = %d, gamma = %lf, mu = %lf\n",bi,gamma,mu);

/* End the DIIS procedure */

					V_elem=0; /* generate the (diagonal) potential energy matrix */
					dim=0;
					for (i=0;i<maxsize;i++) {
						ix=i/(Ny*Nz/4);
						iy=(i-ix*(Ny*Nz/4))/(Nz/2);
						iz=i-iy*(Nz/2)-ix*(Ny*Nz/4);
						if (radius(xx,xy,xz,ix,iy,iz) > rmax) continue;
						sa[K_elem+V_elem]=U*psi2[dim]+0.5*(xx[ix+1]*xx[ix+1]
							+alpha*xy[iy+1]*xy[iy+1]+beta*xz[iz+1]*xz[iz+1])/gamma;
						ia[K_elem+V_elem]=dim+1;
						ja[K_elem+V_elem]=dim+1;
						++V_elem;
						++dim;
					}

					order=dim;
					lda=dim;
					psir=fvector(0,order);
					entries=K_elem+V_elem;
					minval(&entries,ia,ja,sa,&order,&Q,&PINIT,&R,&MMAX,&EPS,&M,D,w,
						&info);
					printf("info = %d\n",info);
					for (i=0;i<R;i++) printf("D[%d]=%lf\n",i,D[i]);
					mut[n]=D[0];
					test=0.0;
					dim=0;
					for (j=0;j<maxsize;j++) {
						jx=j/(Ny*Nz/4);
						jy=(j-jx*(Ny*Nz/4))/(Nz/2);
						jz=j-jy*(Nz/2)-jx*(Ny*Nz/4);
						if (radius(xx,xy,xz,jx,jy,jz) > rmax) continue;
						p[n][dim]=w[dim]; /* store the resulting eigenvector */
						psi[dim]=w[dim]/sqrt(8.0*weightx[jx+1]*weighty[jy+1]*weightz[jz+1])
						 *exp(-0.5*(xx[jx+1]*xx[jx+1]+xy[jy+1]*xy[jy+1]+xz[jz+1]*xz[jz+1]));
						psir[dim]=w[dim]/sqrt(8.0*weightx[jx+1]*weighty[jy+1]*weightz[jz+1])
							*exp(-0.5*(xx[jx+1]*xx[jx+1]+xy[jy+1]*xy[jy+1]+xz[jz+1]*xz[jz+1])
							/gamma);
						psi2[dim]=psi[dim]*psi[dim];
						test += w[dim]*w[dim];
						++dim;
					}
					/* printf("test = %lf\n",test); */

					/* print_psi(psir,xx,xy,xz,rmax,done,Nx,Ny,Nz,maxsize,gamma); */

/* check whether convergence has been reached. The criterion is the norm
of the error vector */

					if (fabs(out[size-1]) < TOL && n > 2) {
						done=1;
						Ubar=0.0;
						for (j=0;j<dim;j++) Ubar += U*w[j]*w[j]*psi2[j];
						En=mu-0.5*Ubar;
						printf("%d  %lf  %lf  %lf  %lf  %d\n",Np,eta0,En,mu,
							fabs(out[size-1]),n);
						fprintf(fout,"%d  %lf  %lf  %lf  %lf  %d\n",Np,eta0,En,mu,
							fabs(out[size-1]),n);
						print_psi(psir,xx,xy,xz,rmax,done,Nx,Ny,Nz,maxsize,gamma);
						free_fvector(psir,0,order);

						/* sprintf(waver,"psi_temp_%d.out",nn);
						fout2=fopen(waver,"w");
						fprintf(fout2,"%d  %d  %d  %d  %lf  %d\n",nn,Nx,Ny,Nz,gamma,dim);
						for (i=0;i<order;i++) fprintf(fout2,"%30.25lf\n",psi[i]);
						fclose(fout2); */
						for (i=0;i<order;i++) w_good[i]=w[i];

/*  Now calculate the excitations! */

/* First obtain spectrum by exact diagonalization on the coarse grid. Once
interpolated, the left and right eigenvectors will be good approximations
(preconditioners) for the fine-grid excitations. */

						if (grid == 1 && bi == bi_max) {
							rvec=fmatrix(0,3*(Emax*Emax-Emin*Emin),0,order);
							/* sprintf(waver,"excite_temp_%5.3lf_%d_%d,%d,%d.out",gamma,nn,
								px,py,pz);
							fin2=fopen(waver,"r");
							if (fin2 != NULL) {
								fscanf(fin2,"%d  %d\n",&e_index,&order);
								for (i=0;i<e_index;i++) {
									for (j=0;j<order;j++) {
										fscanf(fin2,"%lf\n",&rvec[i][j]);
									}
								}
								fclose(fin2);
								break;
							} */

							jobe=1;
							lda=order;
							ldvr=order;
							wr=cvector(0,order);
							vr=cvector(0,ldvr*order);
							work=fvector(0,2*order);
							id=ivector(0,V_elem);
							jd=ivector(0,V_elem);
							sd=fvector(0,V_elem);
							sc=dvector(0,order*order);
							sb=dvector(0,order*order);

							for (i=0;i<V_elem;i++) {
								id[i]=ia[K_elem+i];
								jd[i]=ja[K_elem+i];
								sd[i]=sa[K_elem+i]-mu;
							}
							if (!px && !py && !pz) {
								for (i=0;i<K_elem;i++) {
									sc[order*(ia[i]-1)+ja[i]-1] += sa[i];
									if (ia[i] != ja[i]) sc[order*(ja[i]-1)+ia[i]-1]+=sa[i];
								}
							}
							else {
								dim=0;
								for (i=0;i<maxsize;i++) {
									ix=i/(Ny*Nz/4);
									iy=(i-ix*(Ny*Nz/4))/(Nz/2);
									iz=i-iy*(Nz/2)-ix*(Ny*Nz/4);
									if (radius(xx,xy,xz,ix,iy,iz) > rmax) continue;
									sc[order*dim+dim]=0.5*gamma*(1.0+alpha+beta
										-(xx[ix+1]*xx[ix+1]+alpha*xy[iy+1]*xy[iy+1]
										+beta*xz[iz+1]*xz[iz+1]));
									dum=dim;
									for (j=i;j<maxsize;j++) {
										jx=j/(Ny*Nz/4);
										jy=(j-jx*(Ny*Nz/4))/(Nz/2);
										jz=j-jy*(Nz/2)-jx*(Ny*Nz/4);
										if (radius(xx,xy,xz,jx,jy,jz) > rmax) continue;
										if (iy == jy && iz == jz) {
											sc[order*dim+dum]=0.0;
											for (k=px;k<Nx;k+=2) {
												sc[order*dim+dum] += 2.0*sqrt(weightx[ix+1])
													*sqrt(weightx[jx+1])*gamma*(double)k
													*hx[k+1][ix]*hx[k+1][jx];
												if (dim != dum) sc[order*dum+dim] += 2.0
													*sqrt(weightx[ix+1])*sqrt(weightx[jx+1])*gamma
													*(double)k*hx[k+1][ix]*hx[k+1][jx];
											}
										}
										if (ix == jx && iz == jz) {
											sc[order*dim+dum]=0.0;
											for (k=py;k<Ny;k+=2) {
												sc[order*dim+dum] += 2.0*alpha*sqrt(weighty[iy+1])
													*sqrt(weighty[jy+1])*gamma*(double)k
													*hy[k+1][iy]*hy[k+1][jy];
												if (dim != dum) sc[order*dum+dim] += 2.0*alpha
													*sqrt(weighty[iy+1])*sqrt(weighty[jy+1])*gamma
													*(double)k*hy[k+1][iy]*hy[k+1][jy];
											}
										}
										if (ix == jx && iy == jy) {
											sc[order*dim+dum]=0.0;
											for (k=pz;k<Nz;k+=2) {
												sc[order*dim+dum] += 2.0*beta*sqrt(weightz[iz+1])
													*sqrt(weightz[jz+1])*gamma*(double)k
													*hz[k+1][iz]*hz[k+1][jz];
												if (dim != dum) sc[order*dum+dim] += 2.0*beta
													*sqrt(weightz[iz+1])*sqrt(weightz[jz+1])*gamma
													*(double)k*hz[k+1][iz]*hz[k+1][jz];
											}
										}
										++dum;
									}
									++dim;
								}
							}
							for (i=0;i<order;i++) {
								for (j=0;j<order;j++) {
									for (k=0;k<order;k++) sb[order*i+j] += sc[order*i+k]
										*sc[order*k+j];
									sb[order*i+j] += sc[order*i+j]*(sd[j]+sd[i]
										+2.0*U*psi2[i]);
								}
							}

							for (i=0;i<V_elem;i++) {
								sb[order*(id[i]-1)+jd[i]-1]+=sd[i]*(sd[i]+2.0*U*psi2[i]);
							}
							printf("Bogoliubov matrix has been generated\n");

							/* dgeev(jobvl,jobvr,&order,sb,&lda,wr,wi,vl,&ldvl,
								vr,&ldvr,work,&lwork,&info); */
							sgeev(sb,&lda,&order,wr,vr,&ldvr,work,&info);

							printf("info = %d\n",info);

							eigs=cvector(0,3*(Emax*Emax-Emin*Emin));
							order_old=order;
							eigs_exact=fvector(0,order_old);

							e_index=0;
							for (i=0;i<order;i++) {
								eigs_exact[i]=wr[i].r;
								if (sqrt(wr[i].r) >= Emin && sqrt(wr[i].r) <= Emax) {
									printf("%d  %lf  %lf\n",i,sqrt(wr[i].r),wr[i].i);
									for (j=0;j<order;j++) {
										rvec[e_index][j]=vr[order*i+j].r;
										eigs[e_index]=wr[i];
									}
									++e_index;
								}
							}

							free_cvector(wr,0,order);
							free_cvector(vr,0,ldvr*order);
							free_fvector(work,0,2*order);
							free_ivector(id,0,V_elem);
							free_ivector(jd,0,V_elem);
							free_fvector(sd,0,V_elem);
							free_dvector(sc,0,order*order);
							free_dvector(sb,0,order*order);

						}  /* if statement checking if this is the coarse grid */

/* Second, use the interpolated excitation amplitudes as the filter for the
filter diagonalisation. This is accomplished by projecting the Hamiltonian into
the e_index-dimensional basis spanned by the orthonormal eigenfunctions, then
solving the resulting eigenproblem exactly for the eigenvalues. Note that
using the approximate excitations as preconditioners means that one doesn't
have to re-orthogonalise the space. */

						if (grid > 1 && grid == numgrid) {
							lda=e_index;
							ldvl=1;
							ldvr=1;
							lwork=10*e_index;
							ic=ivector(0,3*K_elem);
							jc=ivector(0,3*K_elem);
							sc=fvector(0,3*K_elem);
							id=ivector(0,V_elem);
							jd=ivector(0,V_elem);
							sd=fvector(0,V_elem);
							ib=ivector(0,3*K_elem+V_elem);
							jb=ivector(0,3*K_elem+V_elem);
							sb=fvector(0,3*K_elem+V_elem);
							diag=zvector(0,V_elem);

							if (!px && !py && !pz) {
								sc_len=0;
								for (i=0;i<K_elem;i++) {
									ic[sc_len]=ia[i];
									jc[sc_len]=ja[i];
									sc[sc_len]=sa[i];
									++sc_len;
									if (ia[i] != ja[i]) {
										ic[sc_len]=ja[i];
										jc[sc_len]=ia[i];
										sc[sc_len]=sa[i];
										++sc_len;
									}
								}
								printf("sc_len = %d, K_elem = %d\n",sc_len,K_elem);
							}
							else {
								dim=0;
								sc_len=0;
								for (i=0;i<maxsize;i++) {
									ix=i/(Ny*Nz/4);
									iy=(i-ix*(Ny*Nz/4))/(Nz/2);
									iz=i-iy*(Nz/2)-ix*(Ny*Nz/4);
									if (radius(xx,xy,xz,ix,iy,iz) > rmax) continue;
									sc[sc_len]=0.5*gamma*(1.0+alpha+beta-(xx[ix+1]*xx[ix+1]
										+alpha*xy[iy+1]*xy[iy+1]+beta*xz[iz+1]*xz[iz+1]));
									ic[sc_len]=dim+1;
									jc[sc_len]=dim+1;
									++sc_len;
									dum=dim;
									for (j=i;j<maxsize;j++) {
										jx=j/(Ny*Nz/4);
										jy=(j-jx*(Ny*Nz/4))/(Nz/2);
										jz=j-jy*(Nz/2)-jx*(Ny*Nz/4);
										if (radius(xx,xy,xz,jx,jy,jz) > rmax) continue;
										if (iy == jy && iz == jz) {
											sc[sc_len]=0.0;
											for (k=px;k<Nx;k+=2) {
												sc[sc_len] += 2.0*sqrt(weightx[ix+1])
													*sqrt(weightx[jx+1])*gamma*(double)k
													*hx[k+1][ix]*hx[k+1][jx];
											}
											ic[sc_len]=dim+1;
											jc[sc_len]=dum+1;
											++sc_len;
											if (ic[sc_len-1] != jc[sc_len-1]) {
												ic[sc_len]=dum+1;
												jc[sc_len]=dim+1;
												sc[sc_len]=sc[sc_len-1];
												++sc_len;
											}
										}
										if (ix == jx && iz == jz) {
											sc[sc_len]=0.0;
											for (k=py;k<Ny;k+=2) {
												sc[sc_len] += 2.0*alpha*sqrt(weighty[iy+1])
													*sqrt(weighty[jy+1])*gamma*(double)k
													*hy[k+1][iy]*hy[k+1][jy];
											}
											ic[sc_len]=dim+1;
											jc[sc_len]=dum+1;
											++sc_len;
											if (ic[sc_len-1] != jc[sc_len-1]) {
												ic[sc_len]=dum+1;
												jc[sc_len]=dim+1;
												sc[sc_len]=sc[sc_len-1];
												++sc_len;
											}
										}
										if (ix == jx && iy == jy) {
											sa[K_elem]=0.0;
											for (k=pz;k<Nz;k+=2) {
												sc[sc_len] += 2.0*beta*sqrt(weightz[iz+1])
													*sqrt(weightz[jz+1])*gamma*(double)k
													*hz[k+1][iz]*hz[k+1][jz];
											}
											ic[sc_len]=dim+1;
											jc[sc_len]=dum+1;
											++sc_len;
											if (ic[sc_len-1] != jc[sc_len-1]) {
												ic[sc_len]=dum+1;
												jc[sc_len]=dim+1;
												sc[sc_len]=sc[sc_len-1];
												++sc_len;
											}
										}
										++dum;
									}
									++dim;
								}
							}
							for (i=0;i<V_elem;i++) {
								id[i]=ia[K_elem+i];
								jd[i]=ja[K_elem+i];
								sd[i]=sa[K_elem+i]-mu;
							}
							sb_len=0;
							for (i=0;i<sc_len;i++) {  /* TV+(V+2V_H)T */
								ib[sb_len]=ic[i];
								jb[sb_len]=jc[i];
								sb[sb_len]=sc[i]*(sd[ic[i]-1]+sd[jc[i]-1]+2.0*U*psi2[jc[i]-1]);
								/* if (ic[i] == jc[i]) diag[ic[i]-1].r += sc[i]; */
								/* if (ib[sb_len] == jb[sb_len]) 
									diag[ib[sb_len]-1].r += sb[sb_len]; */
								++sb_len;
							}
							for (i=0;i<V_elem;i++) {  /* add V*(V+2V_H) */
								ib[sb_len]=id[i];
								jb[sb_len]=jd[i];
								sb[sb_len]=sd[i]*(sd[i]+2.0*U*psi2[i]);
								/* if (ib[sb_len] == jb[sb_len])
									diag[ib[sb_len]-1].r += sb[sb_len]; */
								++sb_len;
							}
							printf("sb_len = %d\n",sb_len);

							/* j=0;
							for (i=0;i<V_elem;i++) {
								diag[i].r=eigs_exact[j];
								if (j<order_old-1) ++j;
								else j=0;
							} */
							for (i=0;i<V_elem;i++) printf("%d  %lf\n",i,diag[i]);

							maxit=100;
							drctv=FALSE_;  /* preconditioner flag */
							nonsym=TRUE_;  /* true if non-symmetric */
							toler=TOL;
							n_eigs=4;  /* number of excitations to calculate */
							scale=1.0;
							cdvd(vex,eigs,&e_index,&toler,&toler,&dim,&n_eigs,&maxit,
								&drctv,&nonsym,&sb_len,ib,jb,sb,&sc_len,ic,jc,sc,diag);
							for (i=0;i<e_index;i++) printf("%lf  %lf\n",eigs[i].r,eigs[i].i);

							free_cvector(eigs,0,3*(Emax*Emax-Emin*Emin));
							free_cvector(vex,0,e_index*dim);
							free_zvector(diag,0,V_elem);

						} /* if statement checking whether we're all done */

						if (n == imax) exit(1);
						break;
					}
					if (n == imax && grid != numgrid) {
						for (i=0;i<order;i++) w[i]=w_good[i];
						bi_low[grid+1]=bi;
						bi=bi_high[grid]+1;
						break;
					}
					if (n == imax && grid == numgrid) exit(1);
				} /* iteration loop for self-consistency */

				free_fmatrix(B,1,imax+1,1,imax+1); /* arrays for the DIIS routine */
				free_fvector(out,0,imax);

				if (bi == bi_max && grid == numgrid) {
					free_ivector(ia,0,entries);
					free_ivector(ja,0,entries);
					free_fvector(sa,0,entries);
					free_fmatrix(p,0,imax+1,0,dim_alloc);
					free_fmatrix(pa,0,imax+1,0,dim_alloc);
					free_fvector(psi,0,dim_alloc);
					free_fvector(psi2,0,dim_alloc);
					free_fvector(mut,0,imax+1);
					free_fvector(xx,1,Nx);
					free_fvector(xy,1,Ny);
					free_fvector(xz,1,Nz);
					free_fvector(weightx,1,Nx);
					free_fvector(weighty,1,Ny);
					free_fvector(weightz,1,Nz);
					free_fvector(w,0,Q*dim_alloc);
					free_fvector(w_good,0,Q*dim_alloc);
					free_fmatrix(hx,0,Nx,0,Nx);
					free_fmatrix(hy,0,Ny,0,Ny);
					free_fmatrix(hz,0,Nz,0,Nz);
				}
			}  /* kinetic contribution loop */
		}  /* grid loop */
	}  /* Np loop */
	fclose(fout);
}

int *ivector(long nl, long nh)
{
  int *v;

  v=(int *)malloc((size_t) ((nh-nl+1+NR_END)*sizeof(int)));
  if (!v) { printf("allocation failure in ivector\n"); exit(1); }
  return v-nl+NR_END;
}

double *fvector(long nl, long nh)
{
  double *v;

  v=(double *)malloc((size_t) ((nh-nl+1+NR_END)*sizeof(double)));
  if (!v) { printf("allocation failure in fvector\n"); exit(1); }
  return v-nl+NR_END;
}

double *dvector(long nl, long nh)
{
  double *v;

  v=(double *)calloc((size_t)1, (size_t) ((nh-nl+1+NR_END)*sizeof(double)));
  if (!v) { printf("allocation failure in dvector\n"); exit(1); }
  return v-nl+NR_END;
}

doublecomplex *cvector(long nl, long nh)
{
  doublecomplex *v;

  v=(doublecomplex *)malloc((size_t) ((nh-nl+1+NR_END)
		*sizeof(doublecomplex)));
  if (!v) { printf("allocation failure in cvector\n"); exit(1); }
  return v-nl+NR_END;
}

doublecomplex *zvector(long nl, long nh)
{
  doublecomplex *v;

  v=(doublecomplex *)calloc((size_t)1, (size_t) ((nh-nl+1+NR_END)
		*sizeof(doublecomplex)));
  if (!v) { printf("allocation failure in zvector\n"); exit(1); }
  return v-nl+NR_END;
}

double **fmatrix(long nrl, long nrh, long ncl, long nch)
/* allocate a double matrix with subscript range m[nrl..nrh][ncl..nch] */
{
	long i, nrow=nrh-nrl+1,ncol=nch-ncl+1;
	double **m;

	/* allocate pointers to rows */
	m=(double **) malloc((size_t)((nrow+NR_END)*sizeof(double*)));
	if (!m) {
		printf("allocation failure 1 in matrix()");
		exit(1);
	}
	m += NR_END;
	m -= nrl;

	/* allocate rows and set pointers to them */
	m[nrl]=(double *) malloc((size_t)((nrow*ncol+NR_END)*sizeof(double)));
	if (!m[nrl]) {
		printf("allocation failure 2 in matrix()");
		exit(1);
	}
	m[nrl] += NR_END;
	m[nrl] -= ncl;

	for(i=nrl+1;i<=nrh;i++) m[i]=m[i-1]+ncol;

	/* return pointer to array of pointers to rows */
	return m;
}

double **hermit(double *x, long N)
{
	int i,j;
	double **ht,p1,p2,p3;
	
	ht=fmatrix(0,N,0,N);

	for (j=0;j<N;j++) {
		ht[0][j]=0.0;
		ht[1][j]=PIM4;
	}
	for (i=0;i<N-1;i++) {
		for (j=0;j<N;j++) {
			ht[i+2][j]=x[j+1]*sqrt((double)(2.0/((double)i+1.0)))*ht[i+1][j]
				-sqrt((double)i/((double)(i+1)))*ht[i][j];
		}
	}
	return(ht);
}

void free_ivector(int *v, long nl, long nh)
{
	free((FREE_ARG) (v+nl-NR_END));
}

void free_fvector(double *v, long nl, long nh)
{
	free((FREE_ARG) (v+nl-NR_END));
}

void free_dvector(double *v, long nl, long nh)
{
	free((FREE_ARG) (v+nl-NR_END));
}

void free_cvector(doublecomplex *v, long nl, long nh)
{
	free((FREE_ARG) (v+nl-NR_END));
}

void free_zvector(doublecomplex *v, long nl, long nh)
{
	free((FREE_ARG) (v+nl-NR_END));
}

void free_fmatrix(double **m, long nrl, long nrh, long ncl, long nch)
{
	free((FREE_ARG) (m[nrl]+ncl-NR_END));
	free((FREE_ARG) (m+nrl-NR_END));
}

double radius(double *xp, double *yp, double *zp, long a, long b, long c)
{
	double ans;

	ans=sqrt(xp[a+1]*xp[a+1]+yp[b+1]*yp[b+1]+zp[c+1]*zp[c+1]);
	return(ans);
}

void print_psi(double *p, double *x, double *y, double *z, double r, int flag,
	int Nx, int Ny, int Nz, int size, double gamma)
{
	int i,ix,iy,iz,dim;

	dim=0;
	for (i=0;i<size;i++) {
		ix=i/(Ny*Nz/4);
		iy=(i-ix*(Ny*Nz/4))/(Nz/2);
		iz=i-iy*(Nz/2)-ix*(Ny*Nz/4);
		if (radius(x,y,z,ix,iy,iz) > r) continue;
		if (iy == Ny/2-1 && iz == Nz/2-1) printf("%lf  %lf  %lf\n",
			x[ix+1]/sqrt(gamma),p[dim],p[dim]*p[dim]);
		++dim;
	}
	dim=0;
	for (i=0;i<size;i++) {
		ix=i/(Ny*Nz/4);
		iy=(i-ix*(Ny*Nz/4))/(Nz/2);
		iz=i-iy*(Nz/2)-ix*(Ny*Nz/4);
		if (radius(x,y,z,ix,iy,iz) > r) continue;
		if (ix == Nx/2-1 && iz == Nz/2-1) {
			if (!flag) printf("%lf  %lf  %lf\n",y[iy+1]/sqrt(gamma),p[dim],
				p[dim]*p[dim]);
			else printf("%lf  %lf  %lf\n",y[iy+1]/sqrt(alpha*gamma),p[dim],
				p[dim]*p[dim]);
		}
		++dim;
	}
	dim=0;
	for (i=0;i<size;i++) {
		ix=i/(Ny*Nz/4);
		iy=(i-ix*(Ny*Nz/4))/(Nz/2);
		iz=i-iy*(Nz/2)-ix*(Ny*Nz/4);
		if (radius(x,y,z,ix,iy,iz) > r) continue;
		if (ix == Nx/2-1 && iy == Ny/2-1) {
			if (!flag) printf("%lf  %lf  %lf\n",z[iz+1]/sqrt(gamma),p[dim],
				p[dim]*p[dim]);
			else printf("%lf  %lf  %lf\n",z[iz+1]/sqrt(beta*gamma),p[dim],
				p[dim]*p[dim]);
		}
		++dim;
	}
}

/* next is the doublecomplex routine to multiply H on a vector. Note that H
is real */

void honv2(int *n, int *vlen, int *iv, int *jv, double *sv, int *tlen, 
	int *it, int *jt, double *st, doublecomplex *vec, doublecomplex *hvec,
	int *num1, int *num2)
{
	int i,j;
	doublecomplex *tt, *uu, *vv;

/* for (j=0;j<(*num1)+(*num2)-1;j++) { */
	for (j=(*num1)-1;j<(*num2);j++) {

/* first generate T*T*x */

		tt=zvector(0,*n);
		uu=zvector(0,*n);
		vv=zvector(0,*n);

		for (i=0;i<*tlen;i++) {
			uu[jt[i]-1].r += st[i]*vec[j*(*n)+it[i]-1].r;
			uu[jt[i]-1].i += st[i]*vec[j*(*n)+it[i]-1].i;
		}
		for (i=0;i<*tlen;i++) {
			tt[jt[i]-1].r += st[i]*uu[it[i]-1].r;
			tt[jt[i]-1].i += st[i]*uu[it[i]-1].i;
		}

	/* now do the rest */
			
		for (i=0;i<*vlen;i++) {
			vv[jv[i]-1].r += sv[i]*vec[j*(*n)+iv[i]-1].r;
			vv[jv[i]-1].i += sv[i]*vec[j*(*n)+iv[i]-1].i;
		}
		for (i=0;i<*n;i++) {
			hvec[j*(*n)+i].r=tt[i].r+vv[i].r;
			hvec[j*(*n)+i].i=tt[i].i+vv[i].i;
		}

		free_zvector(tt,0,*n);
		free_zvector(uu,0,*n);
		free_zvector(vv,0,*n);
	}
}
