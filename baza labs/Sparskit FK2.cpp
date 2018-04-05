#include <stdio.h>
#include <string.h>
#include <math.h>
#include <malloc.h>
#include <ctime>


#if defined(HAVE_UNDERSCORE)
#define etime_      etime
#define ilut_       ilut
#define ilutp_      ilutp
#define ilud_       ilud
#define iludp_      iludp
#define iluk_       iluk
#define ilu0_       ilu0
#define milu0_      milu0
#define cmkreord_   cmkreord
#define sortcol_    sortcol
#define skit_       skit
#define psplot_     psplot
#define cg_         cg
#define cgnr_       cgnr
#define bcg_        bcg
#define dbcg_       dbcg
#define bcgstab_    bcgstab
#define tfqmr_      tfqmr
#define fom_        fom
#define gmres_      gmres
#define fgmres_     fgmres
#define dqgmres_    dqgmres
#define amux_       amux
#define atmux_      atmux
#define lusol_      lusol
#define lutsol_     lutsol
#define csrcoo_     csrcoo
#define ma28ad_     ma28ad
#define ma28cd_     ma28cd
#define dnrm2_      dnrm2
#define flu_        flu
#define pgmres_     pgmres
#define getdia_     getdia
#define amudia_     amudia
#define diamua_     diamua
#define rnrms_      rnrms
#define cnrms_      cnrms
#endif

/* Fortran prototypes */
extern "C" {
  void  ilut_     (int*,double*,int*,int*,int*,double*,
		   double*,int*,int*,int*,double*,int*,int*);
  void  ilutp_    (int*,double*,int*,int*,int*,double*,
		   double*,int*,double*,int*,
		   int*,int*,double*,int*,int*,int*);
  void  ilud_     (int*,double*,int*,int*,double*,
		   double*,double*,int*,
		   int*,int*,double*,int*,int*);
  void  iludp_    (int*,double*,int*,int*,double*,
		   double*,double*,
		   int*,double*,int*,int*,int*,
		   double*,int*,int*,int*);
  void  iluk_     (int*,double*,int*,int*,int*,
		   double*,int*,int*,
		   int*,int*,double*,int*,int*);
  void  ilu0_     (int*,double*,int*,int*,double*,int*,int*,int*,int*);
  void  milu0_    (int*,double*,int*,int*,double*,int*,int*,int*,int*);
  void  cmkreord_ (int*,double*,int*,int*,double*,int*,int*,int*,
		   int*,int*,int*,int*,int*,int*);
  void  sortcol_  (int*,double*,int*,int*,int*,double*);
  void  skit_     (int*,double*,int*,int*,int*,int*,int*);
  void  psplot_   (int*,int*,int*,int*,int*);
  void  cg_       (int*,double*,double*,int*,double*,double*);
  void  cgnr_     (int*,double*,double*,int*,double*,double*);
  void  bcg_      (int*,double*,double*,int*,double*,double*);
  void  dbcg_     (int*,double*,double*,int*,double*,double*);
  void  bcgstab_  (int*,double*,double*,int*,double*,double*);
  void  tfqmr_    (int*,double*,double*,int*,double*,double*);
  void  fom_      (int*,double*,double*,int*,double*,double*);
  void  gmres_    (int*,double*,double*,int*,double*,double*);
  void  fgmres_   (int*,double*,double*,int*,double*,double*);
  void  dqgmres_  (int*,double*,double*,int*,double*,double*);
  void  amux_     (int*,double*,double*,double*,int*,int*);
  void  atmux_    (int*,double*,double*,double*,int*,int*);
  void  lusol_    (int*,double*,double*,double*,int*,int*);
  void  lutsol_   (int*,double*,double*,double*,int*,int*);
  void  csrcoo_   (int*,int*,int*,double*,int*,int*,int*,double*,int*,int*,int*);
  void  ma28ad_   (int*,int*,double*,int*,int*,int*,int*,double*,int*,int*,double*,int*);
  void  ma28cd_   (int*,double*,int*,int*,int*,double*,double*,int*);
  double dnrm2_   (int*,double*,int*);
  void  flu_      (int*,double*,double*,double*,int*,int*,double*,double*,
		   double*,double*,double*);
  void  pgmres_   (int*,int*,double*,double*,double*,double*,
		   int*,int*,double*,int*,int*,
		   double*,int*,int*,int*);
  void getdia_    (int*,int*,int*,double*,int*,int*,int*,double*,int*,int*);
  void diamua_    (int*,int*,double*,int*,int*,double*,double*,int*,int*);
  void amudia_    (int*,int*,double*,int*,int*,double*,double*,int*,int*);
  void rnrms_     (int*,int*,double*,int*,int*,double*);
  void cnrms_     (int*,int*,double*,int*,int*,double*);
}



int main() {
  double   fpar[17];
  double  *a, *alu, *w, *rhs, *sol;
  double   res;
  double   res1=1.;
  int      i, j, k, iwk, ierr, ipar[17];
  int     *ja, *ia, *jw, *jlu, *ju;
  int      its, end;
  int lfil;
  double droptol;
  int m=0;
  int n, nnz;


  /* parameters for changing*/
	char f_name[6] = "4.dat";
  lfil = 80;
  droptol = 1.0/100;


/*cicle po imeni*////////////////////////
/*
int v=0;
f_name[0] = "1.dat";
f_name[1] = "2.dat";
f_name[2] = "3.dat";
f_name[3] = "4.dat";
for (v = 0; v < 4;v++){
*/

  /*Reading Matrix */
  FILE *file;
	i=0;
	float ss;
	int s;
	file=fopen(f_name,"r+");
	fscanf(file, "%d", &n);
	ia = (int*)calloc(n+1,sizeof(int));
	for (i = 0; i < n+1;i++) {
		fscanf(file, "%d",&s);
		ia[i] = s;
	}
	nnz = ia[n]-1;
	ja = (int*)calloc(nnz,sizeof(int));
	a = (double*)calloc(nnz,sizeof(double));
	for (i = 0; i < nnz;i++) {
		fscanf(file, "%d",&s);
		ja[i] = s;
	}
	for (i = 0; i < nnz;i++) {
		fscanf(file, "%f",&ss);
		a[i] = ss;
	}
	fclose(file); 
	

  /*RHS*/
  rhs = (double*) malloc(n * sizeof(double));
  for ( size_t i = 0; i < n; i++) {		
    rhs[i] = sin(double(i));
  }
	



  sol = (double*) calloc(n, sizeof(double));
  end = 0;

	iwk = 2 * (n+1) * (lfil+1); 
  alu = (double*) malloc(iwk * sizeof(double));
  jlu = (int*) malloc(iwk * sizeof(int));
  ju  = (int*) malloc((n+1) * sizeof(int));
    
  w  = (double*) malloc((n+1) * sizeof(double));
  jw = (int*) malloc(2 * (n+1) * sizeof(int));


	/*ILUT */
	clock_t start_t = clock();	
  ilut_(&n, a, ja, ia, &lfil, &droptol, alu, jlu, ju, &iwk, w, jw, &ierr); 
	clock_t end_t = clock();
	double time_ilut = double(end_t - start_t)/ double(CLOCKS_PER_SEC);
  free(w); free(jw); 

  /* Iterations */

  ipar[1] = 0;
  ipar[2] = 2;
  ipar[3] = 1;
  ipar[4] = 7 * n;
  ipar[5] = 10;
  ipar[6] = 500;

  fpar[1] = 1.0E-6;
  fpar[2] = 1.0E-10;
  fpar[11] = 0.0;

  w = (double*) malloc(ipar[4] * sizeof(double));

  its = 0;
  end = 0;
  res = 0.0;

	start_t = clock();
	int iter_number = 0;
  while(1) {
  	bcg_(&n, rhs, sol, &ipar[1], &fpar[1], w); 
    if(!end){
      if(ipar[7] != its){
	//if(its) printf(" %4d  %.7e  %.7e", its, res, res/res1);
	its = ipar[7] ;
      }
      res = fpar[5];
      if(its==1) res1 = fpar[5] ;

      switch(ipar[1]){
      case 1 :
	amux_(&n, &w[ipar[8]-1], &w[ipar[9]-1], a, ja, ia); break;
      case 2 :
	atmux_(&n, &w[ipar[8]-1], &w[ipar[9]-1], a, ja, ia); break;
      case 3 : case 5 :
	lusol_(&n, &w[ipar[8]-1], &w[ipar[9]-1], alu, jlu, ju); break;
      case 4 : case 6 :
	lutsol_(&n, &w[ipar[8]-1], &w[ipar[9]-1], alu, jlu, ju); break;
      case 0 :
	end = 1; break;
      case -1 :
	printf("Iterative solver has iterated too many times\n"); end = 1; break;
      case -2 :
	printf("Iterative solver was not given enough work space\n");
	printf("The work space should at least have %d elements\n", ipar[4]);
	end = 1; break;
      case -3 :
	printf("Iterative solver is facing a break-down\n"); 
	end = 1; break;
      default :
	printf("Iterative solver terminated (code = %d)\n", ipar[1]); 
	end = 1; break;
      }

    }
    if(end) break;
	iter_number++;
  }
	end_t = clock();
	double time_iter = double(end_t - start_t)/ double(CLOCKS_PER_SEC);

	printf("File %s:\n", f_name);
	printf ("for lfil = %d, droptol = %.3f => ILUT time = %.6f, iteration time = %.6f, iteration number = %d \n", lfil, droptol, time_ilut, time_iter, iter_number);

	/*FREE */
	free(a); free(ia); free(ja); free(rhs); free(alu); free(jlu); free(w); free(ju);

  return 0;
  }


