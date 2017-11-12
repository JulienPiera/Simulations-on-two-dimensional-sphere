//************* version mars 2016

typedef struct
{
  int N; //nombre de particules
  int N2; //nombre de particules de type 2 (si polydisperse)
  int nthreads; // nombre de threads utilisés
  long bstep,nstep; // step before and after
  int moviestep; // pas entre enregistrement de 2 images
  int tFs; //temps sur lequel on calcule Fs
  double eta; //packing fraction
  double eta2; //second packing fraction for a polydisperse system (eta2=0 if monodisperse)
  double sigma; //diameter
  double sigma2; //second diameter
  double sigma12; // diametre associé au potentiel entre 1 particule 1 et 1 particule 2
  double temp; //temperature
  double beta; //inverse temperature
  double dt; //pas de temps de la dynamique molÃ©culaire
  double deltaphi; //amplitude delta phi
  double deltaz; //amplitude delta z
  double nrjbas; // lower bound of the energy histogram
  double nrjhaut; // upper bound of the energy histogram
  double deltar;
  double kFs; // vecteur d'onde choisi pour Fs(k,t), en unitÃ©s de pi/sigma
  int kmax; // nombre d'onde max dans le calcul de S(k), entier
  int kcut; // nombre d'onde max où on calcule c(k) pour le F_infini
  int dtmoy; // pas dans le calcul des valeurs moyennes temporelles
  int nconf; // nombre de configs sur lesquelles on calcule S(k)
  int maxpot;
  int maxhisto;
  int maxacos;  
  int maxlegendre;
  int flagposanc; // vaut 0 si on part d'une config alÃ©atoire, vaut 1 si on part de la derniÃ¨re config (positions) d'une ancienne simu
  int flagtherma; // 0 si on calcule l'histogramme depuis le dÃ©but, 1 si on reprend l'histo d'une ancienne simu
  int flagcorr; // 1 si on calcule g(r), 0 sinon
 // int flagcorrdefauts; // 1 si on calcule g_q(r)
  int flagstruc; // 1 si on calcule S(k), 0 sinon
  int flagg6; // 1 si on calcule G6(r), 0 sinon
  int flagFs; // 1 si on calcule Fs, 0 sinon
  int flagFtot; // 1 si on calcule le facteur de structure complet, 0 sinon
  int flagkhi4; // 1 si on calcule khi 4, 0 sinon  
  int flagdeplac; //1 si on calcule <r(t)>, 0 sinon
  int flagdeplac657; //1 si on calcule le déplacement moyen des particules 6
  int flagdevoro; //1 si on calcule les defaults, 0 sinon
  int flagdynam; //1 si un fichier de configuration existe, 0 sinon
  int flagkwave; // 0 si k=2*pi/sigma , 1 si k=2*pi/R
  int flagmovie;
  int flagsavedefauts; // 1 si on enregistre les positions des défauts
  //int flagpara; //1 si on fait du calcul parallèle, 0 sinon
  int flagspectre;
  int bstep_old; // nombre de pas de temps de l'ancienne simu qu'on reprend, avant thermalisation
  int step_old;
} RUN;


typedef struct
{
  int flagr; //=0 avant thermalisation, 1 aprÃ¨s
  double nrj;
  double nrj2;
  double nrjac; //energie accumulÃ©e jusqu'Ã  t
  double nrjac2;
  double cv; //chaleur latente, calculÃ©e Ã  partir des grandeurs ci dessus
  double *histoener;
  double *histovoro;
  double *ndislo; // nbre moyen de dislocations 5-7
  double *n66; // nbre de particules 6 dont les voisins sont 6
} THERMO;


typedef struct
{
 double phi;
 double z;
} POS;



typedef struct
{
 double x;
 double y;
 double z;
 double col;
} POS_cart;


typedef struct
{
 double theta;
 double phi;
} POS_ang;

typedef struct
{
 double Re;
 double Im;
} complexe; 


static const double pi = M_PI;
static const double pi2=2.*M_PI;
// function prototypes
void read_data(int argc, char *argv[], RUN *run);
//void read_data_structure(int argc, char *argv[], double *struc_fact, RUN *run);
void read_structure(double *struc_fact, RUN run);
void save_defauts(int flag, int step, POS_cart *pos2_cart, THERMO *thermo, RUN run);
void dyn_molec_init(POS_cart *pos1_cart, POS_cart *omg1_cart, POS_ang *pos1, POS_ang *vit1, double *arccos, RUN run, gsl_rng *rng);
inline double poten(double x,double *pot,RUN run);
inline double terme_A(double x,double *force_A,RUN run);
inline double distance(POS posi,POS posj, double *arccos, RUN run);
inline double distance_ang(POS_ang posi,POS_ang posj, double *arccos, RUN run);
inline double distance_cart(POS_cart posi,POS_cart posj, double *arccos, RUN run);

void reinit_vit(POS_cart *pos1_cart, POS_cart *omg1_cart,  POS_ang *pos1, POS_ang *vit1, double *arccos, RUN run, gsl_rng *rng);

inline double vitesse(int i, POS_ang *pos1, POS_ang *pos2, POS_ang *pos3, double *arccos, RUN run);
//inline double arccos(double x, double *arccos, RUN run);
inline double arccosinus(double x, double *arccos, RUN run);
inline double calcul_legendre(double x, double *legendre, RUN run);

inline void histo_distances(int i1, double *hist_gr, POS_cart *pos, double *arccos, RUN run); // ajoutÃ©
inline void histo_deplacements(int i0, int iter, double **hist_distances, POS *pos, double *arccos, RUN run);


void calcul_struc_inf(double *fact1, double *fact2, double *struc_fact, double *corr_directe, double *memoire1, const RUN run, gsl_rng *rng);
void calcul_distrib_paire(int flag, int step,  POS_cart *pos2_cart, double *arccos, RUN run, gsl_rng *rng); // calcule g(r)
void calcul_corr_defauts(int flag, int step, POS_cart *pos2_cart, double *arccos, RUN run, gsl_rng *rng);
void spectre_defauts(int flag, int step, POS_cart *pos2_cart, THERMO *thermo, double *arccos, RUN run);
void calcul_nrj_pot(int flag, int step, POS_cart *pos2_cart, double *pot, double *pot2, double *pot12, THERMO *thermo, double *arccos, RUN run);
void calcul_nrj_pot_lissee(int flag, int step, POS_cart *pos2_cart, double *pot, double *pot2, double *pot12, THERMO *thermo, double *arccos, RUN run);
void calcul_images(int flag, int step, POS_cart *pos2_cart, THERMO *thermo, RUN run);
void calcul_deplac(int flag, int step, int dt_Fs, POS_cart *pos2_cart, double *arccos, RUN run);
void calcul_deplac2(int flag, int step, int dt_Fs, POS_cart *pos2_cart, double *arccos, RUN run);
void calcul_deplac657(int flag, int step, int dt_Fs, POS_cart *pos2_cart, double *arccos, THERMO *thermo, RUN run);
void calcul_dist20part(int flag, int step, POS_cart *pos2_cart, double *arccos, RUN run);
void calcul_Fs(int flag, int step, int dt_Fs, POS_cart *pos2_cart, double *Fs, double *arccos, double *legendre, RUN run);
void calcul_Ftot(int flag, int step, int dt_Fs, POS_cart *pos2_cart, double *arccos, double *legendre, RUN run);
void calcul_khi4(int flag, int step, POS_cart *pos2_cart, double *Fs, double *arccos, double *legendre, RUN run);
void calcul_scar(int flag, double *arccos, RUN run, gsl_rng *rng);
void calcul_voisins_defauts(int i, int ndefauts, int *nscar, int *scar_list, POS_cart *pos_defauts, double *arccos, int *voisins, int *col, RUN run, gsl_rng *rng);
POS_cart calcul_pos_scar(int step,int ndefauts, int *nscar, 
		     int *scar_list, POS_cart *pos_defauts, double *arccos, int *voisins, int *col, RUN run, gsl_rng *rng);
void tourne(POS_cart *a, POS_cart *b, POS_cart *c, POS_cart *mem_a, POS_cart *mem_b, POS_cart *mem_c, RUN run);

void initial_pot(double *pot, double *pot2, double *pot12, RUN run );
void initial_force(double *force_A,RUN run);
void initial_config(POS *pos,RUN run,gsl_rng *rng);
void initial_config_anc(POS_cart *pos, THERMO *thermo, RUN run); //rÃ©initialiser Ã  partir d'une autre simu
void arccos_tabu(double *arccos,RUN run); // tabule arccosinus
void arccos_interpol(double *arccos, double *interpol_A, double *interpol_B,RUN run);

void legendre_tabu(double *legendre,RUN run); // tabule polynÃ´me de legendre d'ordre (int)(2*pi/sigma)


inline void force(POS_ang *terme_force, double *force_A, POS_ang *pos2, double *arccos, RUN run);

inline void verlet(POS_cart *pos1_cart, POS_cart *pos2_cart, POS_cart *omg1_cart, 
	     POS_cart *omg2_cart, POS_cart *force1_cart, POS_cart *force2_cart,
		  double *arccos, double *force_A, RUN run);
inline POS_cart somme_vect(POS_cart a, POS_cart b);
inline double prod_scalaire(POS_cart a, POS_cart b);
inline POS_cart prod_vect(POS_cart a, POS_cart b);
inline POS_cart scalaire(double lambda, POS_cart a);
inline complexe somme_complexe(complexe a, complexe b);
inline complexe prod_complexe(complexe a, complexe b);
inline complexe conjug_complexe(complexe a);
inline double module_complexe(complexe a);
inline complexe exp_complexe(complexe a);
inline complexe scalaire_complexe(double lambda, complexe a);
int line_count(FILE *n);

inline POS_cart force_cart(int i, POS_cart *pos_cart,double *arccos, double *force_A, RUN run);


void dyn_molec(THERMO *thermo, POS_cart *pos1_cart, POS_cart *pos2_cart, POS_cart *omg1_cart, POS_cart *omg2_cart, double *pot, double *pot2, double *pot12,
	       double *force_A,double *arccos,double *legendre, RUN run, gsl_rng *rng);

void write_data(POS *pos,THERMO *thermo,RUN run);


void inversion_systeme(POS_cart *omg, double a, double b, double c, double d, double e, double f, double omgx, double omgy, double omgz);


double calcul_nrj(POS *pos, double *pot,double *arccos,RUN run);
double calcul_nrj_ang(POS_ang *pos2, double *pot, double *arccos, RUN run);
// inline double calcul_nrj_cart(POS_cart *pos1_cart,double *pot, double *arccos, double **tab_distances, RUN run);
inline double calcul_nrj_cart(POS_cart *pos1_cart, double *pot, double *pot2, double *pot12, double *arccos, RUN run);
void voronoi(double *x, double *y, double *z, int *voisins, int *col, int n, int *ndislo);
double interpol_lin(double a, double x, double fx, double y, double fy);


void ELAPSEDTIME(const char *s,double duree);
void allocdouble(double **array,int  xsize);
//void allocdouble2(int ***array,int  xsize,int ysize);
void allocint(int **array,int  xsize);

// function definitions
void allocint(int **array,int  xsize)
/* dynamically allocates a one dimensional (*array)[xize]
   with datasize bytes per element */
{
 if (*array != NULL)  /* free previous allocation */
   free(*array);
  /* allocate memory for array of xsize pointers plus data */
 *array=(int *) calloc(xsize, sizeof(int));
 if (*array ==  NULL)
  {fprintf(stderr," memory allocation problems\n");
   exit(-1);}
 return;
}



  void allocdouble(double **array,int  xsize)
/* dynamically allocates a one dimensional (*array)[xize]
   with datasize bytes per element */
{
  if (*array != NULL)  /* free previous allocation */
   free(*array);
  /* allocate memory for array of xsize pointers plus data */
 *array=(double *) calloc(xsize, sizeof(double));
 if (*array ==  NULL)
  {fprintf(stderr," memory allocation problems with allocdouble\n");
   exit(-1);}
 return;
}
/*-----------------------------------------------------------*/
void ELAPSEDTIME(const char *s,double duree)
{
  double fduree_s , duree_mn_temp;
  int duree_h,duree_mn,duree_s;
  fduree_s      = fmod(duree,60.);
  duree_s      = floor(fduree_s);
  duree_mn_temp = ( duree - fduree_s ) / 60.;
  duree_mn     = floor( fmod(duree_mn_temp,60.) );
  duree_h      = ( duree_mn_temp - duree_mn ) / 60.;
  printf("%s %dh  %dmn  %fs\n",s,duree_h,duree_mn,fduree_s);
  return;
}


 void allocdouble2(double ***array,int  xsize,int ysize)
 /* dynamically allocates a two dimensional (*array)[xize][ysize]
    with datasize bytes per element */
 { int i;
  if (*array != NULL)  /* free previous allocation */
    free(*array);
   /* allocate memory for array of ysize pointers plus data */
  *array=(double **) malloc(xsize*sizeof(double *));
  if (*array ==  NULL)
   {
     fprintf(stderr," memory allocation problems with allocdouble2\n");
    exit(-1);
  }
  for(i=0;i<xsize;i++)
  {
  (*array)[i]=(double *) calloc(ysize,sizeof(double));
   }
  return;
 }
 
 void allocint2(int ***array,int  xsize,int ysize)
 /* dynamically allocates a two dimensional (*array)[xize][ysize]
    with datasize bytes per element */
 { int i;
  if (*array != NULL)  /* free previous allocation */
    free(*array);
   /* allocate memory for array of ysize pointers plus data */
  *array=(int **) malloc(xsize*sizeof(int *));
  if (*array ==  NULL)
   {
     fprintf(stderr," memory allocation problems with allocint2\n");
    exit(-1);
  }
  for(i=0;i<xsize;i++)
  {
  (*array)[i]=(int *) calloc(ysize,sizeof(int));
   }
  return;
 }


//********************************************************************************
//*******************************************************************************


           
