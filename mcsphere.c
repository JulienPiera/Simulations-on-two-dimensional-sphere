// on fixe la fraction d'empilement et le nombre de particles
// on calcule le diametre d'une particule
#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <math.h>
#include <string.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_sf.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_minmax.h>
#include <stdbool.h>
#include <omp.h>
#include "mcsphere.h"
#include "initialisation.h"
#include "math.h"
#include "energie.h"
#include "fonctions.h"
//#include "fonctions_par.h"

// 
// parallelisation du Fs et du deplacement


//double prod_scal_max;
//double prod_scal_max2;
//double prod_scal_max12;
double min_acos; // valeur minimale dans tabulation de arccos : min_acos=cos(2.5*sigma)
double max_acos; // valeur maximale dans tabulation de arccos : max_acos=cos(0.5*sigma)
//const double epsilon=1.0;




int main (int argc, char *argv[]){

 const gsl_rng_type *rngtype=gsl_rng_mt19937;

   
   gsl_rng *rng;

   rng=gsl_rng_alloc(rngtype);
   gsl_rng_env_setup();
   //test
   /*
   int l;
   for(l=0;l<10;l++)
   {
   printf("%e\n",gsl_rng_uniform(rng));
   }
   */
   gsl_rng_set(rng,time(NULL)); // gÃ©nÃ©rateur rÃ©initialisÃ©
   RUN run;
   read_data(argc,argv,&run);
   
  
   
   
   double *pot=NULL;
   double *pot2=NULL;
   double *pot12=NULL;
   double *force_A=NULL; // A(r_ij)
 
   double *arccos=NULL;
   double *legendre=NULL;
   double *hist_gr=NULL;
   double *corr_paire=NULL;
   
   //int *col=NULL; // tableau des couleurs des particules (coordinences)
   //int *ndislo=NULL; // pointeur sur un  entier donnant le nombre de dislocations pour une config
   
  /* double *pos_x=NULL;
   double *pos_y=NULL;
   double *pos_z=NULL;*/ // tableaux des positions suivant x, y et z (sert pour la triangulation)

   double *khi4=NULL; 
 
   POS *pos=NULL; // positions pour Monte Carlo
   POS_ang *pos1=NULL; // positions Ã  t-dt
   POS_ang *pos2=NULL; // positions Ã  t
 
   POS_cart *pos1_cart=NULL; //positions suivant x, y, z Ã  t
   POS_cart *pos2_cart=NULL; //positions suivant x, y, z Ã  t+dt
  
    
   
   POS_cart *omg1_cart=NULL; //vitesses angulaires suivant x, y, z Ã  t
   POS_cart *omg2_cart=NULL; //vitesses angulaires suivant x, y, z Ã  t+dt
  
  
 /* POS_cart *force1_cart=NULL; //forces Ã  t
  POS_cart *force2_cart=NULL;*/ //forces Ã  t+dt
   THERMO thermo;
   thermo.histoener=NULL;
   thermo.histovoro=NULL;
   thermo.n66=NULL;

   allocdouble(&pot,run.maxpot);
   allocdouble(&pot2,run.maxpot);
   allocdouble(&pot12,run.maxpot);
   allocdouble(&arccos,run.maxacos);
   allocdouble(&force_A,run.maxpot);
   allocdouble(&legendre,run.maxlegendre);
   int Ntot=run.N+run.N2;
   double rcut=2.5*run.sigma;
   int nbre_cal=(int)(pi/(4*rcut)); //nbre de calottes autour du pole nord choisi de sorte que sin(thetaup)=pi/4
   double thetaup=nbre_cal*rcut;
//    
   printf("rcut=%g\n",rcut);
 //  printf("Nombre cellules=%d\n\n",Ncell);

   
   int kmax=(int) (3/run.deltar);
   //allocdouble(&hist_gr,kmax); //histogramme caractÃ©risant le nbre de particules par calotte sphÃ©rique
   
  
   
  // allocdouble(&khi4,run.nstep);
   
   allocdouble(&corr_paire,kmax); //tableau qui stocke les valeurs de g(r)
   
   allocdouble(&(thermo.histoener),run.maxhisto);
   
   allocdouble(&(thermo.histovoro),10); // histogramme des couleurs
   allocdouble(&(thermo.n66),1); // n66 moyen
    
   //allocint(&col,run.N); // tableau des couleurs
   
   pos=(POS *) malloc(sizeof(POS_ang)*(Ntot));
   pos1=(POS_ang *) malloc(sizeof(POS_ang)*(Ntot));
   pos2=(POS_ang *) malloc(sizeof(POS_ang)*(Ntot));
   
   pos1_cart=(POS_cart *) malloc(sizeof(POS_cart)*(Ntot)); 
   pos2_cart=(POS_cart *) malloc(sizeof(POS_cart)*(Ntot));

  
   omg1_cart=(POS_cart *) malloc(sizeof(POS_cart)*(Ntot)); 
   omg2_cart=(POS_cart *) malloc(sizeof(POS_cart)*(Ntot));

   

   
   int i;
  
   thermo.flagr=0;
   thermo.ndislo=NULL;
   allocdouble(&(thermo.ndislo),1);
   clock_t startwtime=clock();
   arccos_tabu(arccos,run);
    clock_t endwtime =clock();
    printf("arccos time = %e\n",(double)(endwtime-startwtime)/CLOCKS_PER_SEC);
      startwtime=clock();
   legendre_tabu(legendre,run);
    endwtime =clock();
    printf("legendre time = %e\n",(double)(endwtime-startwtime)/CLOCKS_PER_SEC);
    startwtime=clock();
   initial_pot(pot,pot2,pot12,run);
     endwtime =clock();
    printf("pot time = %e\n",(double)(endwtime-startwtime)/CLOCKS_PER_SEC);
    startwtime=clock();
   initial_force(force_A,run);
   endwtime =clock();
    printf("initial force time = %e\n",(double)(endwtime-startwtime)/CLOCKS_PER_SEC);
   
    
  // int dt_Fs=(int) (run.tFs/500);
  // int dt_Fs_court=(int) (run.tFs/1000);
   //int dt_corr=(int) (run.tFs/100); // temps entre deux configurations dans le calcul de S(k) ou g(r)
           
   dyn_molec(&thermo,pos1_cart,pos2_cart,omg1_cart,omg2_cart,pot,pot2,pot12,force_A,arccos,legendre,run,rng);

   thermo.flagr=1;
   dyn_molec(&thermo,pos1_cart,pos2_cart,omg1_cart,omg2_cart,pot,pot2,pot12,force_A,arccos,legendre,run,rng);

  
  
   
   
   /*!on rÃ©Ã©crit les paramÃ¨tres principaux de la simulation Ã  la fin */
   printf("\n");
   printf("Nbre particles %d\n", run.N);
   printf("Nbre particles 2 %d\n", run.N2);
   printf("Conf avant %ld\n", run.bstep);
   printf("Conf apres %ld\n", run.nstep);
   printf("eta %lf\n", run.eta);
   printf("eta2 %lf\n",run.eta2);
   printf("Temperature %lf\n", run.temp);
   printf("sigma %lf\n", run.sigma);
   printf("sigma2 %lf\n", run.sigma2);
   printf("deltar %f\n", run.deltar);
   
   printf("tFs %d\n", run.tFs);
   printf("flagposanc %d\n", run.flagposanc);
   printf("flagkhi4 %d\n", run.flagkhi4);
   printf("flagdeplac %d\n", run.flagdeplac);
   printf("flagdeplac657 %d\n", run.flagdeplac657);
   printf("flagdevoro %d\n", run.flagdevoro);
   printf("kFs %lf\n",run.kFs);
   printf("kmax %d\n",run.kmax);
   printf("dt %e\n", run.dt);
   printf("nconf %d\n", run.nconf);
   
   printf("flagcorr %d\n", run.flagcorr);
   printf("flagstruc %d\n", run.flagstruc);
   printf("flagg6 %d\n", run.flagg6);
   printf("flagFs %d\n", run.flagFs);
   printf("flagFtot %d\n",run.flagFtot);
   printf("flagdynam %d\n", run.flagdynam);
   printf("flagmovie %d\n", run.flagmovie);
   printf("moviestep %d\n", run.moviestep);
   
   //********
   
  
   gsl_rng_free(rng);
   exit(0);
 } // fin main

   
   
/*!Simulations successives en modifiant les donnees ??*/



void dyn_molec(THERMO *thermo, POS_cart *pos1_cart, POS_cart *pos2_cart, POS_cart *omg1_cart, POS_cart *omg2_cart,
	       double *pot, double *pot2, double *pot12, double *force_A,double *arccos,double *legendre,
	       const RUN run, gsl_rng *rng)
{
  
 
        
  
  
  int Ntot=run.N+run.N2; // nbre total de particules
  int nstep;
  clock_t debuttemp , fintemp ;
   int nthreads;
  if (thermo->flagr==0) {nstep=run.bstep;} 
  else
  {
    nstep=run.nstep;
    double minE=run.nrjbas;
    double maxE=run.nrjhaut*(Ntot);
    double invpasE=run.maxhisto/(maxE-minE);
  }
  
 
  int step;
  int i; int j;
  int k=0; //entier comptant le nombre de realisations dans le calcul de Fs
  int kmax=(int) (3/run.deltar);
  int dt_corr=(int) (run.nstep/run.nconf); // temps entre deux configurations dans le calcul de S(k) ou g(r)
  
  POS_ang *pos1=NULL; // positions a t-dt  
  POS_ang *vit1=NULL; // de type pos mais reprÃ©sente d(theta)/dt, d(phi)/dt a t
  POS_ang *vit2=NULL; // d(theta)/dt, d(phi)/dt a t+dt
  pos1=(POS_ang *) malloc(sizeof(POS_ang)*(Ntot));
  vit1=(POS_ang *) malloc(sizeof(POS_ang)*(Ntot));
  vit2=(POS_ang *) malloc(sizeof(POS_ang)*(Ntot));
  
  POS_cart *tmp; //variable intermediaire qui intervient quand on modifie les positions en fin de boucle
  POS_cart *tmp2; //variable intermediaire qui intervient quand on modifie les vitesses en fin de boucle
 
  static POS_cart *force1_cart=NULL; //forces a t
  static POS_cart *force2_cart=NULL; //forces a t+dt  
  force1_cart=(POS_cart *) malloc(sizeof(POS_cart)*(Ntot));
  force2_cart=(POS_cart *) malloc(sizeof(POS_cart)*(Ntot));
 
    
  int n_part_suivies=20; // on suit n_part_suivies particules au cours du temps
  POS_cart *pos_part_suivie; // tableau des positions des particules que l'on souhaite suivre
  POS_cart *pos_part_suivie_ini; // tableau des positions initiales associÃ©es    
  pos_part_suivie=(POS_cart *) malloc(sizeof(POS_cart)*n_part_suivies);
  pos_part_suivie_ini=(POS_cart *) malloc(sizeof(POS_cart)*n_part_suivies);

  
  double *Fs=NULL; // tableau où on enregistre Fs pour le réutiliser dans le calcul de Khi4
  allocdouble(&Fs,1000);
 
   // deplacement et deplacement carre moyens
   double *hist_deplacement=NULL;
   double *hist_deplacement_prime=NULL; // temps courts (idem que pour Fs)
   allocdouble(&hist_deplacement,1000); 
   allocdouble(&hist_deplacement_prime,500); 
   
   double **hist_deplac_part_suivies=NULL;
   allocdouble2(&hist_deplac_part_suivies,n_part_suivies,500);
   
    clock_t startwtime=clock();
  double T_moy=0;
  
  double *npart=NULL;
  allocdouble(&npart,4);
  
  if((thermo->flagr==0)&&(run.flagposanc==0)) // conditions initiales alÃ©atoires
  {
   dyn_molec_init(pos1_cart, omg1_cart, pos1,vit1,arccos,run,rng); // initialisation des pos et des vitesses
  }
  clock_t endwtime=clock();
   printf("initial dynam  = %e\n",(double)(endwtime-startwtime)/CLOCKS_PER_SEC);
  
  if((thermo->flagr==0)&&(run.flagposanc==1)) // nouvelles positions initiales = anciennes positions finales
  {
   initial_config_anc(pos1_cart,thermo,run);
   reinit_vit(pos1_cart,omg1_cart,pos1,vit1,arccos,run,rng); // initialisation des vitesses
  }


#pragma omp parallel 
  debuttemp=clock();
  omp_set_num_threads(run.nthreads);
  
  
      nthreads=omp_get_num_threads();
      printf("nthreads %d\n",nthreads);  
  for(step=0;step<nstep;step++)
  {    
   double energ_cin=0; // Energie cinÃ©tique moyenne par particule
   double energ_pot=0; // Energie potentielle par particule    
    if((step%500==0)&(thermo->flagr==0)) //rescaling rÃ©gulier avant l'Ã©quilibre
    {
     reinit_vit(pos1_cart,omg1_cart,pos1,vit1,arccos,run,rng);
    }   
#ifdef DEBUG   
#pragma omp parallel 
      nthreads=omp_get_num_threads();
      clock_t debutt=clock();
#endif
      if(step%50==0)
    {
     //printf("position%e\n",pos1_cart[0].x);
     //printf("vitesse%e\n",omg1_cart[0].x);
#pragma omp parallel for private(i) reduction(+:energ_cin)
    for(i=0;i<Ntot;i++)
    {
     energ_cin+=prod_scalaire(omg1_cart[i],omg1_cart[i])+prod_scalaire(omg2_cart[i],omg2_cart[i]); 
    }
    energ_cin*=0.5/((double) (Ntot));    
    energ_pot=1/((double) (Ntot))*calcul_nrj_cart(pos1_cart,pot,pot2,pot12,arccos,run);
    T_moy+=1/((double) nstep/50)*energ_cin; // si calculÃ©e tous les 50 pas de temps
    
    } // fin du if(step%50)
      
    if(step%5000==0)
    {
     printf("step/nstep %d / %d\n",step,(int) nstep);
     printf("Energie cinetique=%e\n",energ_cin);
     printf("Energie potentielle=%e\n",energ_pot);
     printf("Energie totale=%e\n",energ_cin+energ_pot);
    }
      
    verlet(pos1_cart,pos2_cart,omg1_cart,omg2_cart,force1_cart,force2_cart,arccos,force_A,run); 
 #ifdef DEBUG      
    clock_t fint=clock(); 
    double  dureetherm =  (double) ( fint - debutt )  / ((double) nthreads *CLOCKS_PER_SEC);
      ELAPSEDTIME("temps estime... ",dureetherm); // temps estimÃ© initialement - temps Ã©coul
 #endif	      

    
    
    if(step%10000==0)
    {
     printf("flagcorr %d\n", run.flagcorr);
     printf("flagg6 %d\n", run.flagg6);
     printf("flagFs %d\n", run.flagFs);
     printf("flagFtot %d\n",run.flagFtot);
     printf("flagdeplac %d\n", run.flagdeplac);
     printf("flagdeplac657 %d\n", run.flagdeplac657);
     printf("flagdevoro %d\n", run.flagdevoro);
     printf("flagdynam %d\n", run.flagdynam);
     printf("N %d\n",run.N);
     printf("N2 %d\n",run.N2);
     printf("eta %e\n",run.eta);
     printf("eta2 %e\n",run.eta2);
     printf("T %e\n",run.temp);
     printf("bstep %ld\n",run.bstep);
     printf("nstep %ld\n",run.nstep);
     printf("sigma %e\n",run.sigma);
     printf("sigma2 %e\n",run.sigma2);
     printf("deltar %f\n", run.deltar);
     printf("nconf %d\n", run.nconf);
     
     printf("\n");
    }
    
    if(step%1000==0)
    {
//      printf("step=%d\n",step);
//      printf("Energie cinetique=%e\n",energ_cin);
//      printf("Energie potentielle=%e\n",energ_pot);
//      printf("Energie totale=%e\n",energ_cin+energ_pot);
     if((step>0)&&(thermo->flagr==0))
     {
       
       FILE *favanc; // fichier on où on print l'avancement pour simus avec slurm	
	char navanc[50];
	if(run.flagposanc==0)
	{
	 sprintf(navanc,"avancNa%dNb%deta%dTemp%dbstep%dstep%d.res",run.N,run.N2,(int) (1000*run.eta),
		 (int) (1000/run.beta),(int) run.bstep, (int) run.nstep);
	}
	else
	{
	 sprintf(navanc,"avancNa%dNb%deta%dTemp%dbstep%dstep%d.res",run.N,run.N2,(int) (1000*run.eta),
		 (int) (1000/run.beta),(int) run.bstep+ (int) run.bstep_old+ (int) run.step_old,(int) run.nstep);
	}
     
     if ((favanc = fopen (navanc,"w")) == NULL)
        { fprintf(stderr,"Probleme d'ouverture de %s\n",navanc);
          exit (1);  
        }
     
     

#pragma omp parallel 
      nthreads=omp_get_num_threads();
      printf( " threads dyn_molec %d\n",nthreads);
      fintemp=clock();   
       double  dureetherm = ( (double) ( fintemp - debuttemp )*(nstep-step)/step ) / ((double) nthreads *CLOCKS_PER_SEC);
       ELAPSEDTIME("temps estime... ",dureetherm); // temps estimÃ© initialement - temps Ã©coulÃ©
       fprintf(favanc,"temps estime %e s\n",dureetherm);
       
       dureetherm /=((double) nstep-step) ;     
       ELAPSEDTIME("temps pour 1 itération... ",dureetherm); // temps estimÃ© initialement - temps Ã©coulÃ©
       fprintf(favanc,"temps pour 1 itération %e s\n",dureetherm);
       
       fprintf(favanc,"step/nstep %d / %d\n",step,nstep);
       fprintf(favanc,"Energie cinetique %e\n",energ_cin);
       fprintf(favanc,"Energie potentielle %e\n",energ_pot);
       fprintf(favanc,"flagr %d\n",thermo->flagr);
       
       printf("\n");  
       
       //************ écriture dans le fichier avanc
     
    
    
    
    fclose(favanc);
       
     } // fin if (step>0)&&...
     
    if((step>0)&&(thermo->flagr==1))
    {
      
      FILE *favanc; // fichier on où on print l'avancement pour simus avec slurm	
	char navanc[50];
	if(run.flagposanc==0)
	{
	 sprintf(navanc,"avancNa%dNb%deta%dTemp%dbstep%dstep%d.res",run.N,run.N2,(int) (1000*run.eta),
		 (int) (1000/run.beta),(int) run.bstep, (int) run.nstep);
	}
	else
	{
	 sprintf(navanc,"posfinNa%dNb%deta%dTemp%dbstep%dstep%d.res",run.N,run.N2,(int) (1000*run.eta),
		 (int) (1000/run.beta),(int) run.bstep+ (int) run.bstep_old+ (int) run.step_old,(int) run.nstep);
	}
     
     if ((favanc = fopen (navanc,"w")) == NULL)
        { fprintf(stderr,"Probleme d'ouverture de %s\n",navanc);
          exit (1);  
        }
      
#pragma omp parallel
     nthreads=omp_get_num_threads();
     printf( "Nombre de threads %d\n",nthreads); 
     if(step<10000)
     {
       fintemp=clock();   
       double  dureetherm = ( (double) ( fintemp - debuttemp )*(10000-step)/step ) / ((double) nthreads *CLOCKS_PER_SEC);
       ELAPSEDTIME("temps estime (Fs court)... ",dureetherm); // temps estimÃ© initialement - temps Ã©coulÃ©
       fprintf(favanc,"temps estime %e\n ",dureetherm);
       
       dureetherm /=((double) 10000-step) ;
       ELAPSEDTIME("temps pour 1 itération... ",dureetherm); // temps estimÃ© initialement - temps Ã©coulÃ©
       fprintf(favanc,"temps pour 1 itération %e\n",dureetherm);
       
       fprintf(favanc,"step/nstep %d / %d\n",step,nstep);
       fprintf(favanc,"Energie cinetique %e\n",energ_cin);
       fprintf(favanc,"Energie potentielle %e\n",energ_pot);
       fprintf(favanc,"flagr %d\n",thermo->flagr);
       
       printf("\n");  
     }
     
     if(step==10000)
     {
      debuttemp=clock(); // on réinitialise le temps pour la nouvelle estimation avec Fs aux temps longs
     }
     if(step>10000)
     {
      fintemp=clock();   
       double  dureetherm = ( (double) ( fintemp - debuttemp )*(nstep-step)/step ) / ((double) nthreads *CLOCKS_PER_SEC);
       ELAPSEDTIME("temps estime (Fs long)... ",dureetherm); // temps estimÃ© initialement - temps Ã©coulÃ©
       fprintf(favanc,"temps estime %e s\n ",dureetherm);
       
       dureetherm =( (double) ( fintemp - debuttemp )/step ) / ((double) nthreads *CLOCKS_PER_SEC); // temps d'une itération
       ELAPSEDTIME("temps pour 1 itération... ",dureetherm); // temps estimÃ© initialement - temps Ã©coulÃ©
       fprintf(favanc,"temps pour 1 itération %e s\n",dureetherm);
       
       fprintf(favanc,"step/nstep %d / %d\n",step,nstep);
       fprintf(favanc,"Energie cinetique %e\n",energ_cin);
       fprintf(favanc,"Energie potentielle %e\n",energ_pot);
       fprintf(favanc,"flagr %d\n",thermo->flagr);
       
       printf("\n");
     }
     
     fclose(favanc);
     
    } // fin if (step>0)&&
    
    
    }   // fin if (step==1000) 
    if(thermo->flagr==1)
    {
   // if(step%500==0)  { fprintf(fenerg,"%e\t %e\n",step*run.dt/run.sigma, energ_cin+energ_pot);}
    }
    
    //********** enregistrement des positions
    
    if (thermo->flagr==1) // aprÃ¨s thermalisation
    {
      int nbre_pts_Fs=50;
      int nbre_pts_Ftot=50;
      int nbre_pts=50; // nombre de configs sur lesquelles on moyenne
      int dt_Fs=(int) (run.tFs/nbre_pts); //pas de temps dans le calcul de Fs (la moyenne temporelle dans le calcul de Fs est faite sur [0,2*tFs])
      //int dt_Fs=(int) (run.nstep/nbre_pts);
      int dt_Ftot=(int) (run.tFs/nbre_pts);
      int dt_Fs_court=1;
      int dt_defauts=(int) (run.nstep/20000); // temps entre l'enregistrement des positions des défauts
      //int dt_deplac=(int) (run.tFs/1000);
      int nlissage=20;
      int nfrac=20000/nlissage;
      int i;
      
      
      if(step==0) /*! initialisation de ttes les fonctions */
      {	
	save_defauts(0,0,pos2_cart,thermo,run);
	spectre_defauts(0,0,pos2_cart,thermo,arccos,run);
	calcul_nrj_pot(0,0,pos2_cart,pot,pot2,pot12,thermo,arccos,run);
	calcul_nrj_pot_lissee(0,0,pos2_cart,pot,pot2,pot12,thermo,arccos,run);
	calcul_dist20part(0,0,pos2_cart,arccos,run);
	calcul_distrib_paire(0,0,pos2_cart,arccos,run,rng);
	calcul_corr_defauts(0,0,pos2_cart,arccos,run,rng);
	//calcul_g6(0,0,pos2_cart,arccos,run,rng);
	calcul_Fs(0,0,dt_Fs,pos2_cart,Fs,arccos,legendre,run);
	calcul_Ftot(0,0,dt_Ftot,pos2_cart,arccos,legendre,run);
	calcul_deplac(0,0,dt_Fs,pos2_cart,arccos,run);
	//calcul_deplac2(0,0,dt_Fs,pos2_cart,arccos,run); // <r^2(t)>
	calcul_deplac657(0,0,dt_Fs,pos2_cart,arccos,thermo,run);
	calcul_khi4(0,0,pos2_cart,Fs,arccos,legendre,run);
	calcul_structure_factor(0,0,pos2_cart,arccos,run,rng);
	
      }
      
      if(((step % dt_Fs_court)==0)&&(step<10000)) // calcul de Fs aux temps courts
      {
       if(run.flagFs==1)
       {
	calcul_Fs(1,step,dt_Fs_court,pos2_cart,Fs,arccos,legendre,run);
       }
       
       
       if(run.flagdeplac==1) /*! si on calcule <r(t)> */
       {
	if(run.flagFs==0)
        {
         calcul_deplac(1,step,dt_Fs_court,pos2_cart,arccos,run);
	// calcul_deplac2(1,step,dt_Fs_court,pos2_cart,arccos,run);
        }
       }
       if(run.flagdeplac657==1) /*! si on calcule <r(t)> */
       {
	if(run.flagFs==0)
        {
       //  calcul_deplac657(1,step,dt_Fs_court,pos2_cart,arccos,thermo,run);
        }
       }
      }
          
       if((run.flagcorr==1)&&(run.flagg6==0)&&(step%dt_corr==0)) // si on calcule déjà g6 on n'a pas besoin de calculer g(r) en plus
       {
        calcul_distrib_paire(1,step,pos2_cart,arccos,run,rng);
	//calcul_g6(1,step,pos2_cart,arccos,run,rng);
       }
       
//        if((run.flagcorrdefauts==1)&&(run.flagg6==0)&&(step%dt_corr==0)) // si on calcule déjà g6 on n'a pas besoin de calculer g(r) en plus
//        {
// 	calcul_corr_defauts(1,step,pos2_cart,arccos,run,rng);
//        }
       
       
       if((run.flagg6==1)&&(step%dt_corr==0))
       {
	//calcul_g6(1,step,pos2_cart,arccos,run,rng); 
       }
       
       if((run.flagstruc==1)&&(step%dt_corr==0))
       {
	calcul_structure_factor(1,step,pos2_cart,arccos,run,rng);
       }
       
       if(run.flagdevoro==1)
       {
	if(step%run.moviestep==0)
	{
	 if(run.flagspectre==1)
	 {
	  spectre_defauts(1,step,pos2_cart,thermo,arccos,run);	 
	  calcul_nrj_pot(1,step,pos2_cart,pot,pot2,pot12,thermo,arccos,run);
	  //calcul_nrj_pot_lissee(1,step,pos2_cart,pot,pot2,pot12,thermo,arccos,run);
	 }
	// energ_pot=1/((double) (Ntot))*calcul_nrj_cart(pos1_cart,pot,pot2,pot12,arccos,run);
	// fprintf(fenerg,"%e\t %e\n",step*run.dt/run.sigma,energ_pot);
	 if(run.flagmovie==1)
	 {
	  calcul_images(1,step,pos2_cart,thermo,run); 
	 }
	}
	
	if(step%nlissage==0)
	{
	 if(run.flagspectre==1)
	 {
	  calcul_nrj_pot_lissee(1,step,pos2_cart,pot,pot2,pot12,thermo,arccos,run); 
	 }
	}
	
       }
	 
       
      if((step % dt_Fs)==0)
      {
      // calcul_ndefauts(1,pos2_cart,npart,arccos,run);
       if(run.flagdeplac==1) //****************** si on calcule <r(t)> et 20 trajectoires
       {
	calcul_dist20part(1,step,pos2_cart,arccos,run); 
	if(run.flagFs==0)
        {
         calcul_deplac(1,step,dt_Fs,pos2_cart,arccos,run);
	// calcul_deplac2(1,step,dt_Fs,pos2_cart,arccos,run);
        }
       }
       if(run.flagdeplac657==1) //****************** si on calcule <r(t)> et 20 trajectoires
       {
	if(run.flagFs==0)
        {
         calcul_deplac657(1,step,dt_Fs,pos2_cart,arccos,thermo,run);
        }
       }
       /*! si on calcule Fs	              */
       if(run.flagFs==1) // calcul de Fs aux temps longs
       {
        calcul_Fs(1,step,dt_Fs,pos2_cart,Fs,arccos,legendre,run);
       }
      } // fin du if((step % dt_Fs)==0)
      
      if((step % dt_Ftot)==0) //***************** si on calcule Ftot
      {
       if(run.flagFtot==1)
       {
	calcul_Ftot(1,step,dt_Ftot,pos2_cart,arccos,legendre,run); 
       }
       
       /*! si on calcule Khi4 */
       
       if(run.flagkhi4==1)
       {
	calcul_khi4(1,step,pos2_cart,Fs,arccos,legendre,run);
       }
       
      }
      
      if((step % dt_defauts==0))
      {
	if(run.flagsavedefauts==1)
	{
	 save_defauts(1,step,pos2_cart,thermo,run);
	}
      }
 
      
      
    } // fin du if(thermo->flagr==1)
    
  
   tmp=pos1_cart;
   pos1_cart=pos2_cart; // quand t -> t+dt r,i(t-dt) -> r,i(t)
   pos2_cart=tmp; // on modifie le pointeur pos2 pour que les pointeurs pos1 et pos2 ne pointent pas sur la mÃªme pas case (pas de souci
   // car les valeurs de pos2 seront redÃ©finies au pas suivant
   
   
   tmp2=omg1_cart;
   omg1_cart=omg2_cart;
   omg2_cart=tmp2;
    
   //if(thermo->flagr==1) printf("OK\n");
    } //fin de la boucle sur step
   
 //  printf("OK\n");
   if ((thermo->flagr==1)&&(run.flagdevoro==1)) // aprÃ¨s thermalisation
   {
   for(i=0;i<10;i++)
	{
	 printf(" %d %e\n",i,thermo->histovoro[i]/((double) run.nstep)*run.moviestep);
	}	
   }
   
   double nbre_defauts=0;
   for(i=0;i<10;i++)
   {
     if(i!=6)
     {
      int charge_i=fabs(6-i); // on pondÃ¨re le nombre de dÃ©fauts de chaque type par sa charge topologique en valeur absolue
      nbre_defauts+=charge_i*thermo->histovoro[i]/((double) run.nstep)*run.moviestep;
     }
   }
   
   printf("Fraction de defauts thermiques = %e\n",1/((double) (Ntot))*(nbre_defauts-12));
   printf("nbre de dislocations 5-7 / 12 = %e\n",(double) 1/12*thermo->ndislo[0]);
   printf("Temperature moyenne = %e\n",T_moy);
   printf("\n");
   
   /*! enregistrement des données dans des fichiers */
   
   int nbre_pts=50; // nombre de configs sur lesquelles on moyenne
   
   //
   if(thermo->flagr==1)
   {
    if(run.flagsavedefauts==1)
    {
     save_defauts(2,step,pos2_cart,thermo,run);
    }
     
    if((run.flagcorr==1)&&(run.flagg6==0))
    {
     calcul_distrib_paire(2,step,pos2_cart,arccos,run,rng); // calcul de g(r) Ã  la fin de la simu
     
    }
    
//     if((run.flagcorrdefauts==1)&&(run.flagg6==0))
//     {
//      calcul_corr_defauts(2,step,pos2_cart,arccos,run,rng);
//     }
    if(run.flagg6==1)
    {
     calcul_g6(2,step,pos2_cart,arccos,run,rng); 
    }
    
     if(run.flagstruc==1)
     {
      calcul_structure_factor(2,step,pos2_cart,arccos,run,rng);
     }
    
    
   if((run.flagdevoro==1)&&(run.flagspectre==1))
   {
   // spectre_defauts(2,step,pos2_cart,thermo,arccos,run); 
    calcul_nrj_pot(2,step,pos2_cart,pot,pot2,pot12,thermo,arccos,run);
    calcul_nrj_pot_lissee(2,step,pos2_cart,pot,pot2,pot12,thermo,arccos,run);
   }
   
   if(run.flagFs==1)
   {
    int nbre_pts=50;
    int dt_Fs=(int) (run.tFs/nbre_pts);
   // int dt_Fs=(int) (run.nstep/nbre_pts);
    calcul_Fs(2,step,dt_Fs,pos2_cart,Fs,arccos,legendre,run);
   }
   if(run.flagFtot==1)
   {
    int nbre_pts_Ftot=50;
    int dt_Ftot=(int) (run.tFs/nbre_pts);
   // int dt_Ftot=(int) (run.nstep/nbre_pts);
    calcul_Ftot(2,step,dt_Ftot,pos2_cart,arccos,legendre,run);
   }
    
   if(run.flagkhi4==1)
   {
    calcul_khi4(2,step,pos2_cart,Fs,arccos,legendre,run);
   }
      	
   if(run.flagdeplac==1)
   {
    calcul_dist20part(2,step,pos2_cart,arccos,run);
    if(run.flagFs==0)
    {
     int nbre_pts_Fs=50;
     int dt_Fs=(int) (run.tFs/nbre_pts);
     calcul_deplac(2,step,dt_Fs,pos2_cart,arccos,run);
    // calcul_deplac2(2,step,dt_Fs,pos2_cart,arccos,run);
    }
   }
   
   if(run.flagdeplac657==1)
   {
    if(run.flagFs==0)
    {
     int nbre_pts=50;
     int dt_Fs=(int) (run.tFs/nbre_pts);
     //int dt_Fs=(int) (run.nstep/nbre_pts);
     calcul_deplac657(2,step,dt_Fs,pos2_cart,arccos,thermo,run);
    }
   }
   
   } // fin if thermo->flagr=1
 
   // fichier oÃ¹ sont enregistrÃ©es les positions finales (Ã  bstep + nstep)
   
       
        FILE *fpos;	
	char npos[50];
	if(run.flagposanc==0)
	{
	 sprintf(npos,"posfinNa%dNb%deta%dTemp%dbstep%dstep%d.res",run.N,run.N2,(int) (1000*run.eta),
		 (int) (1000/run.beta),(int) run.bstep, (int) run.nstep);
	}
	else
	{
	 sprintf(npos,"posfinNa%dNb%deta%dTemp%dbstep%dstep%d.res",run.N,run.N2,(int) (1000*run.eta),
		 (int) (1000/run.beta),(int) run.bstep+ (int) run.bstep_old+ (int) run.step_old,(int) run.nstep);
	}
        if ((fpos = fopen (npos,"w")) == NULL)
        { fprintf(stderr,"Probleme d'ouverture de %s\n",npos);
          exit (1);  
        }
	 
	 int ipar;
	 for(ipar=0;ipar<Ntot;ipar++)
	 {
	   fprintf(fpos,"%.4f\t %.4f\t %.4f\n",pos2_cart[ipar].x,pos2_cart[ipar].y,pos2_cart[ipar].z);
	 }
	 
	fclose(fpos);
	
	/*! fichier de résultats Tmoy, Fraction de disinclinaisons, nbre de dislocations, nbre de disinclinaisons 4 à 7*/
	FILE *fres;	
	char nres[50];
	if(run.flagposanc==0)
	{
	 sprintf(nres,"resultatsNa%dNb%deta%dTemp%dbstep%dstep%d.res",run.N,run.N2,(int) (1000*run.eta),
		(int) (1000/run.beta),(int) run.bstep, (int) run.nstep);
	}
	else
	{
	 sprintf(nres,"resultatsNa%dNb%deta%dTemp%dbstep%dstep%d.res",run.N,run.N2,(int) (1000*run.eta),
		(int) (1000/run.beta),(int) run.bstep+ (int) run.bstep_old+ (int) run.step_old,(int) run.nstep);
	}
        if ((fres = fopen (nres,"w")) == NULL)
        { fprintf(stderr,"Probleme d'ouverture de %s\n",nres);
          exit (1);  
        }
	 
         fprintf(fres,"Nbre particles %d\n", run.N);
         fprintf(fres,"Nbre particles 2 %d\n", run.N2);
         fprintf(fres,"Conf avant %ld\n", run.bstep);
         fprintf(fres,"Conf apres %ld\n", run.nstep);
         fprintf(fres,"eta %lf\n", run.eta);
         fprintf(fres,"eta2 %lf\n",run.eta2);
         fprintf(fres,"Temperature %lf\n", run.temp);
         fprintf(fres,"sigma %lf\n", run.sigma);
         fprintf(fres,"sigma2 %lf\n", run.sigma2);
         fprintf(fres,"deltar %f\n", run.deltar);  
         fprintf(fres,"tFs %d\n", run.tFs);         
         fprintf(fres,"flagkhi4 %d\n", run.flagkhi4);
         fprintf(fres,"kFs %lf\n",run.kFs);
         fprintf(fres,"dt %e\n", run.dt);
	 printf("nconf %d\n", run.nconf);
	 fprintf(fres,"flagposanc %d\n", run.flagposanc);
         fprintf(fres,"flagcorr %d\n", run.flagcorr);
	 fprintf(fres,"flagstruc %d\n", run.flagstruc);
         fprintf(fres,"flagFs %d\n", run.flagFs);
         fprintf(fres,"flagdeplac %d\n", run.flagdeplac);
	 fprintf(fres,"flagdeplac657 %d\n", run.flagdeplac657);
         fprintf(fres,"flagdevoro %d\n", run.flagdevoro);
	 fprintf(fres,"flagsavedefauts %d\n", run.flagsavedefauts);
	 fprintf(fres,"flagspectre %d\n", run.flagspectre);
         fprintf(fres,"flagdynam %d\n", run.flagdynam);
         fprintf(fres,"flagmovie %d\n", run.flagmovie);
	 fprintf(fres,"moviestep %d\n", run.moviestep);
	 fprintf(fres,"Tmoy\t%e\n", T_moy);
	 if ((thermo->flagr==1)&&(run.flagdevoro==1)) // aprÃ¨s thermalisation
         {
	  fprintf(fres,"Frac defauts thermiques\t %e\n", 1/((double) (Ntot))*(nbre_defauts-12));
	  fprintf(fres,"nbre de dislocations 5-7 / 12 \t %e\n",(double) 1/12*thermo->ndislo[0]);	
          for(i=0;i<10;i++)
	  {
	   fprintf(fres, "ndisinc %d %e\n",i,thermo->histovoro[i]/((double) run.nstep)*run.moviestep);
	  }	
         }
	 fclose(fres);
	
 }/*!fin de dyn_molec */

