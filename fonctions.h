
//version 27 avril 2015

double min_acos; // valeur minimale dans tabulation de arccos : min_acos=cos(2.5*sigma)
double max_acos; // valeur maximale dans tabulation de arccos : max_acos=cos(0.5*sigma)



/*! enregistrement positions des défauts*/

void save_defauts(int flag, int step, POS_cart *pos2_cart, THERMO *thermo, RUN run)
{
 int i;
 static int *col=NULL; // tableau des couleurs des particules (coordinences)
 static int *ndislo=NULL; // pointeur sur un  entier donnant le nombre de dislocations pour une config
 static int *voisins=NULL;
 static FILE *fdefauts;
 char ndefauts[200];
 static double *pos_x=NULL;
 static double *pos_y=NULL;
 static double *pos_z=NULL;
 int dt_defauts=(int) (run.nstep/10000);
  
 
 if(flag==0)
 {
  allocdouble(&pos_x,run.N+run.N2);
  allocdouble(&pos_y,run.N+run.N2);
  allocdouble(&pos_z,run.N+run.N2);
  allocint(&col,run.N+run.N2);
  allocint(&ndislo,1);
  allocint(&voisins,(run.N+run.N2)*8);
  
  
   //FILE *fdefauts;
    //char ndefauts[50];
          
  
  
 } // fin if flag==0
 
 if(flag==1)
 {
   int t=step/dt_defauts; // t=0,1,..., 10000
   for(i=0;i<run.N+run.N2;i++)
   {
    pos_x[i]=pos2_cart[i].x;
    pos_y[i]=pos2_cart[i].y;
    pos_z[i]=pos2_cart[i].z;
   }
  
  voronoi(pos_x, pos_y, pos_z,voisins,col,run.N+run.N2,ndislo);
  
  for(i=0;i<run.N+run.N2;i++) // ne pas paralléliser car on ajoute dans un tableau
  {	
   thermo->histovoro[col[i]]++;	 
  }
  
  thermo->ndislo[0]+=ndislo[0]/((double) run.nstep)*dt_defauts;
  /*! on normalise par le nombre d'instants sur lesquels on moyenne  */
  
  if(run.flagposanc==0)
    {
     sprintf(ndefauts,"positions_defauts/posdefautsNa%dNb%deta%dTemp%dbstep%dstep%dt%d",
	     run.N,run.N2,(int) (1000*run.eta),(int) (1000/run.beta),(int) run.bstep,(int) run.nstep,t);
    }
    else  //sinon le nouveau bstep est bstep+step_oldsimu
    {
     sprintf(ndefauts,"positions_defauts/posdefautsNa%dNb%deta%dTemp%dbstep%dstep%dt%d",
	     run.N,run.N2,(int) (1000*run.eta),(int) (1000/run.beta),(int) (run.bstep+run.bstep_old+run.step_old),(int) (run.nstep),t);
    }
     if ((fdefauts = fopen (ndefauts,"w")) == NULL)
    { fprintf(stderr,"Probleme d'ouverture de %s\n",ndefauts);
      exit (1);  
    }
  
  
  
    for(i=0;i<run.N+run.N2;i++)
    {
     if(col[i]!=6)
     {
      fprintf(fdefauts,"%.4f\t %.4f\t %.4f\t %d\n",  pos2_cart[i].x, pos2_cart[i].y, pos2_cart[i].z,col[i]);
     }
    }
    fclose(fdefauts);
 } // fin if flag==1
 
 if(flag==2)
 {
   FILE *ffonction_voro;
  char voro[50];
  if(run.flagposanc==0)
  {

    sprintf(voro,"voroNa%dNb%deta%dTemp%dbstep%dstep%dtFs%d.res",run.N,run.N2,(int) (1000*run.eta),(int) (1000/run.beta),(int) run.bstep,(int) run.nstep,(int) (run.tFs));
  }
  else  //sinon le nouveau bstep est bstep+step_oldsimu
  {
    sprintf(voro,"voroNa%dNb%deta%dTemp%dbstep%dstep%dtFs%d.res",run.N,run.N2,(int) (1000*run.eta),(int) (1000/run.beta),(int) (run.bstep+run.bstep_old+run.step_old),(int) (run.nstep),(int) (run.tFs));
  }
    
    
    
  if ((ffonction_voro = fopen (voro,"w")) == NULL)
    {
     fprintf(stderr,"Probleme d'ouverture de %s\n",voro);
     exit (1);  
    }
   
   double N_tot; /*! nombre total de particules = somme des n_i, n_i couleur de la parti (pour vérifier)*/
     for(i=0;i<10;i++)
     {
      N_tot+=thermo->histovoro[i]/((double) run.nstep)*dt_defauts;
     }
     
     
     
     
   for(i=0;i<10;i++)
   {
    fprintf(ffonction_voro,"%d\t %e\n",  i, thermo->histovoro[i]/N_tot*dt_defauts/((double) run.nstep));
   }
    fprintf(ffonction_voro,"ndislo/12 %e\n",  thermo->ndislo[0]/12);
    fprintf(ffonction_voro,"Ntot (verif) %e\n",  N_tot);
   
    fclose(ffonction_voro); 
 }
 
} /*! fin save defauts*/











  //************************************************
  

  

 
 
 
 
 






/*! calcul de r_i(t+dt) */

inline void verlet(POS_cart *pos1_cart, POS_cart *pos2_cart, POS_cart *omg1_cart, 
	     POS_cart *omg2_cart, POS_cart *force1_cart, POS_cart *force2_cart,
		  double *arccos, double *force_A, RUN run)
{
  
   int nrun=run.N;
   int nrun2=run.N2;
    const double prod_scal_max=cos(2.5*run.sigma);
    const double prod_scal_max2=cos(2.5*run.sigma2);
    const double prod_scal_max12=cos(2.5*run.sigma12);
    
     double upper=2.5*run.sigma;
     double upper2=2.5*run.sigma2;
     int invupper=run.maxpot/upper;
     int invupper2=run.maxpot/upper2;
   double startwtime,endwtime;
  int i;
#pragma omp parallel 
// if (omp_get_thread_num() == 0)   startwtime =  omp_get_wtime();
#pragma omp parallel for default(shared), private(i)  
   for(i=0;i<run.N+run.N2;i++) /*! on met d'abord Ã  jour ttes les positions avant de mettre Ã  jour les vitesses */
  {
   force1_cart[i]=force2_cart[i]; /*! quand t -> t'=t+dt l'ancienne force Ã  t+dt devient la nouvelle force Ã  t */
   POS_cart A=somme_vect(scalaire(run.dt,omg1_cart[i]), scalaire(-0.5*run.dt*run.dt, prod_vect(pos1_cart[i],force1_cart[i])));
   if(prod_scalaire(A,A)>1){printf("i=%d\t A=%e\n",i,prod_scalaire(A,A));}
   pos2_cart[i]=somme_vect(prod_vect(A,pos1_cart[i]), scalaire(sqrt(1-prod_scalaire(A,A)),pos1_cart[i]));   
  }/*! fin boucle sur i */
  
 
  
 int j; 
 
 /*! calcul des forces exercées sur les particules 1 */
#pragma omp parallel for shared(force2_cart,pos2_cart,nrun,nrun2,prod_scal_max,prod_scal_max2),  private(j) 
 for(i=0;i<nrun;i++)
 {
   POS_cart force={0,0,0};
   for(j=1;j<nrun;j++)
   {
     int j1=(j+i) % (nrun);
    // printf("i %d\t j1 %d\n",i,j1);
//       if(j!=i)
//       {
    //double prod_scal=pos2_cart[i].x*pos2_cart[j].x+pos2_cart[i].y*pos2_cart[j].y+pos2_cart[i].z*pos2_cart[j].z;
      double prod_scal=prod_scalaire(pos2_cart[i],pos2_cart[j1]);
       if(prod_scal>prod_scal_max) 
       {	
        if ((prod_scal<=-1) || (prod_scal>=1)) {printf("Prob force 1->1 i j  prod_scal %d %d %e\n",i,j,prod_scal);}      
        double dist=acos(GSL_MIN_DBL(prod_scal,1.0));       
     // double dist=arccosinus(GSL_MIN_DBL(prod_scal,1.0),arccos,run);      
        double r=dist/run.sigma;
        double r2=r*r;
        double r7=r2*r2*r2*r;
        double r13=r7*r2*r2*r2;
        double tmp=(12/run.sigma)*(2.0/r13-1.0/r7)*1/(sin(dist))/*-force_shift*/; //A(r) (faut-il shifter la force ?) 
        force=somme_vect(force, scalaire(tmp, pos2_cart[j1]));
       }
    // }
      
  } // fin boucle sur j

   
   for(j=nrun;j<nrun+nrun2;j++)
   {
     double prod_scal=prod_scalaire(pos2_cart[i],pos2_cart[j]);
     if(prod_scal>prod_scal_max12)
     {	
      if ((prod_scal<=-1) || (prod_scal>=1)) {printf("Prob force 2->1 i j  prod_scal %d %d %e\n",i,j,prod_scal);}      
       double dist_2=acos(GSL_MIN_DBL(prod_scal,1.0));           
       double r=dist_2/run.sigma12;
       double r2=r*r;
       double r7=r2*r2*r2*r;
       double r13=r7*r2*r2*r2;
       double tmp=(12/run.sigma12 )*(2.0/r13-1.0/r7)*1/(sin(dist_2))/*-force_shift*/; //A(r) (faut-il shifter la force ?) 
      force=somme_vect(force, scalaire(tmp, pos2_cart[j]));
   //  }
     }
    } // fin boucle sur j
 force2_cart[i]=force;
 } // fin boucle sur i
 
 //*********** i>N
 #pragma omp parallel for shared(force2_cart,pos2_cart,nrun,nrun2,prod_scal_max12,prod_scal_max2),  private(j) 
 for(i=nrun;i<nrun+nrun2;i++)
 {
   POS_cart force={0,0,0};
   for(j=1;j<nrun;j++)
   {
     //int j1=(j+i) % (nrun);
//       if(j!=i)
//       {
    //double prod_scal=pos2_cart[i].x*pos2_cart[j].x+pos2_cart[i].y*pos2_cart[j].y+pos2_cart[i].z*pos2_cart[j].z;
      double prod_scal=prod_scalaire(pos2_cart[i],pos2_cart[j]);
       if(prod_scal>prod_scal_max12) 
       {	
        if ((prod_scal<=-1) || (prod_scal>=1)) {printf("Prob force 1->1 i j  prod_scal %d %d %e\n",i,j,prod_scal);}      
        double dist=acos(GSL_MIN_DBL(prod_scal,1.0));       
     // double dist=arccosinus(GSL_MIN_DBL(prod_scal,1.0),arccos,run);      
        double r=dist/run.sigma12;
        double r2=r*r;
        double r7=r2*r2*r2*r;
        double r13=r7*r2*r2*r2;
        double tmp=(12/run.sigma12)*(2.0/r13-1.0/r7)*1/(sin(dist))/*-force_shift*/; //A(r) (faut-il shifter la force ?) 
        force=somme_vect(force, scalaire(tmp, pos2_cart[j]));
       }
    // }
      
  } // fin boucle sur j

   
   for(j=nrun+1;j<nrun+nrun2;j++)
   {
      int j1=nrun+(i+j-2*nrun) % (nrun2);
   //   printf("i %d\t j1 %d\n",i,j1);
//      if(j!=i)
//      {
    //double prod_scal=pos2_cart[i].x*pos2_cart[j].x+pos2_cart[i].y*pos2_cart[j].y+pos2_cart[i].z*pos2_cart[j].z;
     double prod_scal=prod_scalaire(pos2_cart[i],pos2_cart[j1]);
     if(prod_scal>prod_scal_max2)
     {	
      if ((prod_scal<=-1) || (prod_scal>=1)) {printf("Prob force 2->1 i j  prod_scal %d %d %e\n",i,j,prod_scal);}      
       double dist_2=acos(GSL_MIN_DBL(prod_scal,1.0));       
     // double dist=arccosinus(GSL_MIN_DBL(prod_scal,1.0),arccos,run);      
       double r=dist_2/run.sigma2;
       double r2=r*r;
       double r7=r2*r2*r2*r;
       double r13=r7*r2*r2*r2;
       double tmp=(12/run.sigma2 )*(2.0/r13-1.0/r7)*1/(sin(dist_2))/*-force_shift*/; //A(r) (faut-il shifter la force ?) 
      force=somme_vect(force, scalaire(tmp, pos2_cart[j1]));
     }
   //  }
    } // fin boucle sur j
 force2_cart[i]=force;
 } // fin boucle sur i

 
 
 
#pragma omp parallel for  shared(pos1_cart,pos2_cart,force1_cart,force2_cart,run), private(i)
  for(i=0;i<nrun+nrun2;i++)
  {
   omg2_cart[i]=somme_vect(omg1_cart[i], scalaire(-0.5*run.dt, somme_vect(prod_vect(pos1_cart[i],force1_cart[i]), prod_vect(pos2_cart[i],force2_cart[i])))); 
  }
  

 if (omp_get_thread_num() == 0)    {endwtime =  omp_get_wtime();
//  printf("verlet time = %f\n",(double) (endwtime-startwtime));
    startwtime=endwtime;}


 
} // fin verlet



/*! g(r) */

void calcul_distrib_paire(int flag, int step, POS_cart *pos2_cart, double *arccos, RUN run, gsl_rng *rng)
{
  int k;
  int kmax=(int) (3/run.deltar); //valeur max de k
  double rho=run.N/(4*pi);
  double rho_2=run.N2/(4*pi);
  double rho_12=rho+rho_2;
  static double *corr_paire=NULL;
  static double *corr_paire_2=NULL;
  static double *corr_paire_12=NULL;
  static double *hist_gr=NULL;
  static double *hist_gr_2=NULL;
  static double *hist_gr_12=NULL;
  const double prod_scal_max=cos(2.5*run.sigma);
  
  if(flag==0)
  {
  allocdouble(&corr_paire,kmax); //tableau qui stocke les valeurs de g(r)
  allocdouble(&hist_gr,kmax); //histogramme caracterisant le nbre de particules par calotte spherique
  if(run.N2>0)
  {
   allocdouble(&hist_gr_2,kmax);
   allocdouble(&corr_paire_2,kmax);
   allocdouble(&hist_gr_12,kmax);
   allocdouble(&corr_paire_12,kmax);
  }
  }
  
  if(flag==1)
  {
   int i1=(int)((run.N) *gsl_rng_uniform(rng)); //a  reinitialiser tous les N pas de temps
   int j; //indice des particules
   int k; //indice des couronnes

   for(j=0;j<run.N;j++)
   {
    for(k=0;k<kmax;k++)
    {
      double prod_scal=prod_scalaire(pos2_cart[i1],pos2_cart[j]);
      if(prod_scal>prod_scal_max)
      {
      double dist=distance_cart(pos2_cart[i1],pos2_cart[j],arccos,run);
     if((j!=i1)&&((dist-k*run.deltar)*(dist
      -(k+1)*run.deltar)<0))
     {
      hist_gr[k]++;
     }
      }
    }  
   } 
   
   if(run.N2>0)
   {
    int i2=run.N+(int)((run.N2) *gsl_rng_uniform(rng));
   for(j=run.N;j<run.N+run.N2;j++)
   {
    for(k=0;k<kmax;k++)
    {
     if((j!=i2)&&((distance_cart(pos2_cart[i2],pos2_cart[j],arccos,run)-k*run.deltar)*(distance_cart(pos2_cart[i2],pos2_cart[j],arccos,run)
      -(k+1)*run.deltar)<0))
     {
      hist_gr_2[k]++;
     }
    }  
   }
   
    for(j=0;j<run.N;j++)
    {
     for(k=0;k<kmax;k++)
     {
      if((distance_cart(pos2_cart[i2],pos2_cart[j],arccos,run)-k*run.deltar)*(distance_cart(pos2_cart[i2],pos2_cart[j],arccos,run)
      -(k+1)*run.deltar)<0)
      {
       hist_gr_12[k]++;
      }
     } 
    }
   
   
   }
   
  } // fin if(flag==1)
  
  if(flag==2)
  {
  FILE *fcorrelation;
  char ncorr[50];
  
  int nstep=run.nstep;
  int bstep=run.bstep;
  
  if(run.flagposanc==0)
  { 
  sprintf(ncorr,"correlationaaNa%dNb%deta%dTemp%dbstep%dstep%d.res",run.N,run.N2,(int) (1000*run.eta),(int) (1000/run.beta),bstep,nstep);
  }
  else
  {
  sprintf(ncorr,"correlationaaNa%dNb%deta%dTemp%dbstep%dstep%d.res",run.N,run.N2,(int) (1000*run.eta),(int) (1000/run.beta),(int) (run.bstep_old+run.step_old+bstep),(int) (nstep));  
  }
  
  if ((fcorrelation = fopen (ncorr,"w")) == NULL)
  { 
   fprintf(stderr,"Probleme d'ouverture de %s\n",ncorr);
   exit (1);  
  }
 
  for(k=0;k<kmax;k++)
  {
   corr_paire[k]=hist_gr[k]/(rho*run.nconf*2*pi*(cos(k*run.deltar)-cos((k+1)*run.deltar)));   
   fprintf(fcorrelation, "%e %e\n", run.deltar*(k+0.5)/run.sigma, corr_paire[k]);
  }
  
  fclose(fcorrelation);
  
  if(run.N2>0)
  {
   FILE *fcorrelation_2;
   char ncorr_2[50];
  
   int nstep=run.nstep;
   int bstep=run.bstep;
  
  if(run.flagposanc==0)
  { 
  sprintf(ncorr_2,"correlationbbNa%dNb%deta%dTemp%dbstep%dstep%d.res",run.N,run.N2,(int) (1000*run.eta),(int) (1000/run.beta),bstep,nstep);
  }
  else
  {
  sprintf(ncorr_2,"correlationbbNa%dNb%deta%dTemp%dbstep%dstep%d.res",run.N,run.N2,(int) (1000*run.eta),(int) (1000/run.beta),(int) (run.bstep_old+run.step_old+bstep),(int) (nstep));  
  }
  
  if ((fcorrelation_2 = fopen (ncorr_2,"w")) == NULL)
  { 
   fprintf(stderr,"Probleme d'ouverture de %s\n",ncorr_2);
   exit (1);  
  }
 
  for(k=0;k<kmax;k++)
  {
   corr_paire_2[k]=hist_gr_2[k]/(rho_2*run.nconf*2*pi*(cos(k*run.deltar)-cos((k+1)*run.deltar)));   
   fprintf(fcorrelation_2, "%e %e\n", run.deltar*(k+0.5)/run.sigma2, corr_paire_2[k]);
  }
  
  fclose(fcorrelation_2); 
  
  /*! gAB(r) */
  FILE *fcorrelation_12;
   char ncorr_12[50];

  
  if(run.flagposanc==0)
  { 
  sprintf(ncorr_12,"correlationabNa%dNb%deta%dTemp%dbstep%dstep%d.res",run.N,run.N2,(int) (1000*run.eta),(int) (1000/run.beta),bstep,nstep);
  }
  else
  {
  sprintf(ncorr_12,"correlationabNa%dNb%deta%dTemp%dbstep%dstep%d.res",run.N,run.N2,(int) (1000*run.eta),(int) (1000/run.beta),(int) (run.bstep_old+run.step_old+bstep),(int) (nstep));  
  }
  
  if ((fcorrelation_12 = fopen (ncorr_12,"w")) == NULL)
  { 
   fprintf(stderr,"Probleme d'ouverture de %s\n",ncorr_12);
   exit (1);  
  }
 
  for(k=0;k<kmax;k++)
  {
   corr_paire_12[k]=hist_gr_12[k]/(rho*run.nconf*2*pi*(cos(k*run.deltar)-cos((k+1)*run.deltar)));   
   fprintf(fcorrelation_12, "%e %e\n", run.deltar*(k+0.5)/run.sigma12, corr_paire_12[k]);
  }
  
  fclose(fcorrelation_12); 
  
  
  }// fin if N2>0
  
  } //fin if flag ==2
 
} //fin calcul_distrib_paire



/*! fonction de corrélation des défauts g_q(r) */

void calcul_corr_defauts(int flag, int step, POS_cart *pos2_cart, double *arccos, RUN run, gsl_rng *rng)
{
  int i,k;
  int kmax=(int) (3/run.deltar); //valeur max de k
  static double *corr_defauts=NULL;
  static double *hist_corr_defauts=NULL;
  static double *pos_x=NULL;
  static double *pos_y=NULL;
  static double *pos_z=NULL;
  static int *col=NULL; // tableau des couleurs des particules (coordinences)
  static int *ndislo=NULL; // pointeur sur un  entier donnant le nombre de dislocations pour une config
  static int *voisins=NULL;
  
  if(flag==0)
  {
   allocdouble(&corr_defauts,kmax); //tableau qui stocke les valeurs de g(r)
   allocdouble(&hist_corr_defauts,kmax); //histogramme caractÃ©risant le nbre de particules par calotte 
   allocdouble(&pos_x,run.N+run.N2);
   allocdouble(&pos_y,run.N+run.N2);
   allocdouble(&pos_z,run.N+run.N2);
   allocint(&col,run.N+run.N2);
   allocint(&ndislo,1);
   allocint(&voisins,(run.N+run.N2)*8);
  }
  if(flag==1)
  {
   int i1=(int)((run.N) *gsl_rng_uniform(rng)); //Ã  rÃ©initialiser tous les N pas de temps
   int j; //indice des particules
   int k; //indice des couronnes

   
   for(i=0;i<run.N+run.N2;i++)
   {
    pos_x[i]=pos2_cart[i].x;
    pos_y[i]=pos2_cart[i].y;
    pos_z[i]=pos2_cart[i].z;
   } 
   voronoi(pos_x, pos_y, pos_z,voisins,col,run.N+run.N2,ndislo);
   
   for(j=0;j<run.N;j++)
   {
    for(k=0;k<kmax;k++)
    {
     double dist=distance_cart(pos2_cart[i1],pos2_cart[j],arccos,run);
     if((j!=i1)&&((dist-k*run.deltar)*(dist-(k+1)*run.deltar)<0))
     {
      double qi=6-col[i1];
      double qj=6-col[j];
      hist_corr_defauts[k]+=qi*qj;
     }
    }  
   } 
   
  } // fin if(flag==1)
  
  if(flag==2)
  {
  FILE *fcorrdefauts;
  char ncorrdefauts[50];
  
  int nstep=run.nstep;
  int bstep=run.bstep;
  
  if(run.flagposanc==0)
  { 
  sprintf(ncorrdefauts,"correlationdefautsaaNa%dNb%deta%dTemp%dbstep%dstep%d.res",run.N,run.N2,(int) (1000*run.eta),(int) (1000/run.beta),bstep,nstep);
  }
  else
  {
  sprintf(ncorrdefauts,"correlationdefautsaaNa%dNb%deta%dTemp%dbstep%dstep%d.res",run.N,run.N2,(int) (1000*run.eta),(int) (1000/run.beta),(int) (run.bstep_old+run.step_old+bstep),(int) (nstep));  
  }
  
  if ((fcorrdefauts = fopen (ncorrdefauts,"w")) == NULL)
  { 
   fprintf(stderr,"Probleme d'ouverture de %s\n",ncorrdefauts);
   exit (1);  
  }
 
  for(k=0;k<kmax;k++)
  {
   corr_defauts[k]=hist_corr_defauts[k]/(run.nconf*2*pi*(cos(k*run.deltar)-cos((k+1)*run.deltar)));   // normalisation taille de couronne
   fprintf(fcorrdefauts, "%e %e\n", run.deltar*(k+0.5)/run.sigma, fabs(corr_defauts[k]));
  }
  
  fclose(fcorrdefauts);
  
  
  } //fin if flag ==2
  
  
}



/*! fonction calculant S(k) */


void calcul_structure_factor(int flag, int step, POS_cart *pos2_cart, double *arccos, RUN run, gsl_rng *rng)
{
  int k;
  int k_wave=2*(int)(run.kFs*pi/run.sigma);
  int i;
  static int nconf;
  int kmax=run.kmax; //valeur max de k
  double rho=run.N/(4*pi);
  static double *structure_factor=NULL;
 // static double *hist_gr=NULL;
  
  if(flag==0)
  {
   allocdouble(&structure_factor,kmax); //tableau qui stocke les valeurs de g(r)
   nconf=0;
  }
  
  if(flag==1)
  {
    for(k=1;k<kmax;k++) // on calcule tous les S(k)
    {
     double rho_k=0; // densité
     for(i=0;i<run.N;i++)
     {
      double z=pos2_cart[i].z;
      rho_k+=gsl_sf_legendre_Pl(k,z);
     }
     structure_factor[k]+=rho_k*rho_k;
     
    }  
     nconf++;
    // printf("nconf %d\n", nconf);  
  } // fin if(flag==1)  
  
  if(flag==2)
  {
  FILE *fstructure;
  char nstructure[50];
  
  int nstep=run.nstep;
  int bstep=run.bstep;
  
  if(run.flagposanc==0)
  { 
  sprintf(nstructure,"strucfactorN%deta%dTemp%dbstep%dstep%dkmax%d.res",run.N,(int) (1000*run.eta),(int) (1000/run.beta),bstep,nstep,run.kmax);
  }
  else
  {
  sprintf(nstructure,"strucfactorN%deta%dTemp%dbstep%dstep%dkmax%d.res",run.N,(int) (1000*run.eta),(int) (1000/run.beta),
	  (int) (run.bstep_old+run.step_old+bstep),(int) (nstep),run.kmax);  
  }
  
  if ((fstructure = fopen (nstructure,"w")) == NULL)
  { 
   fprintf(stderr,"Probleme d'ouverture de %s\n",nstructure);
   exit (1);  
  }
 
  for(k=1;k<kmax;k++)
  {
   structure_factor[k]=(2*k+1)/((double) run.N)*1/((double) nconf)*structure_factor[k]; // on normalise S(k)=(2k+1)/N*1/Nconf*...
   //double c_k=4*3.14159/(double run.N)*(1-1/structure_factor[k]); //c(k) corrélation directe
   fprintf(fstructure, "%d\t %e\n", k, structure_factor[k]);
  }
  
  fclose(fstructure);
  
  
  
  } //fin if flag ==2
  
  
} // fin calcul_structure_static

// //************************ calcul grain boundary scar
// 
// void calcul_scar(int flag, double *arccos, RUN run, gsl_rng *rng)
// {
//  static int *col=NULL; // tableau des couleurs des particules (coordinences)
//  static int *voisins=NULL;
//  static int *ndislo=NULL; // pointeur sur un  entier donnant le nombre de dislocations pour une config
//  static int *scar_size=NULL; // tableau des tailles de la scar à différents instants
//  static int *scar_list1=NULL; // tableau indiquant quelles particules sont dans la scar (scar_list[i]=1 si i dans la scar)
//  static int *scar_list2=NULL;
//  static int *nscar;
//  static int *scar_col=NULL;
//  static double *pos_x=NULL;
//  static double *pos_y=NULL;
//  static double *pos_z=NULL;
// 
//    
//  const double prod_scal_max_scar=cos(3*run.sigma);
//  //printf("prod_scal_max_scar %e\n", prod_scal_max_scar);
//  
// //  static double *pos_defauts_x=NULL; // positions de tous les défauts à t, suivant x 
// //  static double *pos_defauts_y=NULL; // positions de tous les défauts à t, suivant y 
// //  static double *pos_defauts_z=NULL; // positions de tous les défauts à t, suivant z 
//  
// //  static double **scar_pos_x=NULL; // positions des particules de la scar suivant x
// //  static double **scar_pos_y=NULL;
// //  static double **scar_pos_z=NULL;
//  
//  POS_cart *pos_defauts=NULL; // positions de tous les defauts
//  POS_cart *pos_scar=NULL; // positions des défauts dans 1 scar
//  int t;
//  
//  if(flag==0)
//   {
//    
//    allocdouble(&pos_x,run.N+run.N2);
//    allocdouble(&pos_y,run.N+run.N2);
//    allocdouble(&pos_z,run.N+run.N2);
// //    allocdouble(&pos_defauts_x,run.N+run.N2);
// //    allocdouble(&pos_defauts_y,run.N+run.N2);
// //    allocdouble(&pos_defauts_z,run.N+run.N2);
//    
//    allocint(&col,run.N+run.N2);
//    allocint(&ndislo,1);
//    allocint(&voisins,(run.N+run.N2)*8);
//    allocint(&scar_size,500);
//    allocint(&nscar,1);
// //    allocint(&scar_list1,run.N+run.N2);
// //    allocint(&scar_list2,run.N+run.N2);
//    
//    
//    
//    
// //    allocdouble2(&scar_pos_x,run.N+run.N2,500);
// //    allocdouble2(&scar_pos_y,run.N+run.N2,500);
// //    allocdouble2(&scar_pos_z,run.N+run.N2,500);
//    
//    //******************* on cherche une scar dans tout le système
//    
//    int i;
//    int j;
//    //******************************* à t=0
//  
//  for(t=1;t<50;t++)
//  {
//    
//   
//    FILE *fdefauts;
//    char ndefauts[100];
//    
//   if(run.flagposanc==0)
//     {
//      sprintf(ndefauts,"positions_defauts/posdefautsNa%dNb%deta%dTemp%dbstep%dstep%dtFs%dt%d",
// 	     run.N,run.N2,(int) (1000*run.eta),(int) (1000/run.beta),(int) run.bstep,(int) run.nstep,(int) (run.tFs),t);
//     }
//     else  //sinon le nouveau bstep est bstep+step_oldsimu
//     {
//      sprintf(ndefauts,"positions_defauts/posdefautsNa%dNb%deta%dTemp%dbstep%dstep%dtFs%dt%d",
// 	     run.N,run.N2,(int) (1000*run.eta),(int) (1000/run.beta),(int) (run.bstep+run.bstep_old+run.step_old),(int) (run.nstep),(int) (run.tFs),t);
//     }
//      if ((fdefauts = fopen (ndefauts,"r")) == NULL)
//     { fprintf(stderr,"Probleme d'ouverture de %s\n",ndefauts);
//       exit (1);  
//     }
//   
//   
//   
//   //**************** calcul du nbre total de défauts comme le nombre de lignes dans fdefauts
//   int ndefauts_tot=line_count(fdefauts)-1;
//   
// 
//     
//   fclose(fdefauts);
//   
//   //ndefauts_tot=5; // test
//   printf("t %d\n",t);
//   printf("ndefauts %d\n",ndefauts_tot);
//  //**************** 
//     
//   pos_defauts=(POS_cart *) malloc(sizeof(POS_cart)*(ndefauts_tot));
//   allocint(&scar_list1,ndefauts_tot);
//   allocint(&scar_col,ndefauts_tot);
//   
//   if ((fdefauts = fopen (ndefauts,"r")) == NULL)
//     { fprintf(stderr,"Probleme d'ouverture de %s\n",ndefauts);
//       exit (1);  
//     }
//   
//    for(i=0;i<ndefauts_tot;i++)
//    {
//     //fscanf(fdefauts,"%le\t %le\n",&pos_init_phi[i],&pos_init_z[i]);
//     fscanf(fdefauts,"%le\t %le\t %le\t %d\n",&pos_defauts[i].x,&pos_defauts[i].y,&pos_defauts[i].z,&scar_col[i]);
//     //printf("pos.phi=%e pos.z=%e\n",pos[i].phi,pos[i].z);2
//    }
//   fclose(fdefauts);
//    
//    
//  //  voronoi(pos_x, pos_y, pos_z,voisins,col,run.N+run.N2,ndislo);
//    
// //***************** recherche d'une scar
//    
// 
//    nscar[0]=0;
//   // int scar_size=nscar[0];
//    
//    if(t==0)
//    {
//     while(nscar[0]<5)
//    {
//     nscar[0]=0;
//     int i1=(int)((ndefauts_tot) *gsl_rng_uniform(rng)); //attention il faudra toujours suivre la même scar !!
//     for(j=0;j<ndefauts_tot;j++){scar_list1[j]=0;} // on réinitialise tant qu'on n'a pas le bon i1
//    
//      calcul_voisins_defauts(i1,ndefauts_tot,nscar,scar_list1,pos_defauts,arccos,voisins,col,run,rng);
//     // scar_size=nscar[0];
//      printf("nscar %d\n\n", nscar[0]);
//    }// fin while(nscar<8)
//    }
//    
//    if(t>=1)
//    {
//      
//    }
//    
//    
//   
//   
//   POS_cart posm=calcul_pos_scar(t,ndefauts_tot,nscar, 
// 		     scar_list1,pos_defauts,arccos,voisins,col,run,rng);
//   printf("x_scar %e\t y_scar %e\t z_scar %e\n", posm.x, posm.y, posm.z);
//   printf("\n");
//   
//   
//   
//   } // fin for(t=...)
//   } // fin if flag==0
//   
//   
// //   if(flag==1)
// //   {
// //   
// //    
// // 
// //   } // fin if flag==1
//   
//   }// fin for(t=0;500;t++)
//  
// }
// 
// void calcul_voisins_defauts(int i, int ndefauts, int *nscar, 
//  int *scar_list, POS_cart *pos_defauts, double *arccos, int *voisins, int *col, RUN run, gsl_rng *rng) 
// // fonction renvoyant la liste des indices des defauts dans un rayon de 4 sigma de i
// {
//  const double prod_scal_max_scar=cos(8*run.sigma);
//  int j;
//  
//  printf("on ajoute i=%d a la scar \n",i);
//  
//  //***** test ecriture recursive
//  
//  scar_list[i]=1; // on "marque" i
//  nscar[0]++;
// // printf("nscar %d\n\n", nscar);
//  
//  for(j=0;j<ndefauts;j++)
//  {
//   if(j!=i)
//   {
//     double prod_scal=prod_scalaire(pos_defauts[i],pos_defauts[j]);
//     if(prod_scal<prod_scal_max_scar) // on teste si dist(i,j)<3*sigma
//     {  
//     double dist=distance_cart(pos_defauts[i],pos_defauts[j],arccos,run);
//    // double dist=1;
//          if(dist<8*run.sigma) // test si j est un "fils" de i
//          {
//           if(scar_list[j]==0)
// 	  {
// 	   calcul_voisins_defauts(j,ndefauts,nscar,scar_list,pos_defauts,arccos,voisins,col,run,rng);
// 	  }
//          }
//     }
//   }
//  }
// } // fin calcul_voisins_defauts
// 
// 
// //************ calcul de fonctions caractérisant la scar
// 
// POS_cart calcul_pos_scar(int step,int ndefauts, int *nscar, 
// 		     int *scar_list, POS_cart *pos_defauts, double *arccos, int *voisins, int *col, RUN run, gsl_rng *rng)
// {
//   
//  POS_ang posm_ang;
//  POS_cart posm;
//  int i;
//  
//  for(i=0;i<nscar[0];i++)
//  {
//   POS_ang posi;
//   posi.theta=acos(pos_defauts[i].z);
//   posi.phi=atan(pos_defauts[i].y/pos_defauts[i].x);
//   posm_ang.theta+=posi.theta;
//   posm_ang.phi+=posi.phi;
//  }
//  posm_ang.phi=1/((double) nscar[0])*posm_ang.phi;
//  posm_ang.theta=1/((double) nscar[0])*posm_ang.theta;
//  
//  double theta=posm_ang.theta;
//  double phi=posm_ang.phi;
//  
//  posm.x=sin(theta)*cos(phi);
//  posm.y=sin(theta)*sin(phi);
//  posm.z=cos(theta);
//  
//  return(posm);
// }



/*! calcul G6(r) (ne marche pas) */

void calcul_g6(int flag, int step, POS_cart *pos2_cart, double *arccos, RUN run, gsl_rng *rng)
{
  double rho=(run.N+run.N2)/(4*pi);
  int kmax=(int) (3/run.deltar); //valeur max de k
  POS_cart ez;
  ez.x=0; ez.y=0; ez.z=1;
  static double *corr_paire=NULL;
  //static complexe *g6=NULL;
  static double *corr_g6=NULL;
  static double *hist_gr=NULL;
  //static complexe *hist_g6=NULL;
  static double *hist_g6_Re=NULL;
  static double *hist_g6_Im=NULL;
  static double *pos_x=NULL;
  static double *pos_y=NULL;
  static double *pos_z=NULL;
  static int *col=NULL; // tableau des couleurs des particules (coordinences)
  static int *ndislo=NULL; // pointeur sur un  entier donnant le nombre de dislocations pour une config
  static int *voisins=NULL;
  
  complexe i_complexe;
  i_complexe.Re=0;
  i_complexe.Im=1;

  
  if(flag==0)
  {
   allocdouble(&hist_gr,kmax); //histogramme caracterisant le nbre de particules par calotte spherique
   allocdouble(&corr_paire,kmax);
   allocdouble(&hist_g6_Re,kmax);
   allocdouble(&hist_g6_Im,kmax);
   allocdouble(&corr_g6,kmax); 
   int k;
//    for(k=0;k<kmax;k++)
//    {
//     hist_g6[k].Re=hist_g6_Re[k];
//     hist_g6[k].Im=hist_g6_Im[k];
//    }
   allocdouble(&pos_x,run.N+run.N2);
   allocdouble(&pos_y,run.N+run.N2);
   allocdouble(&pos_z,run.N+run.N2);
   allocint(&col,run.N+run.N2);
   allocint(&ndislo,1);
   allocint(&voisins,(run.N+run.N2)*8);
  }
  
  if(flag==1)
  {
   int i;
   int i1=(int)((run.N) *gsl_rng_uniform(rng)); //Ã  rÃ©initialiser tous les N pas de temps
   int j; //indice des particules
   int l; //indice des voisins de j dans la couronne [r,r+dr]
   int k; //indice des couronnes
   
   for(i=0;i<run.N+run.N2;i++)
   {
    pos_x[i]=pos2_cart[i].x;
    pos_y[i]=pos2_cart[i].y;
    pos_z[i]=pos2_cart[i].z;
   } 
   voronoi(pos_x, pos_y, pos_z,voisins,col,run.N+run.N2,ndislo);
   
   complexe psi_6_0;
   psi_6_0.Re=0.; psi_6_0.Im=0.;
   for(l=0;l<8;l++)
      {
       int vois=voisins[8*i+l];
       POS_cart r_il=somme_vect(pos2_cart[l],scalaire(-1,pos2_cart[i]));
       double theta_il=acos(prod_scalaire(r_il,ez)/sqrt(prod_scalaire(r_il,r_il)));
       psi_6_0=somme_complexe(psi_6_0,exp_complexe(scalaire_complexe((double) 6*theta_il,i_complexe)));
      }
      psi_6_0=scalaire_complexe(1/((double) col[i]),psi_6_0); // on divise par le nbre de voisins de j
   
   for(j=0;j<run.N+run.N2;j++)
   {
    for(k=0;k<kmax;k++)
    {
      int npart_cour=0;
      complexe histo_g6;
      histo_g6.Re=0; histo_g6.Im=0;
      
     if((j!=i1)&&((distance_cart(pos2_cart[i1],pos2_cart[j],arccos,run)-k*run.deltar)*(distance_cart(pos2_cart[i1],pos2_cart[j],arccos,run)
      -(k+1)*run.deltar)<0))
     { 
      hist_gr[k]++;
      complexe psi_6_r;
      psi_6_r.Re=0; psi_6_r.Im=0;
      for(l=0;l<8;l++)
      {
       int vois=voisins[8*j+l];
       if(vois>0)
       {
        POS_cart r_jl=somme_vect(pos2_cart[l],scalaire(-1,pos2_cart[j]));
        double theta_jl=acos(prod_scalaire(r_jl,ez)/sqrt(prod_scalaire(r_jl,r_jl)));
	psi_6_r=somme_complexe(psi_6_r,exp_complexe(scalaire_complexe(6*theta_jl,i_complexe)));
       }
      }
      psi_6_r=scalaire_complexe(1/((double) col[j]),psi_6_r); // on divise par le nbre de voisins de j
      histo_g6=somme_complexe(histo_g6,prod_complexe(psi_6_0,psi_6_r));
      npart_cour++;
     }
//      hist_g6_Re[k]=somme_complexe(hist_g6_Re[k],scalaire_complexe(1/((double) npart_cour),histo_g6)).Re;
//      hist_g6_Im[k]=somme_complexe(hist_g6_Im[k],scalaire_complexe(1/((double) npart_cour),histo_g6)).Im;
    }    
   } // fin if(flag==1)
   
  if(flag==2)
  {
  FILE *fg6;
  char ng6[50];
  
  int nstep=run.nstep;
  int bstep=run.bstep;
  
  if(run.flagposanc==0)
  { 
  sprintf(ng6,"g6aaNa%dNb%deta%dTemp%dbstep%dstep%d.res",run.N,run.N2,(int) (1000*run.eta),(int) (1000/run.beta),bstep,nstep);
  }
  else
  {
  sprintf(ng6,"g6aaNa%dNb%deta%dTemp%dbstep%dstep%d.res",run.N,run.N2,(int) (1000*run.eta),(int) (1000/run.beta),(int) (run.bstep_old+run.step_old+bstep),(int) (nstep));  
  }
  
  if ((fg6 = fopen (ng6,"w")) == NULL)
  { 
   fprintf(stderr,"Probleme d'ouverture de %s\n",ng6);
   exit (1);  
  }
 
  for(k=0;k<kmax;k++)
  {
   //corr_paire[k]=hist_gr[k]/(rho*run.nstep*2*pi*(cos(k*run.deltar)-cos((k+1)*run.deltar)));  
   corr_paire[k]=hist_gr[k]/(rho*run.nstep*2*pi*(cos(k*run.deltar)-cos((k+1)*run.deltar)));
   complexe g6;
   g6.Re=hist_g6_Re[k];
   g6.Im=hist_g6_Im[k];
   
   corr_g6[k]=module_complexe(scalaire_complexe(1/((double) run.nstep*corr_paire[k]),g6));
   fprintf(fg6, "%e %e\n", run.deltar*(k+0.5)/run.sigma, corr_g6[k]);
  }
  
  fclose(fg6);
    
  } //fin if flag ==2
  
}
} // fin calcul_g6


/*! Enregistrement des pos pour film */

void calcul_images(int flag, int step, POS_cart *pos2_cart, THERMO *thermo, RUN run)
{
 int i; int j;
 static double *pos_x=NULL;
 static double *pos_y=NULL;
 static double *pos_z=NULL;
 static int *col=NULL; // tableau des couleurs des particules (coordinences)
 static int *ndislo=NULL; // pointeur sur un  entier donnant le nombre de dislocations pour une config
 static int *voisins=NULL;
  
 if(flag==0)
 {
    allocdouble(&pos_x,run.N+run.N2); 
    allocdouble(&pos_y,run.N+run.N2); 
    allocdouble(&pos_z,run.N+run.N2); 
    allocint(&col,run.N+run.N2);
    allocint(&ndislo,1);
    allocint(&voisins,(run.N+run.N2)*8);
  //  printf("intialisation calcul_defauts OK\n");
 }
 
 if(flag==1)
 {
  
#pragma omp parallel for private(i)
  for(i=0;i<run.N+run.N2;i++)
  {
   pos_x[i]=pos2_cart[i].x;
   pos_y[i]=pos2_cart[i].y;
   pos_z[i]=pos2_cart[i].z;
  }
#pragma omp barrier
  for(i=0;i<run.N+run.N2;i++)
   {
    for(j=0;j<8;j++)
    {
     voisins[8*i+j]=0;
    }
   }

  voronoi(pos_x, pos_y, pos_z,voisins,col,run.N+run.N2,ndislo);
  
  
  
  if(run.flagmovie==1)
  {
  if(step%run.moviestep==0)
  {
  
  FILE *fd;
  char movie[100];
  sprintf(movie,"films/movNa%dNb%deta%dTemp%dbstep%dstep%ima%d",run.N,run.N2,(int) (1000*run.eta),(int) (1000/run.beta),(int) run.bstep, (int) run.nstep,step/run.moviestep);
      
  if ((fd = fopen (movie,"w")) == NULL)
    { fprintf(stderr,"Probleme d'ouverture de %s\n",movie);
      exit (1);  
    }
  for(i=0;i<run.N+run.N2;i++)
  {
    fprintf(fd,"%.4f\t %.4f\t %.4f\t %d\n",pos_x[i],pos_y[i],pos_z[i],col[i]);
  }
  fclose(fd);
  }
  }
 } // fin if flag==1
 
} // fin calcul_images


/*! calcul des défauts topologiques, P_6.. */

void spectre_defauts(int flag, int step, POS_cart *pos2_cart, THERMO *thermo, double *arccos, RUN run)
{
 int i, j, k, t, n6, ndefauts, n_dislo_appar, f;
 static double *pos_x=NULL;
 static double *pos_y=NULL;
 static double *pos_z=NULL;
 static int *col=NULL; // tableau des couleurs des particules (coordinences)
 static int *ndislo=NULL; // pointeur sur un  entier donnant le nombre de dislocations pour une config
 static int *voisins=NULL;
 
 static double *frac6_fonc_t=NULL; // fraction de particules 6 au cours de t
 static double *frac_dislo_fonc_t=NULL; // fraction de 5 ou 7 associés en dislocations parmi les défauts
 static double *spectre_frac6=NULL;
 static int *n66=NULL;
 int nfrac=1000; // nombre de points pour frac6 et spectre_frac6
 int nstep=(int) (run.nstep);
 const double prod_scal_max=cos(2.5*run.sigma);
 
 if(flag==0)
 {
  allocdouble(&pos_x,run.N+run.N2); 
  allocdouble(&pos_y,run.N+run.N2); 
  allocdouble(&pos_z,run.N+run.N2); 
  allocint(&col,run.N+run.N2);
  allocint(&ndislo,1);
  allocint(&n66,1);
  allocint(&voisins,(run.N+run.N2)*8);
  allocdouble(&frac6_fonc_t,nfrac);
  allocdouble(&frac_dislo_fonc_t,nfrac);
  allocdouble(&spectre_frac6,nfrac);
  t=0;
  //  printf("intialisation calcul_defauts OK\n");
 }
 
 if(flag==1)
 {
  t=step/run.moviestep;
  
#pragma omp parallel for private(i)
  for(i=0;i<run.N+run.N2;i++)
  {
   pos_x[i]=pos2_cart[i].x;
   pos_y[i]=pos2_cart[i].y;
   pos_z[i]=pos2_cart[i].z;
  }
#pragma omp barrier
  for(i=0;i<run.N+run.N2;i++)
   {
    for(j=0;j<8;j++)
    {
     voisins[8*i+j]=0;
    }
   }

  voronoi(pos_x, pos_y, pos_z,voisins,col,run.N+run.N2,ndislo);

  n6=0;
  #pragma omp parallel for private(i) reduction(+:n6)
  for(i=0;i<run.N+run.N2;i++)
  {
   if(col[i]==6) n6++; 
  }
  #pragma omp barrier
  
  ndefauts=run.N-n6;

  
  n_dislo_appar=0;
  
  
   
   frac6_fonc_t[t]=n6/((double) run.N);
  if(ndislo[0]>0) frac_dislo_fonc_t[t]=n_dislo_appar/(2*(double) ndislo[0]); // il faut diviser par 2 car chaque dislo appariee comptee 2 fois
 //  printf("t %d\t frac6 %e\t frac57 %e\n",t, frac6_fonc_t[t], frac_dislo_fonc_t[t]);
 //  t++;

   
 } // fin if flag == 1
 
 if(flag==2)
 {
 
/*! séries temporelles des fractions de 6 et défauts */ // commenté car erreur segmentation
 
 FILE *ffonction_defauts;
  char defauts[50];
 if(run.flagposanc==0)
  {

    sprintf(defauts,"defautsNa%dNb%deta%dTemp%dbstep%dstep%dtFs%d.res",run.N,run.N2,
    (int) (1000*run.eta),(int) (1000/run.beta),(int) run.bstep,(int) run.nstep,(int) (run.tFs));
  }
  else  //sinon le nouveau bstep est bstep+step_oldsimu
  {
    sprintf(defauts,"defautsNa%dNb%deta%dTemp%dbstep%dstep%dtFs%d.res",run.N,run.N2,
    (int) (1000*run.eta),(int) (1000/run.beta),(int) (run.bstep+run.bstep_old+run.step_old),(int) (run.nstep),(int) (run.tFs));
  }
 
    
  if ((ffonction_defauts = fopen (defauts,"w")) == NULL)
    {
     fprintf(stderr,"Probleme d'ouverture de %s\n",defauts);
     exit (1);  
    }
    
  
  printf("nstep %d\n",nstep);
  //int dt_defauts2=nstep/nfrac;
  int dt_defauts2=run.tFs/nfrac;
  double delta_f=1/(nfrac*dt_defauts2);
 // long dt_defauts2=run.nstep/nfrac;
 for(t=0;t<nfrac;t++)
 {
  fprintf(ffonction_defauts, "%e\t %e\t %e\n",t*run.dt/run.sigma*dt_defauts2, frac6_fonc_t[t], frac_dislo_fonc_t[t]); 
 }
 fclose(ffonction_defauts); 
 
 //********************* TFourier
 
// int delta_t=nfrac/run.N;
 //int delta_t=dt_defauts;

complexe sum;
complexe i0; i0.Re=0; i0.Im=1;
 int kt, kf;

 
for(kf=0;kf<nfrac;kf++) 
{
  //double f=pi2/(run.N*delta_t)*i;
  
  
 double sum_cos; // sum_cos=Somme(0<=kt<=nfrac){P_6(kt)cos(2*pi*knu*kt/nfrac)}
 double sum_sin; // sum_sin=Somme(0<=kt<=nfrac){P_6(kt)sin(2*pi*knu*kt/nfrac)}
 for(kt=0;kt<nfrac;kt++)
 {
  sum_cos+=frac6_fonc_t[kt]*cos(pi2*kt*kf/nfrac);
  sum_sin+=frac6_fonc_t[kt]*sin(pi2*kt*kf/nfrac);
 }
 spectre_frac6[kf]=1/((double) nfrac*nfrac)*(sum_cos*sum_cos+sum_sin*sum_sin);
 
}  // fin for kf

 
  FILE *fspectre_defauts;
  char spectre_defauts[50];
 if(run.flagposanc==0)
  {

    sprintf(spectre_defauts,"spectre_defautsNa%dNb%deta%dTemp%dbstep%dstep%dtFs%d.res",run.N,run.N2,(int) (1000*run.eta),(int) (1000/run.beta),(int) run.bstep,(int) run.nstep,(int) (run.tFs));
  }
  else  //sinon le nouveau bstep est bstep+step_oldsimu
  {
    sprintf(spectre_defauts,"spectre_defautsNa%dNb%deta%dTemp%dbstep%dstep%dtFs%d.res",run.N,run.N2,(int) (1000*run.eta),(int) (1000/run.beta),(int) (run.bstep+run.bstep_old+run.step_old),(int) (run.nstep),(int) (run.tFs));
  }
 
    
  if ((fspectre_defauts = fopen (spectre_defauts,"w")) == NULL)
    {
     fprintf(stderr,"Probleme d'ouverture de %s\n",spectre_defauts);
     exit (1);  
    }
    
    nstep=(int) (run.nstep);
    nfrac=1000; 
  //  dt_defauts2=nstep/nfrac;
    dt_defauts2=run.tFs/nfrac;
    delta_f=1/((double) nfrac*dt_defauts2);
    printf("nfrac %d\n", nfrac);
    printf("dt_defauts %d\n", dt_defauts2);
    printf("sigma %e\n", run.sigma);
    printf("dt %e\n", run.dt);
    printf("freq_max %e\n",run.sigma/(run.dt*dt_defauts2));
    
 // int dt_defauts=run.nstep/1000;
 for(kf=0;kf<nfrac;kf++)
 {
  
 // fprintf(fspectre_defauts, "%e\t %e\n",1/((double) nfrac*run.dt*delta_t*dt_defauts)*run.sigma*kf, spectre_frac6[kf]);
   fprintf(fspectre_defauts, "%e\t %e\n",1/((double) nfrac*run.dt*dt_defauts2)*run.sigma*kf, spectre_frac6[kf]); 
 }

 fclose(fspectre_defauts); 

 
 } // fin if flag == 2 
} // fin calcul_defauts



/*! deplacement moyen total */

void calcul_deplac(int flag, int step, int dt_Fs, POS_cart *pos2_cart,double *arccos, RUN run)
{
 int i;
 int j;
 int nbre_box=50; // nombre de points calculés
 int dt_deplac=(int) (run.tFs/nbre_box);
 int itt1=step/dt_deplac;
 int itt0=itt1% nbre_box;
 int itt2;
 int itt3;
 int tFs=run.tFs;
 //int dt_Fs=tFs/500;
 int dt_Fs_court=1;
 
 POS_cart pos1;
 POS_cart pos2;
 
 static double **hist_pos_x_deplac_long=NULL; //histo des coordonnees phi (pour calcul des grandeurs dynamiques)
 static double **hist_pos_y_deplac_long=NULL;
 static double **hist_pos_z_deplac_long=NULL;
 static double **hist_pos_x_deplac_court=NULL; //histo des coordonnees phi (pour calcul des grandeurs dynamiques)
 static double **hist_pos_y_deplac_court=NULL;
 static double **hist_pos_z_deplac_court=NULL;
 
 static double *pos1_x=NULL;
 static double *pos1_y=NULL;
 static double *pos1_z=NULL;
 static double *pos2_x=NULL;
 static double *pos2_y=NULL;
 static double *pos2_z=NULL;
 
 static double *hist_deplac_long=NULL;
 static double *hist_deplac_court=NULL;
 static int compteur_deplac_court;
 static int compteur_deplac_long;

 
 if(flag==0)
 {
 allocdouble2(&hist_pos_x_deplac_long,run.N+run.N2,nbre_box);
 allocdouble2(&hist_pos_y_deplac_long,run.N+run.N2,nbre_box);
 allocdouble2(&hist_pos_z_deplac_long,run.N+run.N2,nbre_box);
 allocdouble2(&hist_pos_x_deplac_court,run.N+run.N2,nbre_box);
 allocdouble2(&hist_pos_y_deplac_court,run.N+run.N2,nbre_box);
 allocdouble2(&hist_pos_z_deplac_court,run.N+run.N2,nbre_box);
 allocdouble(&hist_deplac_long,nbre_box);
 allocdouble(&hist_deplac_court,nbre_box);
 
 allocdouble(&pos1_x,run.N+run.N2);
 allocdouble(&pos1_y,run.N+run.N2);
 allocdouble(&pos1_z,run.N+run.N2);
 allocdouble(&pos2_x,run.N+run.N2);
 allocdouble(&pos2_y,run.N+run.N2);
 allocdouble(&pos2_z,run.N+run.N2);
 
 compteur_deplac_court=0;
 compteur_deplac_long=0;
 }
 
 if(flag==1)
 {
   
  /*! calcul temps longs */
   
   itt1=step/dt_Fs;
   itt0=itt1% nbre_box;
   
   if(dt_Fs==(int) (run.tFs/nbre_box))
   //if(dt_Fs==(int) (run.nstep/nbre_box))
   {
   
#pragma omp parallel for default(shared) private (i) 
 for(i=0;i<run.N+run.N2;i++) //il faut enregistrer toutes les positions pour le calcul de Fs(k,t) !!!
{
 hist_pos_x_deplac_long[i][itt0]=pos2_cart[i].x;
 hist_pos_y_deplac_long[i][itt0]=pos2_cart[i].y;
 hist_pos_z_deplac_long[i][itt0]=pos2_cart[i].z;  
}
#pragma omp barrier

 itt1++; itt0=itt1% nbre_box;
 
 if(itt1>=nbre_box)
 {int intdeb=  nbre_box/ dt_Fs; compteur_deplac_long++;// compteur_deplac6_long++;
   for(itt2=intdeb;itt2<nbre_box;itt2++) // on ne recalcule pas aux temps courts
   {
     
     itt3=(itt0+itt2) % nbre_box;
     double hist_deplac=0;

     
//if(run.flagpara==1)
//{ 
#pragma omp parallel for default(shared) private (i)
  for(i=0;i<run.N+run.N2;i++) //il faut enregistrer toutes les positions pour le calcul de Fs(k,t) !!!
  {   
    pos1_x[i]=hist_pos_x_deplac_long[i][itt0];
    pos1_y[i]=hist_pos_y_deplac_long[i][itt0];
    pos1_z[i]=hist_pos_z_deplac_long[i][itt0];	
    pos2_x[i]=hist_pos_x_deplac_long[i][itt3];
    pos2_y[i]=hist_pos_y_deplac_long[i][itt3];
    pos2_z[i]=hist_pos_z_deplac_long[i][itt3];     
   } // fin boucle sur i
#pragma omp barrier   


   
#pragma omp parallel for default(shared) private (i,pos1,pos2)  reduction(+:hist_deplac) 
   for(i=0;i<run.N+run.N2;i++)
   {    
    pos1.x=pos1_x[i];
    pos1.y=pos1_y[i];
    pos1.z=pos1_z[i];
    pos2.x=pos2_x[i];
    pos2.y=pos2_y[i];
    pos2.z=pos2_z[i];
    double dist=distance_cart(pos1,pos2,arccos,run);
    hist_deplac+=dist;
   }
   
#pragma omp barrier
 hist_deplac_long[itt2]+=hist_deplac;
    } // fin boucle sur itt2 (t)
    
    

      
    
  } // fin if(itt1>=500)
   } // fin if(dt_Fs=...)
   
   

   
 } // fin if flag==1
 if(flag==2)
 {
  FILE *ffonction_deplacement;
  char ndeplacement[50];
  if(run.flagposanc==0)
  {

    sprintf(ndeplacement,"deplacNa%dNb%deta%dTemp%dbstep%dstep%dtFs%d.res",run.N,run.N2,(int) (1000*run.eta),(int) (1000/run.beta),(int) run.bstep,(int) run.nstep,(int) (run.tFs));
  }
  else  //sinon le nouveau bstep est bstep+step_oldsimu
  {
    sprintf(ndeplacement,"deplacNa%dNb%deta%dTemp%dbstep%dstep%dtFs%d.res",run.N,run.N2,(int) (1000*run.eta),(int) (1000/run.beta),(int) (run.bstep+run.bstep_old+run.step_old),(int) (run.nstep),(int) (run.tFs));
  }
    
    
    
  if ((ffonction_deplacement = fopen (ndeplacement,"w")) == NULL)
    { 
     fprintf(stderr,"Probleme d'ouverture de %s\n",ndeplacement);
     exit (1);  
    }
  
  
   int t;
   //int dt_Fs=(int) (run.tFs/nbre_box);
   int dt_Fs_court=1;
   int t_min=nbre_box*dt_Fs_court/dt_Fs;
  
  
   int t_deplac_prime=run.nstep;
   int dt_deplac=(int) (run.tFs/nbre_box);
  
   double deplac_0_long=(double) (compteur_deplac_long) *(run.N+run.N2)*run.sigma;
   double deplac_0_court=(double) (compteur_deplac_court) *(run.N+run.N2)*run.sigma;

    

   for(t=0;t<nbre_box;t++)
       {
	 fprintf(ffonction_deplacement,"%e\t %e\n",  t*run.dt/run.sigma*dt_deplac, hist_deplac_long[t]/deplac_0_long);
	 //fprintf(ffonction_deplacement,"%e\t %e\n",  t*run.dt/run.sigma*dt_Fs_court, hist_deplac_court[t]/((double) 500*run.N)); 
       }     
    fclose(ffonction_deplacement);
 
 } // 
 
  
} // fin calcul_deplac


// void calcul_deplac2(int flag, int step, int dt_Fs, POS_cart *pos2_cart,double *arccos, RUN run)
// {
//  int i;
//  int j;
//  int nbre_box=50; // nombre de points calculés
//  int dt_deplac=(int) (run.tFs/nbre_box);
//  //int itt1=step/dt_deplac;
//  int itt1=step/dt_Fs;
//  int itt0=itt1% nbre_box;
//  int itt2;
//  int itt3;
//  int tFs=run.tFs;
//  //int dt_Fs=tFs/500;
//  int dt_Fs_court=1;
//  
//  POS_cart pos1;
//  POS_cart pos2;
//  
//  static double **hist_pos_x_deplac_long=NULL; //histo des coordonnÃ©es phi (pour calcul des grandeurs dynamiques)
//  static double **hist_pos_y_deplac_long=NULL;
//  static double **hist_pos_z_deplac_long=NULL;
//  static double **hist_pos_x_deplac_court=NULL; //histo des coordonnÃ©es phi (pour calcul des grandeurs dynamiques)
//  static double **hist_pos_y_deplac_court=NULL;
//  static double **hist_pos_z_deplac_court=NULL;
//  
//  static double *pos1_x=NULL;
//  static double *pos1_y=NULL;
//  static double *pos1_z=NULL;
//  static double *pos2_x=NULL;
//  static double *pos2_y=NULL;
//  static double *pos2_z=NULL;
//  
//  static double *hist_deplac_long=NULL;
//  static double *hist_deplac_court=NULL;
//  static int compteur_deplac_court;
//  static int compteur_deplac_long;
// 
//  
//  if(flag==0)
//  {
//  allocdouble2(&hist_pos_x_deplac_long,run.N+run.N2,nbre_box);
//  allocdouble2(&hist_pos_y_deplac_long,run.N+run.N2,nbre_box);
//  allocdouble2(&hist_pos_z_deplac_long,run.N+run.N2,nbre_box);
//  allocdouble2(&hist_pos_x_deplac_court,run.N+run.N2,nbre_box);
//  allocdouble2(&hist_pos_y_deplac_court,run.N+run.N2,nbre_box);
//  allocdouble2(&hist_pos_z_deplac_court,run.N+run.N2,nbre_box);
//  allocdouble(&hist_deplac_long,nbre_box);
//  allocdouble(&hist_deplac_court,nbre_box);
//  
//  allocdouble(&pos1_x,run.N+run.N2);
//  allocdouble(&pos1_y,run.N+run.N2);
//  allocdouble(&pos1_z,run.N+run.N2);
//  allocdouble(&pos2_x,run.N+run.N2);
//  allocdouble(&pos2_y,run.N+run.N2);
//  allocdouble(&pos2_z,run.N+run.N2);
//  
//  compteur_deplac_court=0;
//  compteur_deplac_long=0;
//  }
//  
//  if(flag==1)
//  {
//    
//   //*********************** calcul temps longs 
//    
//    itt1=step/dt_Fs;
//    itt0=itt1% nbre_box;
//    
//    if(dt_Fs==(int) (run.tFs/nbre_box))
//    //if(dt_Fs==(int) (run.nstep/nbre_box))
//    {
//    
// #pragma omp parallel for default(shared) private (i) 
//  for(i=0;i<run.N+run.N2;i++) //il faut enregistrer toutes les positions pour le calcul de Fs(k,t) !!!
// {
//  hist_pos_x_deplac_long[i][itt0]=pos2_cart[i].x;
//  hist_pos_y_deplac_long[i][itt0]=pos2_cart[i].y;
//  hist_pos_z_deplac_long[i][itt0]=pos2_cart[i].z;  
// }
// #pragma omp barrier
// 
//  itt1++; itt0=itt1% nbre_box;
//  
//  if(itt1>=nbre_box)
//  {int intdeb=  nbre_box/ dt_Fs; compteur_deplac_long++;// compteur_deplac6_long++;
//    for(itt2=intdeb;itt2<nbre_box;itt2++) // on ne recalcule pas aux temps courts
//    {
//      
//      itt3=(itt0+itt2) % nbre_box;
//      double hist_deplac=0;
// 
//      
// //if(run.flagpara==1)
// //{ 
// #pragma omp parallel for default(shared) private (i)
//   for(i=0;i<run.N+run.N2;i++) //il faut enregistrer toutes les positions pour le calcul de Fs(k,t) !!!
//   {   
//     pos1_x[i]=hist_pos_x_deplac_long[i][itt0];
//     pos1_y[i]=hist_pos_y_deplac_long[i][itt0];
//     pos1_z[i]=hist_pos_z_deplac_long[i][itt0];	
//     pos2_x[i]=hist_pos_x_deplac_long[i][itt3];
//     pos2_y[i]=hist_pos_y_deplac_long[i][itt3];
//     pos2_z[i]=hist_pos_z_deplac_long[i][itt3];     
//    } // fin boucle sur i
// #pragma omp barrier   
// 
// 
//    
// #pragma omp parallel for default(shared) private (i,pos1,pos2)  reduction(+:hist_deplac) 
//    for(i=0;i<run.N+run.N2;i++)
//    {    
//     pos1.x=pos1_x[i];
//     pos1.y=pos1_y[i];
//     pos1.z=pos1_z[i];
//     pos2.x=pos2_x[i];
//     pos2.y=pos2_y[i];
//     pos2.z=pos2_z[i];
//     double dist=distance_cart(pos1,pos2,arccos,run);
//     hist_deplac+=dist*dist;
//    }
//    
// #pragma omp barrier
//  hist_deplac_long[itt2]+=hist_deplac;
//     } // fin boucle sur itt2 (t)
//     
//     
// 
//       
//     
//   } // fin if(itt1>=500)
//    } // fin if(dt_Fs=...)
//    
//    
// 
//    
//  } // fin if flag==1
//  if(flag==2)
//  {
//   FILE *ffonction_deplacement;
//   char ndeplacement[50];
//   if(run.flagposanc==0)
//   {
// 
//     sprintf(ndeplacement,"deplac2Na%dNb%deta%dTemp%dbstep%dstep%dtFs%d.res",run.N,run.N2,(int) (1000*run.eta),(int) (1000/run.beta),(int) run.bstep,(int) run.nstep,(int) (run.tFs));
//   }
//   else  //sinon le nouveau bstep est bstep+step_oldsimu
//   {
//     sprintf(ndeplacement,"deplac2Na%dNb%deta%dTemp%dbstep%dstep%dtFs%d.res",run.N,run.N2,(int) (1000*run.eta),(int) (1000/run.beta),(int) (run.bstep+run.bstep_old+run.step_old),(int) (run.nstep),(int) (run.tFs));
//   }
//     
//     
//     
//   if ((ffonction_deplacement = fopen (ndeplacement,"w")) == NULL)
//     { 
//      fprintf(stderr,"Probleme d'ouverture de %s\n",ndeplacement);
//      exit (1);  
//     }
//   
//   
//    int t;
//    //int dt_Fs=(int) (run.tFs/nbre_box);
//    int dt_Fs_court=1;
//    int t_min=nbre_box*dt_Fs_court/dt_Fs;
//    //int t_min=nbre_box/10;
//   
//    int t_deplac_prime=run.nstep;
//    int dt_deplac=(int) (run.tFs/nbre_box);
//   
//    double deplac_0_long=(double) (compteur_deplac_long) *(run.N+run.N2)*run.sigma*run.sigma;
//    double deplac_0_court=(double) (compteur_deplac_court) *(run.N+run.N2)*run.sigma*run.sigma;
// 
//     
// 
// //    for(t=0;t<nbre_box;t++)
// //        {
// // 	 fprintf(ffonction_deplacement,"%e\t %e\n",  t*run.dt/run.sigma*dt_deplac, hist_deplac_long[t]/deplac_0_long);
// // 	 //fprintf(ffonction_deplacement,"%e\t %e\n",  t*run.dt/run.sigma*dt_Fs_court, hist_deplac_court[t]/((double) 500*run.N)); 
// //        }     
//    //*****************
//    
//    
//       for(t=0;t<nbre_box;t++)
//       {
//        fprintf(ffonction_deplacement,"%e\t %e\n",  t*run.dt/run.sigma*1, hist_deplac_court[t]/deplac_0_court); 
//       } 
//        
//       for(t=t_min;t<nbre_box;t++)
//       {
//        fprintf(ffonction_deplacement,"%e\t %e\n",  t*run.dt/run.sigma*dt_deplac, hist_deplac_long[t]/deplac_0_long);
//       } 
//    
//     fclose(ffonction_deplacement);
//  
//  } // 
//  
//   
// } // fin calcul_deplac2


/*! r6, r5, r7, r66 */
void calcul_deplac657(int flag, int step, int dt_Fs, POS_cart *pos2_cart, double *arccos, THERMO *thermo, RUN run)
{
 int i;
 int j;
 int nbre_box=50; // nombre de points calculés
 int nbre_pts=50; // nombre de points sur lesquels on moyenne
 int dt_deplac=(int) (run.tFs/nbre_pts);
 //int dt_deplac=(int) (run.nstep/nbre_pts);
 int itt1=step/dt_deplac; // t'
 int itt0=itt1% nbre_pts; // t' modulo nbre_box
 int itt2; // t 
 int itt3; // t'+t modulo nbre_box
 int itt4; // t''
 int itt5; //t'' modulo nbre_box
 int tFs=run.tFs;
 //int dt_Fs=tFs/500;
 int dt_Fs_court=1;
 
 POS_cart pos1;
 POS_cart pos2;
 
 static double **hist_pos_x_deplac_long=NULL; //histo des coordonnÃ©es phi (pour calcul des grandeurs dynamiques)
 static double **hist_pos_y_deplac_long=NULL;
 static double **hist_pos_z_deplac_long=NULL;
 static double **hist_pos_x_deplac_court=NULL; //histo des coordonnÃ©es phi (pour calcul des grandeurs dynamiques)
 static double **hist_pos_y_deplac_court=NULL;
 static double **hist_pos_z_deplac_court=NULL;
 
 static double *pos1_x=NULL; // pos à t'
 static double *pos1_y=NULL; 
 static double *pos1_z=NULL;
 static double *pos2_x=NULL; // pos à t'+t
 static double *pos2_y=NULL;
 static double *pos2_z=NULL;
 static double *pos3_x=NULL; // pos à t'', pour r6, r7, r5
 static double *pos3_y=NULL;
 static double *pos3_z=NULL;
 
 static int *col=NULL; // tableau des couleurs des particules (coordinences)
 static int *ndislo=NULL; // pointeur sur un  entier donnant le nombre de dislocations pour une config
 static int *voisins=NULL;
 static int *denom=NULL; // tableau des Delta_i(t,t')
 
 static double *hist_deplac6_long=NULL;
 static double *hist_deplac66_long=NULL; // deplacement des particules à 6 voisins et dont les voisins ont aussi 6 voisins (pour être suffisamment loin des 5 et 7)
 static double *hist_deplac5_long=NULL;
 static double *hist_deplac7_long=NULL;
 
 static int compteur_deplac_court;
 static int compteur_deplac_long;
 
 if(flag==0)
 {
 allocdouble2(&hist_pos_x_deplac_long,run.N+run.N2,nbre_box);
 allocdouble2(&hist_pos_y_deplac_long,run.N+run.N2,nbre_box);
 allocdouble2(&hist_pos_z_deplac_long,run.N+run.N2,nbre_box);
 allocdouble2(&hist_pos_x_deplac_court,run.N+run.N2,nbre_box);
 allocdouble2(&hist_pos_y_deplac_court,run.N+run.N2,nbre_box);
 allocdouble2(&hist_pos_z_deplac_court,run.N+run.N2,nbre_box);
 allocdouble(&hist_deplac6_long,nbre_box);
 allocdouble(&hist_deplac66_long,nbre_box);
 allocdouble(&hist_deplac5_long,nbre_box);
 allocdouble(&hist_deplac7_long,nbre_box);
 
 allocdouble(&pos1_x,run.N+run.N2);
 allocdouble(&pos1_y,run.N+run.N2);
 allocdouble(&pos1_z,run.N+run.N2);
 allocdouble(&pos2_x,run.N+run.N2);
 allocdouble(&pos2_y,run.N+run.N2);
 allocdouble(&pos2_z,run.N+run.N2);
 allocdouble(&pos3_x,run.N+run.N2);
 allocdouble(&pos3_y,run.N+run.N2);
 allocdouble(&pos3_z,run.N+run.N2);
 allocint(&col,run.N+run.N2);
 allocint(&ndislo,1);
 allocint(&voisins,(run.N+run.N2)*8);
 allocint(&denom,run.N+run.N2);
 compteur_deplac_court=0;
 compteur_deplac_long=0;
 }
 
 if(flag==1)
 {
   
   //********************** on réinitialise denom[i]
// #pragma omp parallel for default(shared) private (i)
//    for(i=0;i<run.N+run.N2;i++)
//    {
//     denom[i]=0; 
//    }
// #pragma omp barrier
   
   
   
  //*********************** calcul temps longs 
   
   itt1=step/dt_Fs;
   itt0=itt1% nbre_pts;
   
   if(dt_Fs==(int) (run.tFs/nbre_pts))
   //if(dt_Fs==(int) (run.nstep/nbre_pts))
   {
   
#pragma omp parallel for default(shared) private (i) 
 for(i=0;i<run.N+run.N2;i++) /*!il faut enregistrer toutes les positions pour le calcul de Fs(k,t) !*/
 {
  hist_pos_x_deplac_long[i][itt0]=pos2_cart[i].x;
  hist_pos_y_deplac_long[i][itt0]=pos2_cart[i].y;
  hist_pos_z_deplac_long[i][itt0]=pos2_cart[i].z;  
 }
#pragma omp barrier

 itt1++; itt0=itt1% nbre_pts;
 
  // Delta(t,t') par lequel on normalise delta_6,i(t'')
 
if(itt1>=nbre_pts)
 {int intdeb=nbre_pts/dt_Fs;
 compteur_deplac_long++;// compteur_deplac6_long++;
   

   
   //**************** Calcul des déplacements moyens
   
   for(itt2=1;itt2<nbre_box;itt2++) /*! boucle sur t */
   {
     
     itt3=(itt0+itt2) % nbre_box;
     double hist_deplac=0;
     double hist_deplac6=0;
     double hist_deplac5=0;
     double hist_deplac7=0;
     double hist_deplac66=0;
     
        //*************** Calcul des Dénom Delta_i(t,t')
  // itt3=(itt0+itt2) % 100;
//    for(itt4=itt0;itt4<itt0+itt2;itt4++) // boucle sur t'',t'<=t''<=t'+t, on ne recalcule pas aux temps courts
//    {
//      itt5=itt4 % nbre_box;
//      
// #pragma omp parallel for default(shared) private (i)
//      for(i=0;i<run.N+run.N2;i++)
//      {
//       pos2_x[i]=hist_pos_x_deplac_long[i][itt5];
//       pos2_y[i]=hist_pos_y_deplac_long[i][itt5];
//       pos2_z[i]=hist_pos_z_deplac_long[i][itt5];   
//      }
// #pragma omp barrier 
//      for(i=0;i<run.N+run.N2;i++)
//      {
//       for(j=0;j<8;j++)
//       {
//        voisins[8*i+j]=0;
//       }
//      }
//    voronoi(pos2_x, pos2_y, pos2_z,voisins,col,run.N+run.N2,ndislo);
//    
//    int denom_i;
// #pragma omp parallel for default(shared) private (i)  reduction(+:denom_i) 
//    for(i=0;i<run.N+run.N2;i++)
//    {  
//     denom_i=0;
//     int delta6=0;
//     int delta5=0;
//     int delta7=0;
//     int delta66=1;
//     
//     if(col[i]==6){delta6=1;}
//     denom_i+=delta6;
//     if(col[i]==5){delta5=1;}
//     denom_i+=delta5;
//     if(col[i]==7){delta7=1;} 
//     denom_i+=delta7;
//     denom[i]+=denom_i;
//    }
// #pragma omp barrier
//   // printf("Delta_i=1(t,t') %d\n", denom[0]);
//    
//    } // fin for itt4
   //***************** Fin calcul des Dénom Delta_i(t,t')
 
  
   //***************** Calcul de r6, r5, r7, r66(t)
  
#pragma omp parallel for default(shared) private (i)
  for(i=0;i<run.N+run.N2;i++) //il faut enregistrer toutes les positions pour le calcul de Fs(k,t) !!!
  {   
    pos1_x[i]=hist_pos_x_deplac_long[i][itt0]; // pos à t'
    pos1_y[i]=hist_pos_y_deplac_long[i][itt0];
    pos1_z[i]=hist_pos_z_deplac_long[i][itt0];	
    pos2_x[i]=hist_pos_x_deplac_long[i][itt3]; // pos à t'+t
    pos2_y[i]=hist_pos_y_deplac_long[i][itt3];
    pos2_z[i]=hist_pos_z_deplac_long[i][itt3]; 
  }
#pragma omp barrier
  
  for(itt4=itt0;itt4<itt0+itt2;itt4++) // boucle sur t''
  {
  
  itt5=itt4 % nbre_box;
  
#pragma omp parallel for default(shared) private (i)
  for(i=0;i<run.N+run.N2;i++) //il faut enregistrer toutes les positions pour le calcul de Fs(k,t) !!!
  {  
    pos3_x[i]=hist_pos_x_deplac_long[i][itt5]; // pos à t'', calculée en plus de celles à t' et t'+t pour la pondération
    pos3_y[i]=hist_pos_y_deplac_long[i][itt5];
    pos3_z[i]=hist_pos_z_deplac_long[i][itt5];    
  } // fin boucle sur i
#pragma omp barrier

   for(i=0;i<run.N+run.N2;i++)
   {
    for(j=0;j<8;j++)
    {
     voisins[8*i+j]=0;
    }
   }
   voronoi(pos3_x, pos3_y, pos3_z,voisins,col,run.N+run.N2,ndislo);
  
  
#pragma omp parallel for default(shared) private (i,pos1,pos2)  reduction(+:hist_deplac6, hist_deplac5,  hist_deplac7, hist_deplac66) 
   for(i=0;i<run.N+run.N2;i++)
   {    
    pos1.x=pos1_x[i];
    pos1.y=pos1_y[i];
    pos1.z=pos1_z[i];
    pos2.x=pos2_x[i];
    pos2.y=pos2_y[i];
    pos2.z=pos2_z[i];
    double dist=distance_cart(pos1,pos2,arccos,run); // dist(t',t'+t)
    double delta6=0;
    double delta5=0;
    double delta7=0;
    double delta66=1; // on rend delta66 nul dès qu'un voisin a un nbre de voisins différent de 6

    
    /*! nouvelle def de r6, r5, r7, r66 */
    
    if(col[i]==6){delta6=1;} 
    if(col[i]==5){delta5=1;}
    if(col[i]==7){delta7=1;}
    if(col[i]==6)
    {
     for(j=0;j<8;j++)
     {
      int vois=voisins[8*i+j];
     // printf("i %d\t j %d\t vois %d\n",i,j,vois);
      if((vois>0)&&(col[vois]!=6)){delta66=0;} // on teste si vois est bien l'indice d'un voisin (sinon vois=0) et si vois n'a pas 6 voisins
     }
    }
    
    //if(denom[i]!=0)
    //{
//      hist_deplac6+=dist*delta6/( (double) denom[i]);
//      hist_deplac5+=dist*delta5/( (double) denom[i]);
//      hist_deplac7+=dist*delta7/( (double) denom[i]);
//      hist_deplac66+=dist*delta66/( (double) denom[i]);
     hist_deplac6+=dist*delta6;
     hist_deplac5+=dist*delta5;
     hist_deplac7+=dist*delta7;
     hist_deplac66+=dist*delta66;
    //}
   } // fin boucle sur i   
#pragma omp barrier
 
 
 
     } // fin boucle sur itt4 (t'')
     hist_deplac6_long[itt2]+=((double) hist_deplac6)/itt2; // on divise par t (pondération : 1/t*sum_{t''=0,t}\delta_n,i(t''))
 //printf("hist_deplac6 %e\n",hist_deplac6_long[itt2]);
 hist_deplac5_long[itt2]+=((double) hist_deplac5)/itt2;
 hist_deplac7_long[itt2]+=((double) hist_deplac7)/itt2;
 hist_deplac66_long[itt2]+=((double) hist_deplac66)/itt2;
    } // fin boucle sur itt2 (t)

      
    
  } // fin if(itt1>=500)
   } // fin if(dt_Fs=...)
   
   

   
 } // fin if flag==1
 
 
 if(flag==2)
 {

 /*! Enregistrement dans un fichier <r6(t)> */   
    FILE *ffonction_deplacement657;
  char ndeplacement657[50];
  if(run.flagposanc==0)
  {

    sprintf(ndeplacement657,"deplac657Na%dNb%deta%dTemp%dbstep%dstep%dtFs%d.res",run.N,run.N2,(int) (1000*run.eta),(int) (1000/run.beta),(int) run.bstep,(int) run.nstep,(int) (run.tFs));
  }
  else  //sinon le nouveau bstep est bstep+step_oldsimu
  {
    sprintf(ndeplacement657,"deplac657Na%dNb%deta%dTemp%dbstep%dstep%dtFs%d.res",run.N,run.N2,(int) (1000*run.eta),(int) (1000/run.beta),(int) (run.bstep+run.bstep_old+run.step_old),(int) (run.nstep),(int) (run.tFs));
  }
    
    
    
  if ((ffonction_deplacement657 = fopen (ndeplacement657,"w")) == NULL)
    { 
     fprintf(stderr,"Probleme d'ouverture de %s\n",ndeplacement657);
     exit (1);  
    }
  
  
   int t;
   int dt_Fs=(int) (run.tFs/nbre_pts);
  //int dt_Fs=(int) (run.nstep/nbre_pts);
   int dt_Fs_court=1;
   int t_min=nbre_pts*dt_Fs_court/dt_Fs;
  
  
   int t_deplac_prime=run.nstep;
   int dt_deplac=(int) (run.tFs/nbre_pts); // pour définir le temps réel il
   
 double n6=thermo->histovoro[6]/((double) run.nstep)*run.moviestep;
 double n5=thermo->histovoro[5]/((double) run.nstep)*run.moviestep;
 double n7=thermo->histovoro[7]/((double) run.nstep)*run.moviestep;
 double n66=thermo->n66[0]/((double) run.nstep)*run.moviestep;
 
 printf("n6 %e\n",n6);
 printf("n5 %e\n",n5);
 printf("n7 %e\n",n7);
 printf("n66 %e\n",n66);
 
  double deplac6_0_long=(double) (compteur_deplac_long)*n6*run.sigma; 
  double deplac5_0_long=(double) (compteur_deplac_long)*n5*run.sigma;   
  double deplac7_0_long=(double) (compteur_deplac_long)*n7*run.sigma;
  double deplac66_0_long=(double) (compteur_deplac_long)*n66*run.sigma; 
    
//     for(t=0;t<100;t++)
//        {
// 	 fprintf(ffonction_deplacement657,"%e\t %e\t %e\t %e\t %e\n",  t*run.dt/run.sigma*dt_deplac, hist_deplac6_long[t]/deplac6_0_long,
// 	   hist_deplac5_long[t]/deplac5_0_long, hist_deplac7_long[t]/deplac7_0_long, hist_deplac66_long[t]/deplac66_0_long);
// 	 //fprintf(ffonction_deplacement,"%e\t %e\n",  t*run.dt/run.sigma*dt_Fs_court, hist_deplac_court[t]/((double) 500*run.N)); 
//        } // ancienne def
//   
     
       for(t=0;t<nbre_box;t++)
       {
	 fprintf(ffonction_deplacement657,"%e\t %e\t %e\t  %e\t %e\n",  t*run.dt/run.sigma*dt_deplac, hist_deplac6_long[t]/deplac6_0_long,
	   hist_deplac5_long[t]/deplac5_0_long, hist_deplac7_long[t]/deplac7_0_long, hist_deplac66_long[t]/deplac66_0_long);
	 //fprintf(ffonction_deplacement,"%e\t %e\n",  t*run.dt/run.sigma*dt_Fs_court, hist_deplac_court[t]/((double) 500*run.N)); 
       } // ancienne def
  
   
    fclose(ffonction_deplacement657);
   
 }
 
  
} // fin calcul_deplac657

/*! calcul des trajectoires de 20 particules + nombre de voisins */

void calcul_dist20part(int flag, int step, POS_cart *pos2_cart, double *arccos, RUN run)
{
 
 int i;
// int dt_deplac=(int) (run.tFs/1000);
 int itt1;
 int dt_simu=(int) (run.nstep/500);
 int n_part_suivies=20;
 
 /*static double *pos_x_part_suivie;
 static double *pos_y_part_suivie;
 static double *pos_z_part_suivie;*/
 static double **hist_deplac_part_suivies;
 static POS_cart *pos_part_suivie; // tableau des positions des particules que l'on souhaite suivre
 static POS_cart *pos_part_suivie_ini; // tableau des positions initiales associÃ©es    
 static double **nvoisins;
 static int *col=NULL; // tableau des couleurs des particules (coordinences)
 static int *ndislo=NULL; // pointeur sur un  entier donnant le nombre de dislocations pour une config
 static int *voisins=NULL;
 static double *pos_x;
 static double *pos_y;
 static double *pos_z;
 
 if(flag==0)
 {
/*  allocdouble(&pos_x_part_suivie,run.N);
  allocdouble(&pos_y_part_suivie,run.N);
  allocdouble(&pos_z_part_suivie,run.N);*/
  allocint(&col,run.N+run.N2);
  allocint(&ndislo,1);
  allocdouble(&pos_x,run.N+run.N2);
  allocdouble(&pos_y,run.N+run.N2);
  allocdouble(&pos_z,run.N+run.N2);
  allocint(&voisins,(run.N+run.N2)*8);
  allocdouble2(&hist_deplac_part_suivies,n_part_suivies,500);
  allocdouble2(&nvoisins,n_part_suivies,500);
  pos_part_suivie=(POS_cart *) malloc(sizeof(POS_cart)*n_part_suivies);
  pos_part_suivie_ini=(POS_cart *) malloc(sizeof(POS_cart)*n_part_suivies);
  
  for(i=0;i<n_part_suivies;i++)
  {     
   pos_part_suivie_ini[i]=pos2_cart[i];	  
  }
 }
 
 if(flag==1)
 {
  itt1=step/dt_simu; //itt1=0,1,...,nstep/dt_Fs      
  if((itt1<500)&&(itt1>=0))
  {	 
   for(i=0;i<n_part_suivies;i++)
   {
    pos_part_suivie[i]=pos2_cart[i];
    double dist=distance_cart(pos_part_suivie_ini[i],pos_part_suivie[i],arccos,run); 
    hist_deplac_part_suivies[i][itt1]=dist;
   }
  }
  
  #pragma omp parallel for private(i)
  for(i=0;i<run.N+run.N2;i++)
  {
   pos_x[i]=pos2_cart[i].x;
   pos_y[i]=pos2_cart[i].y;
   pos_z[i]=pos2_cart[i].z;
  }
  
  voronoi(pos_x, pos_y, pos_z,voisins,col,run.N+run.N2,ndislo);
  for(i=0;i<n_part_suivies;i++)
  {
   nvoisins[i][itt1]=col[i];
  }
  
  }// fin if flag == 1
 
 
 if(flag==2)
 {
  FILE *fdist;
	char ndist[50];
        if(run.flagposanc==0)
	{
	 sprintf(ndist,"dist20partNa%dNb%deta%dTemp%dbstep%dstep%d.res",run.N,run.N2,(int) (1000*run.eta),(int) (1000/run.beta),(int) run.bstep, (int) run.nstep);
	}
	else
        {
	 sprintf(ndist,"dist20partNa%dNb%deta%dTemp%dbstep%dstep%d.res",run.N,run.N2,(int) (1000*run.eta),(int) (1000/run.beta),(int) run.bstep+ (int) run.bstep_old+ (int) run.step_old,(int) run.nstep);
	}
        if ((fdist = fopen (ndist,"w")) == NULL)
        { fprintf(stderr,"Probleme d'ouverture de %s\n",ndist);
          exit (1);  
        }
	  
	  
	 int t;
         for(t=0;t<500;t++)
	 {
	  fprintf(fdist, "%e\t",t*run.dt/run.sigma*dt_simu);
	  for(i=0;i<n_part_suivies;i++)
	  {
	  fprintf(fdist,"%e\t", hist_deplac_part_suivies[i][t]/(run.sigma));
	  }
	  fprintf(fdist,"\n");
	 }
	 fclose(fdist);
	 
  FILE *fdefauts;
	char ndefauts[50];
        if(run.flagposanc==0)
	{
	 sprintf(ndefauts,"ndefauts20partNa%dNb%deta%dTemp%dbstep%dstep%d.res",run.N,run.N2,(int) (1000*run.eta),(int) (1000/run.beta),(int) run.bstep, (int) run.nstep);
	}
	else
        {
	 sprintf(ndefauts,"ndefauts20partNa%dNb%deta%dTemp%dbstep%dstep%d.res",run.N,run.N2,(int) (1000*run.eta),(int) (1000/run.beta),(int) run.bstep+ (int) run.bstep_old+ (int) run.step_old,(int) run.nstep);
	}
        if ((fdefauts = fopen (ndefauts,"w")) == NULL)
        { fprintf(stderr,"Probleme d'ouverture de %s\n",ndefauts);
          exit (1);  
        }
        
        for(t=0;t<500;t++)
	 {
	  fprintf(fdefauts, "%e\t",t*run.dt/run.sigma*dt_simu);
	  for(i=0;i<n_part_suivies;i++)
	  {
	  fprintf(fdefauts,"%e\t", nvoisins[i][t]);
	  }
	  fprintf(fdefauts,"\n");
	 }
	 fclose(fdefauts);
 } // fin if flag == 2
  
} // fin calcul_dist20part


///********************************************

/*! Fs(k,t) */
void calcul_Fs(int flag, int step, int dt_Fs, POS_cart *pos2_cart, double *Fs, double *arccos, double *legendre, RUN run)
{
 int i;
 //int dt_Fs=(int) (run.tFs/500);
 int nbre_pts=100;
 int itt1=step/dt_Fs;
 int itt0=itt1% nbre_pts;
 int itt2;
 int itt3;
 static double **hist_pos_x=NULL; //histo des coordonnees (pour calcul des grandeurs dynamiques)
 static double **hist_pos_y=NULL;
 static double **hist_pos_z=NULL;
 
 static double **hist_pos_x_court=NULL; //histo des coordonnees phi (pour calcul des grandeurs dynamiques)
 static double **hist_pos_y_court=NULL;
 static double **hist_pos_z_court=NULL; 
 static double *hist_Fs_court=NULL;  // temps courts
 static double *hist_Fs_long=NULL;
 static double *hist_deplacement_court=NULL;
 static double *hist_deplacement_long=NULL;
 
 static double *hist_Fs_court_2=NULL;  //temps courts pour les particules de type 2
 static double *hist_Fs_long_2=NULL;
 static double *hist_deplacement_court_2=NULL;
 static double *hist_deplacement_long_2=NULL;

 
 double startwtime,endwtime;
 static int compteurFs_court;
 static int compteurFs_long;
 POS_cart pos1;
 POS_cart pos2;
 
 if(flag==0)
 {
 allocdouble2(&hist_pos_x,run.N+run.N2,nbre_pts);
 allocdouble2(&hist_pos_y,run.N+run.N2,nbre_pts);
 allocdouble2(&hist_pos_z,run.N+run.N2,nbre_pts);  
 allocdouble2(&hist_pos_x_court,run.N+run.N2,nbre_pts);
 allocdouble2(&hist_pos_y_court,run.N+run.N2,nbre_pts);
 allocdouble2(&hist_pos_z_court,run.N+run.N2,nbre_pts);
 allocdouble(&hist_Fs_long,nbre_pts);
 allocdouble(&hist_Fs_court,nbre_pts);
 allocdouble(&hist_deplacement_court,nbre_pts);
 allocdouble(&hist_deplacement_long,nbre_pts);
 if(run.N2>0)
 {
  allocdouble(&hist_Fs_long_2,nbre_pts);
  allocdouble(&hist_Fs_court_2,nbre_pts);
  allocdouble(&hist_deplacement_court_2,nbre_pts);
  allocdouble(&hist_deplacement_long_2,nbre_pts);
 }
 
 
 int i;
 int j; 
 compteurFs_court=0;
 compteurFs_long=0;

//  for(i=0;i<nbre_pts;i++)
//  {
//    hist_Fs_long[i]=0;
//    hist_Fs_court[i]=0;
//    hist_Fs_long_2[i]=0;
//    hist_Fs_court_2[i]=0;
//    hist_deplacement_court[i]=0;
//    hist_deplacement_long[i]=0;
//    hist_deplacement_court_2[i]=0;
//    hist_deplacement_long_2[i]=0;
//  }
 } 
 
 if(flag==1)
 {
   itt1=step/dt_Fs;
   itt0=itt1% nbre_pts;

   if(dt_Fs==(int) (run.tFs/nbre_pts))
   {
#pragma omp parallel for default(shared)   private(i)   
  for(i=0;i<run.N+run.N2;i++) //il faut enregistrer toutes les positions pour le calcul de Fs(k,t) !!!
 {
  hist_pos_x[i][itt0]=pos2_cart[i].x;
  hist_pos_y[i][itt0]=pos2_cart[i].y;
  hist_pos_z[i][itt0]=pos2_cart[i].z;
 }
 
 itt1++;  itt0=itt1% nbre_pts;
   if (itt1>=nbre_pts)
       {
	 int intdeb=  nbre_pts/ dt_Fs;
	 compteurFs_long++;
  for (itt2=intdeb;itt2<nbre_pts;itt2++)
	 {   
	      itt3=(itt0+itt2) % nbre_pts;
	      double hist_cur=0;
	      double hist_cur_tot=0;
	      double hist_deplac=0;
	      
	   //************** Calcul Fs_11
#pragma omp parallel for default(shared) private (i,pos1,pos2) reduction(+:hist_cur) reduction(+:hist_deplac)
  for(i=0;i<run.N;i++) //il faut enregistrer toutes les positions pour le calcul de Fs(k,t) !!!
 	 {	 
           pos1.x=hist_pos_x[i][itt0];
	   pos1.y=hist_pos_y[i][itt0];
           pos1.z=hist_pos_z[i][itt0];	
	   pos2.x=hist_pos_x[i][itt3];
           pos2.y=hist_pos_y[i][itt3];
           pos2.z=hist_pos_z[i][itt3];   
	   double dist=distance_cart(pos1,pos2,arccos,run);
	   hist_cur+=calcul_legendre(dist,legendre,run); 
	   hist_deplac+=dist;
	 } //boucle sur i
#pragma omp barrier
	 hist_Fs_long[itt2]+=hist_cur;
	 hist_deplacement_long[itt2]+=hist_deplac;
	 
// if(run.flagFtot==1)
// {
//   int j;
// #pragma omp parallel for default(shared) private (i,pos1,pos2) reduction(+:hist_cur_tot)
//  for(i=0;i<run.N;i++)
//  {
//   for(j=0;j<run.N;j++)
//   {
//    pos1.x=hist_pos_x[i][itt0];
//    pos1.y=hist_pos_y[i][itt0];
//    pos1.z=hist_pos_z[i][itt0];	
//    pos2.x=hist_pos_x[j][itt3];
//    pos2.y=hist_pos_y[j][itt3];
//    pos2.z=hist_pos_z[j][itt3];
//    double dist=distance_cart(pos1,pos2,arccos,run);
//    hist_cur_tot+=calcul_legendre(dist,legendre,run);
//   }
//  }   
// #pragma omp barrier
//  hist_Ftot_long[itt2]+=hist_cur_tot;
// } // fin if flagtot=1

	  
	 }
	 
	 
	  //************* calcul Fs_22
	  if(run.N2>0)
	  {
	    for (itt2=intdeb;itt2<nbre_pts;itt2++)
	 {   
	      itt3=(itt0+itt2) % nbre_pts;
	    
	    
	  double hist_cur_2=0;
	  double hist_deplac_2=0;
#pragma omp parallel for default(shared) private (i,pos1,pos2) reduction(+:hist_cur_2) reduction(+:hist_deplac_2)
  for(i=run.N;i<run.N+run.N2;i++) //il faut enregistrer toutes les positions pour le calcul de Fs(k,t) !!!
 	 {	 
           pos1.x=hist_pos_x[i][itt0];
	   pos1.y=hist_pos_y[i][itt0];
           pos1.z=hist_pos_z[i][itt0];	
	   pos2.x=hist_pos_x[i][itt3];
           pos2.y=hist_pos_y[i][itt3];
           pos2.z=hist_pos_z[i][itt3];   
	   double dist=distance_cart(pos1,pos2,arccos,run);
	   hist_cur_2+=calcul_legendre(dist,legendre,run); 
	   hist_deplac_2+=dist;
	 } //boucle sur i
#pragma omp barrier
	  hist_Fs_long_2[itt2]+=hist_cur_2;
	  hist_deplacement_long_2[itt2]+=hist_deplac_2;
	  }
	  
	 } // fin boucle sur itt2
       } // fin if(itt1>=nbre_pts)
       

       
    //   if (omp_get_thread_num() == 0)    endwtime =  omp_get_wtime(); 
//if (omp_get_thread_num() == 0) printf("calcul Fs time  = %e\n",(double)(endwtime-startwtime));
   } // fin si if dt_Fs=tFs/nbre_pts
   if(dt_Fs==1)
   {
      if(step<10000) // on ne calcule le Fs aux temps courts que sur un intervalle de temps petit devant le temps total de la simu
      {      
       itt1=step/dt_Fs;
       itt0=itt1% nbre_pts; // itt0=0,1,...499    
#pragma omp parallel for default(shared) private (i)         
       for(i=0;i<run.N+run.N2;i++) //il faut enregistrer toutes les positions pour le calcul de Fs(k,t) !!!
	 {
	  hist_pos_x_court[i][itt0]=pos2_cart[i].x;
          hist_pos_y_court[i][itt0]=pos2_cart[i].y;
          hist_pos_z_court[i][itt0]=pos2_cart[i].z;
	 }
	 
	itt1++;  itt0=itt1% nbre_pts;    
      if (itt1>=nbre_pts)
       {
        compteurFs_court++;	 
	 for (itt2=0;itt2<nbre_pts;itt2++)
	 {  
	   
	   itt3=(itt0+itt2) % nbre_pts;
	   
	   //**************** calcul Fs_11
	   double hist_cur=0;
	   double hist_cur_tot=0;
	   double hist_deplac=0;
	   
#pragma omp parallel for default(shared) private (i,pos1,pos2) reduction(+:hist_cur) reduction(+:hist_deplac)
           for(i=0;i<run.N;i++) //il faut enregistrer toutes les positions pour le calcul de Fs(k,t) !!!
	   {      
	     pos1.x=hist_pos_x_court[i][itt0];
	     pos1.y=hist_pos_y_court[i][itt0];
             pos1.z=hist_pos_z_court[i][itt0];
	     pos2.x=hist_pos_x_court[i][itt3];
             pos2.y=hist_pos_y_court[i][itt3];
             pos2.z=hist_pos_z_court[i][itt3];      
	     double dist=distance_cart(pos1,pos2,arccos,run);
	     hist_cur+=calcul_legendre(dist,legendre,run);
	     hist_deplac+=dist;
	   }  //fin boucle i	 
#pragma omp barrier  	 
	 hist_Fs_court[itt2]+=hist_cur;
	 hist_deplacement_court[itt2]+=hist_deplac;
	 

	 
	 } // fin for itt2
	 
	 //************ calcul Fs_22
	 if(run.N2>0)
	  {
	    for (itt2=0;itt2<nbre_pts;itt2++)
	 {  
	   
	   itt3=(itt0+itt2) % nbre_pts;
	  double hist_cur_2=0;
	  double hist_deplac_2=0;
#pragma omp parallel for default(shared) private (i,pos1,pos2) reduction(+:hist_cur_2) reduction(+:hist_deplac_2)
  for(i=run.N;i<run.N+run.N2;i++) //il faut enregistrer toutes les positions pour le calcul de Fs(k,t) !!!
 	 {	 
           pos1.x=hist_pos_x_court[i][itt0];
	   pos1.y=hist_pos_y_court[i][itt0];
           pos1.z=hist_pos_z_court[i][itt0];	
	   pos2.x=hist_pos_x_court[i][itt3];
           pos2.y=hist_pos_y_court[i][itt3];
           pos2.z=hist_pos_z_court[i][itt3];   
	   double dist=distance_cart(pos1,pos2,arccos,run);
	   hist_cur_2+=calcul_legendre(dist,legendre,run); 
	   hist_deplac_2+=dist;
	 } //boucle sur i
#pragma omp barrier
	  hist_Fs_court_2[itt2]+=hist_cur_2;
	  hist_deplacement_court_2[itt2]+=hist_deplac_2;
	  } // fin boucle itt2
	 
	 
	 } //fin if N2>0
       } // fin if itt1
 
      
      } // fin du if(step<20000)
   // if (omp_get_thread_num() == 0)    endwtime =  omp_get_wtime(); 
//if (omp_get_thread_num() == 0) printf("calcul Fs  2 time  = %e\n",(double)(endwtime-startwtime));    
 }
// if (omp_get_thread_num() == 0)    endwtime =  omp_get_wtime(); 
// if (omp_get_thread_num() == 0) printf("calcul Fs  2 time  = %e\n",(double)(endwtime-startwtime)); 
//  
 } //fin du if (flag == 1)
 if(flag==2)
 {
  FILE *ffonction_Fs;
  char nFs[50];
  if(run.flagposanc==0)
  {
    sprintf(nFs,"FsaaNa%dNb%deta%dTemp%dbstep%dstep%dtFs%d.res",run.N,run.N2,(int) (1000*run.eta),(int) (1000/run.beta),(int) run.bstep,(int) run.nstep,(int) (run.tFs));
  }
  else  //sinon le nouveau bstep est bstep+step_oldsimu
  {
    sprintf(nFs,"FsaaNa%dNb%deta%dTemp%dbstep%dstep%dtFs%d.res",run.N,run.N2,(int) (1000*run.eta),(int) (1000/run.beta),(int) (run.bstep+run.bstep_old+run.step_old),(int) (run.nstep),(int) (run.tFs));
  }
            
  if ((ffonction_Fs = fopen (nFs,"w")) == NULL)
    { fprintf(stderr,"Probleme d'ouverture de %s\n",nFs);
      exit (1);  
    }
    
    //******* si on calcule Ftot
    
//      FILE *ffonction_Ftot;
//   char nFtot[50];
//   if(run.flagposanc==0)
//   {
//     sprintf(nFtot,"FsaaNa%dNb%deta%dTemp%dbstep%dstep%dtFs%d.res",run.N,run.N2,(int) (1000*run.eta),(int) (1000/run.beta),(int) run.bstep,(int) run.nstep,(int) (run.tFs));
//   }
//   else  //sinon le nouveau bstep est bstep+step_oldsimu
//   {
//     sprintf(nFtot,"FsaaNa%dNb%deta%dTemp%dbstep%dstep%dtFs%d.res",run.N,run.N2,(int) (1000*run.eta),(int) (1000/run.beta),(int) (run.bstep+run.bstep_old+run.step_old),(int) (run.nstep),(int) (run.tFs));
//   }
//             
//   if ((ffonction_Ftot = fopen (nFtot,"w")) == NULL)
//     { fprintf(stderr,"Probleme d'ouverture de %s\n",nFtot);
//       exit (1);  
//     }
  
    
    
    
  
   int t;
   int dt_Fs=(int) (run.tFs/nbre_pts);
   int dt_Fs_court=1;
   int t_min=nbre_pts*dt_Fs_court/dt_Fs;
   double Fs_0_long=(double) (compteurFs_long) *(run.N);
   //double Fs_0_court=hist_Fs_court[0];
   double Fs_0_court=(double) (compteurFs_court) *(run.N);
   
  // printf("Fs_0_long=%e\n",Fs_0_long);
  // printf("Fs_0_court=%e\n",Fs_0_court);
   
  for(t=0;t<nbre_pts;t++)
  {
   //printf("%e\n",Fs_0_court);
   fprintf(ffonction_Fs,"%e\t %e\n",  t*run.dt/run.sigma*dt_Fs_court, hist_Fs_court[t]/Fs_0_court);
  }
  
  for(t=t_min;t<nbre_pts;t++)
  {
   fprintf(ffonction_Fs, "%e\t %e\n", t*run.dt/run.sigma*dt_Fs, hist_Fs_long[t]/Fs_0_long);
  }
  fclose(ffonction_Fs);
  
//   if(run.flagFtot==1)
//   {
//    double Ftot_0_long=(double) (compteurFtot_long) *(run.N*run.N);
//    //double Fs_0_court=hist_Fs_court[0];
//    double Ftot_0_court=(double) (compteurFtot_court) *(run.N*run.N);
//    
//    for(t=0;t<nbre_pts;t++)
//    {
//     fprintf(ffonction_Ftot,"%e\t %e\n",  t*run.dt/run.sigma*dt_Fs_court, hist_Ftot_court[t]/Ftot_0_court);
//    }
//   
//   for(t=t_min;t<nbre_pts;t++)
//   {
//    fprintf(ffonction_Ftot, "%e\t %e\n", t*run.dt/run.sigma*dt_Fs, hist_Ftot_long[t]/Ftot_0_long); 
//   }
//   
//   fclose(ffonction_Ftot);
//   
//   }
   
   
   double Fs_0_long_2=(double) (compteurFs_long) *(run.N2);
   //double Fs_0_court=hist_Fs_court[0];
   double Fs_0_court_2=(double) (compteurFs_court) *(run.N2);
  
  // printf("Fs_0_long=%e\n",Fs_0_long);
  // printf("Fs_0_court=%e\n",Fs_0_court);
   if(run.N2>0)
   {
     
    FILE *ffonction_Fs_2;
    char nFs_2[50];
    if(run.flagposanc==0)
    {
     sprintf(nFs_2,"FsbbNa%dNb%deta%dTemp%dbstep%dstep%dtFs%d.res",run.N,run.N2,(int) (1000*run.eta),(int) (1000/run.beta),(int) run.bstep,(int) run.nstep,(int) (run.tFs));
    }
    else  //sinon le nouveau bstep est bstep+step_oldsimu
    {
     sprintf(nFs_2,"FsbbNa%dNb%deta%dTemp%dbstep%dstep%dtFs%d.res",run.N,run.N2,(int) (1000*run.eta),(int) (1000/run.beta),(int) (run.bstep+run.bstep_old+run.step_old),(int) (run.nstep),(int) (run.tFs));
    }
            
  if ((ffonction_Fs_2 = fopen (nFs_2,"w")) == NULL)
    { fprintf(stderr,"Probleme d'ouverture de %s\n",nFs_2);
      exit (1);  
    }   
     
     
    for(t=0;t<nbre_pts;t++)
    {
   //printf("%e\n",Fs_0_court);
     fprintf(ffonction_Fs_2,"%e\t %e\n",  t*run.dt/run.sigma*dt_Fs_court, hist_Fs_court_2[t]/Fs_0_court_2);
    }
  
  
  
    for(t=t_min;t<nbre_pts;t++)
    {
     fprintf(ffonction_Fs_2, "%e\t %e\n", t*run.dt/run.sigma*dt_Fs, hist_Fs_long_2[t]/Fs_0_long_2);
    }
   
     fclose(ffonction_Fs_2);
    } // fin if N2>0
   
   // on enregistre en plus Fs(t) dans un tableau car on s'en ressert pour khi_4(t)
   if(run.flagkhi4==1)
   {
   for(t=0;t<nbre_pts;t++)
   {
    Fs[t]=hist_Fs_court[t]/Fs_0_court;
   }
   for(t=nbre_pts;t<1000-t_min;t++)
   {
    Fs[t]=hist_Fs_long[t+t_min-nbre_pts]/Fs_0_long;
   }
   }
   
   /*! enregistrement distances */
   
   FILE *ffonction_deplacement;
  char ndeplacement[50];
  if(run.flagposanc==0)
  {

    sprintf(ndeplacement,"deplacaaNa%dNb%deta%dTemp%dbstep%dstep%dtFs%d.res",run.N,run.N2,(int) (1000*run.eta),(int) (1000/run.beta),(int) run.bstep,(int) run.nstep,(int) (run.tFs));
  }
  else  //sinon le nouveau bstep est bstep+step_oldsimu
  {
    sprintf(ndeplacement,"deplacaaNa%dNb%deta%dTemp%dbstep%dstep%dtFs%d.res",run.N,run.N2,(int) (1000*run.eta),(int) (1000/run.beta),(int) (run.bstep+run.bstep_old+run.step_old),(int) (run.nstep),(int) (run.tFs));
  }
    
    
    
  if ((ffonction_deplacement = fopen (ndeplacement,"w")) == NULL)
    { 
     fprintf(stderr,"Probleme d'ouverture de %s\n",ndeplacement);
     exit (1);  
    }

    
     
    
   int t_deplac_prime=run.nstep;
   int dt_deplac=run.tFs/nbre_pts;
   double deplac_0_long=Fs_0_long*run.sigma;
   double deplac_0_court=Fs_0_court*run.sigma;
  
     
     //  double denom_1=((double) run.nstep/dt_deplac)*run.sigma; // normalisation de <r(t)> pour petits temps
       
       
       for(t=0;t<nbre_pts;t++)
       {
	 fprintf(ffonction_deplacement,"%e\t %e\n",  t*run.dt/run.sigma*dt_Fs_court, hist_deplacement_court[t]/deplac_0_court); 
       }
       
       for(t=t_min+1;t<nbre_pts;t++)
       {
        fprintf(ffonction_deplacement, "%e\t %e\n", t*run.dt/run.sigma*dt_deplac, hist_deplacement_long[t]/deplac_0_long);
       }
     
     
    fclose(ffonction_deplacement);
   
    if(run.N2>0)
    {
     int dt_Fs=(int) (run.tFs/nbre_pts);
     int dt_Fs_court=1;
     int t_min=nbre_pts*dt_Fs_court/dt_Fs;
     double Fs_0_long_2=(double) (compteurFs_long) *(run.N2);
     double Fs_0_court_2=(double) (compteurFs_court) *(run.N2);
     double deplac_0_long_2=Fs_0_long_2*run.sigma2;
     double deplac_0_court_2=Fs_0_court_2*run.sigma2;
  
     
     //  double denom_1=((double) run.nstep/dt_deplac)*run.sigma; // normalisation de <r(t)> pour petits temps
             
     
     
     
    
    
    FILE *ffonction_deplacement_2;
     char ndeplacement_2[50];
     if(run.flagposanc==0)
     {
      sprintf(ndeplacement_2,"deplacbbNa%dNb%deta%dTemp%dbstep%dstep%dtFs%d.res",run.N,run.N2,(int) (1000*run.eta),(int) (1000/run.beta),(int) run.bstep,(int) run.nstep,(int) (run.tFs));
     }
     else  //sinon le nouveau bstep est bstep+step_oldsimu
     {
      sprintf(ndeplacement_2,"deplacbbNa%dNb%deta%dTemp%dbstep%dstep%dtFs%d.res",run.N,run.N2,(int) (1000*run.eta),(int) (1000/run.beta),(int) (run.bstep+run.bstep_old+run.step_old),(int) (run.nstep),(int) (run.tFs));
     }
    
    
    
     if ((ffonction_deplacement_2 = fopen (ndeplacement_2,"w")) == NULL)
     { 
      fprintf(stderr,"Probleme d'ouverture de %s\n",ndeplacement_2);
      exit (1);  
     } 
    
     for(t=0;t<nbre_pts;t++)
     {
      fprintf(ffonction_deplacement_2,"%e\t %e\n",  t*run.dt/run.sigma*dt_Fs_court, hist_deplacement_court_2[t]/deplac_0_court_2); 
     }
       
     for(t=t_min+1;t<nbre_pts;t++)
     {
      fprintf(ffonction_deplacement_2, "%e\t %e\n", t*run.dt/run.sigma*dt_deplac, hist_deplacement_long_2[t]/deplac_0_long_2);
     }
    
    fclose(ffonction_deplacement_2);
    
    } // fin if N2>0
   
 } // fin if(flag==2)
 
/*
free(hist_pos_x);
free(hist_pos_y);
free(hist_pos_z);
free(hist_pos_x_court);
free(hist_pos_y_court);
free(hist_pos_z_court);
free(hist_Fs);
free(hist_Fs_prime);*/


} // fin calcul_Fs

//**************************

/*! F_tot(k,t) */
void calcul_Ftot(int flag, int step, int dt_Ftot, POS_cart *pos2_cart, double *arccos, double *legendre, RUN run)
{
 int i;
 int nbre_pts=100;
 //int dt_Fs=(int) (run.tFs/500);
 int itt1=step/dt_Ftot;
 int itt0=itt1% nbre_pts;
 int itt2;
 int itt3;
 static double **hist_pos_x=NULL; //histo des coordonnees (pour calcul des grandeurs dynamiques)
 static double **hist_pos_y=NULL;
 static double **hist_pos_z=NULL;
 
 static double **hist_pos_x_court=NULL; //histo des coordonnÃ©es phi (pour calcul des grandeurs dynamiques)
 static double **hist_pos_y_court=NULL;
 static double **hist_pos_z_court=NULL; 
 static double *hist_Ftot_court=NULL;
 static double *hist_Ftot_long=NULL;


 
 double startwtime,endwtime;
 static int compteurFtot_court;
 static int compteurFtot_long;
 POS_cart pos1;
 POS_cart pos2;
 
 if(flag==0)
 {
 allocdouble2(&hist_pos_x,run.N+run.N2,nbre_pts);
 allocdouble2(&hist_pos_y,run.N+run.N2,nbre_pts);
 allocdouble2(&hist_pos_z,run.N+run.N2,nbre_pts);  
 allocdouble2(&hist_pos_x_court,run.N+run.N2,nbre_pts);
 allocdouble2(&hist_pos_y_court,run.N+run.N2,nbre_pts);
 allocdouble2(&hist_pos_z_court,run.N+run.N2,nbre_pts);
 allocdouble(&hist_Ftot_long,nbre_pts);
 allocdouble(&hist_Ftot_court,nbre_pts);

 
 
 int i;
 int j; 
 compteurFtot_court=0;
 compteurFtot_long=0;

 } // fin if flag=0
 
 if(flag==1)
 {
   itt1=step/dt_Ftot;
   itt0=itt1% nbre_pts;

   if(dt_Ftot==(int) (run.tFs/nbre_pts))
   {
#pragma omp parallel for default(shared)   private(i)   
  for(i=0;i<run.N+run.N2;i++) //il faut enregistrer toutes les positions pour le calcul de Fs(k,t) !!!
 {
  hist_pos_x[i][itt0]=pos2_cart[i].x;
  hist_pos_y[i][itt0]=pos2_cart[i].y;
  hist_pos_z[i][itt0]=pos2_cart[i].z;
 }
 
 itt1++;  itt0=itt1% nbre_pts;
   if (itt1>=nbre_pts)
   {
	 int intdeb=  nbre_pts/ dt_Ftot;
	 compteurFtot_long++;
  for (itt2=intdeb;itt2<nbre_pts;itt2++)
  {   
	      itt3=(itt0+itt2) % nbre_pts;
	      double hist_cur_tot=0;

	 

  int j;
//#pragma omp parallel for default(shared) private (i,j,pos1,pos2) reduction(+:hist_cur_tot)
 for(i=0;i<run.N;i++)
 {
   for(j=0;j<run.N;j++)
   {
   pos1.x=hist_pos_x[i][itt0];
   pos1.y=hist_pos_y[i][itt0];
   pos1.z=hist_pos_z[i][itt0];	
   pos2.x=hist_pos_x[j][itt3];
   pos2.y=hist_pos_y[j][itt3];
   pos2.z=hist_pos_z[j][itt3];
   double dist=distance_cart(pos1,pos2,arccos,run);
   hist_cur_tot+=calcul_legendre(dist,legendre,run);
  }
 }   
//#pragma omp barrier
 hist_Ftot_long[itt2]+=hist_cur_tot;

	  
  } // fin for itt2=indeb;..
	 
	 
	  
    } // fin if(itt1>=500)
       

       
    //   if (omp_get_thread_num() == 0)    endwtime =  omp_get_wtime(); 
//if (omp_get_thread_num() == 0) printf("calcul Fs time  = %e\n",(double)(endwtime-startwtime));
   } // fin si if dt_Fs=tFs/500
   
// if (omp_get_thread_num() == 0)    endwtime =  omp_get_wtime(); 
// if (omp_get_thread_num() == 0) printf("calcul Fs  2 time  = %e\n",(double)(endwtime-startwtime)); 
//  
 } //fin du if (flag == 1)
 
//************************

 if(flag==2)
 {

    //******* si on calcule Ftot
    
  FILE *ffonction_Ftot;
  char nFtot[50];
  if(run.flagposanc==0)
  {
    sprintf(nFtot,"FtotNa%dNb%deta%dTemp%dbstep%dstep%dtFs%d.res",run.N,run.N2,(int) (1000*run.eta),(int) (1000/run.beta),(int) run.bstep,(int) run.nstep,(int) (run.tFs));
  }
  else  //sinon le nouveau bstep est bstep+step_oldsimu
  {
    sprintf(nFtot,"FtotNa%dNb%deta%dTemp%dbstep%dstep%dtFs%d.res",run.N,run.N2,(int) (1000*run.eta),(int) (1000/run.beta),(int) (run.bstep+run.bstep_old+run.step_old),(int) (run.nstep),(int) (run.tFs));
  }
            
  if ((ffonction_Ftot = fopen (nFtot,"w")) == NULL)
    { fprintf(stderr,"Probleme d'ouverture de %s\n",nFtot);
      exit (1);  
    }
  
    
    
    
  
   int t;
   int dt_Ftot=(int) (run.tFs/nbre_pts);
//    int dt_Fs_court=1;
//    int t_min=500*dt_Fs_court/dt_Fs;
//    double Fs_0_long=(double) (compteurFs_long) *(run.N);
//    //double Fs_0_court=hist_Fs_court[0];
//    double Fs_0_court=(double) (compteurFs_court) *(run.N);


   double Ftot_0_long=(double) (compteurFtot_long) *(run.N);
  
   
  for(t=0;t<nbre_pts;t++)
  {
   fprintf(ffonction_Ftot,"%e\t %e\n", t*run.dt/run.sigma*dt_Ftot, hist_Ftot_long[t]/Ftot_0_long);
  }
  
  
  fclose(ffonction_Ftot);
   
 } // fin if(flag==2)
 
 
 
} // fin calcul_Ftot



void calcul_khi4(int flag, int step, POS_cart *pos2_cart, double *Fs, double *arccos, double *legendre, RUN run)
{
 int i;
 int j;
 int nbre_pts=50;
 int dt_khi4=(int) (run.tFs/nbre_pts);
 int itt1=step/dt_khi4;
 int itt0=itt1% nbre_pts;
 int itt2;
 int itt3;
 int t;
 
 static double **hist_pos_x=NULL; //histo des coordonnées phi (pour calcul des grandeurs dynamiques)
 static double **hist_pos_y=NULL;
 static double **hist_pos_z=NULL;

 static double *hist_Fs=NULL; // on a besoin de Fs pour Khi4 : khi4(t)=<fs(t)^2>-Fs(t)^2
 static double *hist_khi4=NULL; 
 
 static int compteur_Fs;
 static int compteur_khi4;

 POS_cart pos1;
 POS_cart pos2;
 POS_cart pos1_prime; // seconde position à t' pour calcul de Khi4
 POS_cart pos2_prime; // seconde poition à t'+t pour calcul de Khi4
 
 if(flag==0)
 {
 allocdouble2(&hist_pos_x,run.N+run.N2,nbre_pts);
 allocdouble2(&hist_pos_y,run.N+run.N2,nbre_pts);
 allocdouble2(&hist_pos_z,run.N+run.N2,nbre_pts);  

 allocdouble(&hist_Fs,nbre_pts);
 allocdouble(&hist_khi4,nbre_pts);
 
 compteur_Fs=0;
 compteur_khi4=0;
 }

 if(flag==1)
 {
   itt1=step/dt_khi4;
   itt0=itt1%nbre_pts;
 #pragma omp parallel for default(shared) private (i) 
 for(i=0;i<run.N+run.N2;i++)
 {
  hist_pos_x[i][itt0]=pos2_cart[i].x;
  hist_pos_y[i][itt0]=pos2_cart[i].y;
  hist_pos_z[i][itt0]=pos2_cart[i].z; 
 }
  
  itt1++; itt0=itt1%nbre_pts;
  if(itt1>=nbre_pts)
  {
    int intdeb=nbre_pts/dt_khi4;
    compteur_khi4++; 
    compteur_Fs++;
    
    for(itt2=intdeb;itt2<nbre_pts;itt2++)
    {
     itt3=(itt0+itt2) % nbre_pts;
     double hist_cur=0; // Fs
     double hist_cur2=0; // khi4
    
   
 for(i=0;i<run.N+run.N2;i++) //il faut enregistrer toutes les positions pour le calcul de Fs(k,t) !!!
 {	   
//        if (itt1>=nbre_pts)
//        {// iitnumber++;
       
        pos1.x=hist_pos_x[i][itt0];
	pos1.y=hist_pos_y[i][itt0];
        pos1.z=hist_pos_z[i][itt0];

	pos2.x=hist_pos_x[i][itt3];
        pos2.y=hist_pos_y[i][itt3];
        pos2.z=hist_pos_z[i][itt3];
	   
	double dist=distance_cart(pos1,pos2,arccos,run);
	double Pk_dist=calcul_legendre(dist,legendre,run);
        hist_cur+=Pk_dist;
	hist_cur2+=Pk_dist*Pk_dist; 
	//   hist_khi4[itt2]+=Pk_dist*Pk_dist; // terme pour i=j

	   
         for(j=i+1;j<run.N+run.N2;j++) // i et j jouent des rôles symétriques (ne pas oublier le terme pour i=j)
	 {
           
	   
	   pos1_prime.x=hist_pos_x[j][itt0];
	   pos1_prime.y=hist_pos_y[j][itt0];
           pos1_prime.z=hist_pos_z[j][itt0];
 
	   pos2_prime.x=hist_pos_x[j][itt3];
           pos2_prime.y=hist_pos_y[j][itt3];
           pos2_prime.z=hist_pos_z[j][itt3];
	      
	   double dist_prime=distance_cart(pos1_prime,pos2_prime,arccos,run);
	   double Pk_dist_bis=calcul_legendre(dist_prime,legendre,run);
 
	   hist_cur2+=2*Pk_dist*Pk_dist_bis;
	     // hist_khi4[itt2]+=2*Pk_dist*Pk_dist_bis; // terme pour i<j
	      
		 // on calcule en fait <fs(t)> et <fs(t)^2> pour l'instant
	   
	      
	  } // fin boucle sur j
	 } // fin boucle sur i	
	 hist_Fs[itt2]+=hist_cur;
	 hist_khi4[itt2]+=hist_cur2;
	} // fin boucle sur itt2
 } // fin if itt1>=nbre_pts
 

 } // fin if flag==1
 
 if(flag==2)
 {
   FILE *ffonction_Khi4;
  char nKhi4[50];
  if(run.flagposanc==0)
  {

    sprintf(nKhi4,"khi4Na%dNb%deta%dTemp%dbstep%dstep%dtFs%d.res",run.N,run.N2,(int) (1000*run.eta),(int) (1000/run.beta),(int) run.bstep,(int) run.nstep,(int) (run.tFs));
  }
  else  //sinon le nouveau bstep est bstep+step_oldsimu
  {
    sprintf(nKhi4,"khi4Na%dNb%deta%dTemp%dbstep%dstep%dtFs%d.res",run.N,run.N2,(int) (1000*run.eta),(int) (1000/run.beta),(int) (run.bstep+run.bstep_old+run.step_old),(int) (run.nstep),(int) (run.tFs));
  }
    
    
    
  if ((ffonction_Khi4 = fopen (nKhi4,"w")) == NULL)
    { fprintf(stderr,"Probleme d'ouverture de %s\n",nKhi4);
      exit (1);  
    }
    
  int dt_khi4=(int) (run.tFs/nbre_pts);
  double Khi4_0=(double) (compteur_khi4)*(run.N+run.N2)*(run.N+run.N2);
  double Fs_0=(double) (compteur_Fs)*(run.N+run.N2);
     
  for(t=0;t<nbre_pts;t++)
  {
   fprintf(ffonction_Khi4,"%e\t %e\n",  t*run.dt/run.sigma*dt_khi4, (run.N+run.N2)*(hist_khi4[t]/Khi4_0-(hist_Fs[t]/Fs_0)*(hist_Fs[t]/Fs_0))); // khi4(t)=<fs^2(t)>-<fs(t)>^2
  }
  fclose(ffonction_Khi4); 
 } // fin if flag=2
        
 
} // fin calcul_khi4

    


  




