double min_acos; // valeur minimale dans tabulation de arccos : min_acos=cos(2.5*sigma)
double max_acos; // valeur maximale dans tabulation de arccos : max_acos=cos(0.5*sigma)

void read_data(int argc, char *argv[], RUN *run)
{
        FILE            *debut;
        if (argc != 2)
                        {
                        printf("usage: mcsphere  <fichier de demarrage>\n");
                        exit (1);
                        }

    /* ouverture du fichier d'entree */
        if ((debut = fopen (argv[1],"r")) == NULL)
                 { fprintf(stderr,"Probleme d'ouverture de %s\n",argv[1]);
                   exit (1);
                 }
         fscanf(debut,"%d%*[^\n]", &(run->N)); printf("Nbre particles %d\n", run->N);
	 fscanf(debut,"%d%*[^\n]", &(run->N2)); printf("Nbre particles 2 %d\n", run->N2);
	 fscanf(debut,"%d%*[^\n]", &(run->nthreads)); printf("Nbre de threads %d\n", run->nthreads);
         fscanf(debut,"%ld%*[^\n]", &(run->bstep)); printf("Conf avant %ld\n", run->bstep);
	 fscanf(debut,"%ld%*[^\n]", &(run->nstep)); printf("Conf apres %ld\n", run->nstep);
	 fscanf(debut,"%lf%*[^\n]", &(run->eta)); printf("eta %lf\n", run->eta);
	 //fscanf(debut,"%lf%*[^\n]", &(run->eta2));
	 fscanf(debut,"%lf%*[^\n]", &(run->temp)); printf("Temperature %lf\n", run->temp);
	 run->beta=1/run->temp;
	 fscanf(debut,"%lf%*[^\n]", &(run->dt)); printf("dt %lf\n", run->dt);
	 //fscanf(debut,"%lf%*[^\n]", &(run->beta)); printf("beta %lf\n", run->beta);
         //run->sigma=2*acos(1-2*run->eta/(double) run->N);  printf("diameter %lf\n", run->sigma);
	 run->sigma=2*acos(1-pi*run->eta/(2*(double) run->N));  printf("diameter %lf\n", run->sigma);
	 //run->sigma2=2*acos(1-pi*run->eta2/(2*(double) run->N2));  printf("diameter2 %lf\n", run->sigma2);
	 run->sigma2=0.7*run->sigma; printf("diameter2 %lf\n", run->sigma2);
	 run->eta2=2*run->N2/pi*(1-cos(run->sigma2/2));
	 run->sigma12=0.85*run->sigma; printf("diameter12 %lf\n", run->sigma12);
	 //run->sigma=1;
	 printf("eta %lf\n", run->eta);
	 printf("eta2 %lf\n", run->eta2);
	 
	 
	 fscanf(debut,"%lf%*[^\n]", &(run->deltaphi)); //printf("delta phi %lf\n", run->deltaphi);
         fscanf(debut,"%lf%*[^\n]", &(run->deltaz)); //printf("delta z %lf\n", run->deltaz);
         fscanf(debut,"%lf%*[^\n]", &(run->nrjbas)); //printf("lower bound energy %f\n", run->nrjbas);
         fscanf(debut,"%lf%*[^\n]", &(run->nrjhaut)); //printf("upper bound energy %f\n", run->nrjhaut);
	 fscanf(debut,"%lf%*[^\n]", &(run->deltar)); printf("deltar %f\n", run->deltar);
	 fscanf(debut,"%lf%*[^\n]", &(run->kFs)); printf("kFs %lf\n", run->kFs);
	 fscanf(debut,"%d%*[^\n]", &(run->kmax)); printf("kmax %d\n", run->kmax);
	 fscanf(debut,"%d%*[^\n]", &(run->kcut)); printf("kcut %d\n", run->kcut);
	 fscanf(debut,"%d%*[^\n]", &(run->tFs)); printf("tFs %d\n", run->tFs);
	 fscanf(debut,"%d%*[^\n]", &(run->dtmoy)); printf("dtmoy %d\n", run->dtmoy);
	 fscanf(debut,"%d%*[^\n]", &(run->nconf)); printf("nconf %d\n", run->nconf);
	 fscanf(debut,"%d%*[^\n]", &(run->moviestep)); printf("movie step %d\n", run->moviestep);
	 fscanf(debut,"%d%*[^\n]", &(run->flagmovie)); printf("flagmovie %d\n", run->flagmovie);

	 
	 
	 fscanf(debut,"%d%*[^\n]", &(run->flagposanc)); printf("flagposanc %d\n", run->flagposanc);
	 fscanf(debut,"%d%*[^\n]", &(run->flagtherma)); printf("flagtherma %d\n", run->flagtherma);
	 fscanf(debut,"%d%*[^\n]", &(run->flagcorr)); printf("flagcorr %d\n", run->flagcorr);
	// fscanf(debut,"%d%*[^\n]", &(run->flagcorrdefauts)); printf("flagcorrdefauts %d\n", run->flagcorrdefauts);
	 fscanf(debut,"%d%*[^\n]", &(run->flagstruc)); printf("flagstruc %d\n", run->flagstruc);
	 fscanf(debut,"%d%*[^\n]", &(run->flagg6)); printf("flagg6 %d\n", run->flagg6);
	 fscanf(debut,"%d%*[^\n]", &(run->flagFs)); printf("flagFs %d\n", run->flagFs);
	 fscanf(debut,"%d%*[^\n]", &(run->flagFtot)); //printf("flagFtot %d\n", run->flagFtot);
	 fscanf(debut,"%d%*[^\n]", &(run->flagkhi4)); //printf("flagkhi4 %d\n", run->flagkhi4);
	 fscanf(debut,"%d%*[^\n]", &(run->flagdeplac)); printf("flagdeplac %d\n", run->flagdeplac);
	 fscanf(debut,"%d%*[^\n]", &(run->flagdeplac657)); printf("flagdeplac657 %d\n", run->flagdeplac657);
	 fscanf(debut,"%d%*[^\n]", &(run->flagdevoro)); printf("flagdevoro %d\n", run->flagdevoro);
	 fscanf(debut,"%d%*[^\n]", &(run->flagsavedefauts)); printf("flagsavedefauts %d\n", run->flagsavedefauts);
	 fscanf(debut,"%d%*[^\n]", &(run->flagspectre)); printf("flagspectre %d\n", run->flagspectre);
	 fscanf(debut,"%d%*[^\n]", &(run->flagdynam)); printf("flagdynam %d\n", run->flagdynam);
	 fscanf(debut,"%d%*[^\n]", &(run->flagkwave)); printf("flagkwave %d\n", run->flagkwave);
	 fscanf(debut,"%d%*[^\n]", &(run->bstep_old));
	 fscanf(debut,"%d%*[^\n]", &(run->step_old));
	 run->maxpot=2000000; printf("maxpot %d\n", run->maxpot);
	 run->maxhisto=run->N*10;
	 run->maxacos=4000000; printf("maxacos %d\n", run->maxacos);
	 run->maxlegendre=1000000;
	 //prod_scal_max=cos(2.5*run->sigma);
	 min_acos=cos(2.5*run->sigma);
	 //min_acos=0;
	 max_acos=cos(0.5*run->sigma);
}  


//*******************************************


/*! Positions initiales */
void  initial_config(POS *pos,RUN run,gsl_rng *rng) //1ere et 3eme variables definies dans la fonction, mais type defini en dehors
{
  int i;
  for (i=0;i<run.N+run.N2;i++)
  {
    pos[i].phi=pi2*gsl_rng_uniform(rng);
    pos[i].z=2*gsl_rng_uniform(rng)-1.0;
  }
#if DEBUG
{
FILE *flog;
char nlog[]="configini.res";
if ((flog = fopen (nlog,"w")) == NULL)
    { 
     fprintf(stderr,"Probleme d'ouverture de %s\n",nlog);2
     exit (1);  
    }
   int i;
  for (i=0;i<run.N+run.N2;i++)
  { 
   fprintf(flog,"%e %e %e\n",sqrt(1-pos[i].z*pos[i].z)*cos(pos[i].phi),
		sqrt(1-pos[i].z*pos[i].z)*sin(pos[i].phi),pos[i].z);  
  }
  fclose(flog);
}
#endif
} // fin initial_config



//*****************************************************************

/*!Reprise d'une ancienne config finale comme nouvelle config initiale */

void initial_config_anc(POS_cart *pos, THERMO *thermo, RUN run)
{
  FILE *fposanc;
  char nposanc[50];
  int i; int ihist;
  sprintf(nposanc,"posfinNa%dNb%deta%dTemp%dbstep%dstep%d.res",run.N,run.N2,(int) (1000*run.eta), (int) (1000/run.beta),run.bstep_old,run.step_old);
  if((fposanc=fopen(nposanc,"r"))==NULL)
  {
    printf("Problème d'ouverture de posfin");
    exit(1);
  }
 
   for(i=0;i<run.N+run.N2;i++)
  {
    //fscanf(fposanc,"%le\t %le\n",&pos_init_phi[i],&pos_init_z[i]);
    fscanf(fposanc,"%le\t %le\t %le\n",&pos[i].x,&pos[i].y,&pos[i].z);
    //printf("pos.phi=%e pos.z=%e\n",pos[i].phi,pos[i].z);2
  }
  fclose(fposanc);
 
  
  if(run.flagtherma==1) //on lit l'histo des energies, nrjac et nrjac2
  {
  FILE *fhistoanc;
  char nhistoanc[50];
  sprintf(nhistoanc,"histoancNa%dNb%deta%dTemp%dbstep%dstep%d.res",run.N,run.N2,(int) (1000*run.eta),(int) (1000/run.beta),run.bstep_old,run.step_old);
  if((fhistoanc=fopen(nhistoanc,"r"))==NULL)
  {
    printf("Probleme d'ouverture de histofin");
    exit(1);
  }
  

  for(ihist=0;ihist<run.maxhisto;ihist++)
  {
    fscanf(fhistoanc,"%le\n",&thermo->histoener[ihist]);
  }
  fclose(fhistoanc);
  
   FILE *fenergacc;
	 char nenergacc[50];
	 sprintf(nenergacc,"energaccNa%dNb%deta%dTemp%dbstep%dstep%d.res",run.N,run.N2,(int) (1000*run.eta),(int) (1000/run.beta),run.bstep_old,run.step_old);
	 if((fenergacc=fopen(nenergacc,"r"))==NULL)
         {
          printf("Probleme d'ouverture de energacc");
          exit(1);
         }
         fscanf(fenergacc,"%le\t%le\n",&thermo->nrjac,&thermo->nrjac2);
  fclose(fenergacc);
  }
} // fin initial_config_anc


//******************************************************************
//******************************************************************

/*! Initialisation des pos et des vitesses */

void dyn_molec_init(POS_cart *pos1_cart, POS_cart *omg1_cart, POS_ang *pos1, POS_ang *vit1, double *arccos, RUN run, gsl_rng *rng)
{
  double sigma_min;
  if(run.sigma<run.sigma2){sigma_min=run.sigma;}else{sigma_min=run.sigma2;}
  double dt=run.dt;
  //initialement, pos1 -> positions Ã  t=-dt, pos2 -> positions Ã  t=0
  int i; int j;
  double sumvphi=0; // somme des composantes des vitesses des part suivant phi, Ã  retirer de chaque vphi,i
  double sumvtheta=0;
  double d_min_1=10; // initialise a une valeur grande
  double d_min_2=10; // initialise a une valeur grande
  double d_min_12=10; // initialise a une valeur grande

  for(i=0;i<run.N+run.N2;i++)
  {
    double u=gsl_rng_uniform(rng);
    double v=gsl_rng_uniform(rng); // u et v sont des variables uniformes sur [0,1]
    vit1[i].phi=sqrt(-2*run.temp*log(u))*cos(pi2*v);
    vit1[i].theta=sqrt(-2*run.temp*log(u))*sin(pi2*v); // vphi et vz sont des variables gaussiennes indep N(0,T)
    sumvphi+=vit1[i].phi;
    sumvtheta+=vit1[i].theta;   
  }
  
 if (run.flagdynam==0) 
 {
   int ii=0;
   i=0;
   int k=0;
   bool overlap=false;

// on remplace le dÃÂ©pot alÃÂ©atoire par une configuration de type RSA
  
  /*! remplissage par les particules 1 */
  while(ii<run.N+run.N2)
  {  
   k++;
   overlap=false;
   pos1[ii].phi=pi2*gsl_rng_uniform(rng);
   pos1[ii].theta=pi*gsl_rng_uniform(rng);
   i=0;
   if (ii != 0)
   {
    do
    {
//      if (distance_ang(pos1[i],pos1[ii],arccos,run) < run.sigma* 0.85){overlap=true;}
//      ++i;
     if ((i<run.N)&&(ii<run.N)&&(distance_ang(pos1[i],pos1[ii],arccos,run) < run.sigma* 0.85)){overlap=true;}
     if ((i<run.N)&&(ii>=run.N)&&(distance_ang(pos1[i],pos1[ii],arccos,run) < run.sigma12* 0.9)){overlap=true;}
     if ((i>=run.N)&&(ii<run.N)&&(distance_ang(pos1[i],pos1[ii],arccos,run) < run.sigma12* 0.9)){overlap=true;}
     if ((i>=run.N)&&(ii>=run.N)&&(distance_ang(pos1[i],pos1[ii],arccos,run) < run.sigma2* 0.92)){overlap=true;}
     ++i;
    } while(i<ii && overlap == false);
   }
        if (overlap == false) {++ii;/*printf("nombre particules %d\n",ii)*/;}
   } // fin while
    
    if(k==1000)
    {
     printf("Impossible de placer toutes les particules 1\t, k=%d\n essais", k);
     exit(1);
    }
 printf("remplissage des particules fini\n");
  
 

   for(i=0;i<run.N-1;i++) // on recalcule d_min
   {  
    for(j=i+1;j<run.N;j++)
    {  
     if(distance_ang(pos1[i],pos1[j],arccos,run)<d_min_1){d_min_1=distance_ang(pos1[i],pos1[j],arccos,run);}    
    }
   }
   
   for(i=0;i<run.N;i++) // on recalcule d_min
   {  
    for(j=run.N;j<run.N+run.N2;j++)
    {  
     if(distance_ang(pos1[i],pos1[j],arccos,run)<d_min_12){d_min_12=distance_ang(pos1[i],pos1[j],arccos,run);}    
    }
   }
   
   for(i=run.N;i<run.N+run.N2-1;i++) // on recalcule d_min
   {  
    for(j=i+1;j<run.N+run.N2;j++)
    {  
     if(distance_ang(pos1[i],pos1[j],arccos,run)<d_min_2){d_min_2=distance_ang(pos1[i],pos1[j],arccos,run);}    
    }
   }
   
   
  
 } // fin if(run.flagdynam==0)
 else
 {
  FILE *fconf;
  char nconf[50];
  sprintf(nconf,"confiniNa%dNb%d.res",run.N,run.N2);
if ((fconf = fopen (nconf,"r")) == NULL)
    { 
     fprintf(stderr,"Probleme d'ouverture de %s\n",nconf);
     exit(1);  
    }
    int i;
    for (i=0;i<run.N+run.N2;i++)
    { 
     fscanf(fconf,"%le %le %le\n",&(pos1_cart[i].x),&(pos1_cart[i].y),&(pos1_cart[i].z)); 
    }
  fclose(fconf);
 }

  /*!pour eviter que les particules n'aient de mouvement global, on retire a chaque composante de la vitesse de i la somme
    des comp sur ttes les particules*/
    
    for(i=0;i<run.N+run.N2;i++)
    {
     vit1[i].phi-=1/((double) run.N+run.N2)*sumvphi;
     vit1[i].theta-=1/((double) run.N+run.N2)*sumvtheta;
     //printf("vphi=%e  vtheta=%e\n",vit1[i].phi,vit1[i].theta);
    }
 if(run.flagdynam==0) 
 {
    for(i=0;i<run.N+run.N2;i++)
    {
     pos1_cart[i].x=sin(pos1[i].theta)*cos(pos1[i].phi);
     pos1_cart[i].y=sin(pos1[i].theta)*sin(pos1[i].phi);
     pos1_cart[i].z=cos(pos1[i].theta);  
     
     omg1_cart[i].x=cos(pos1[i].theta)*cos(pos1[i].phi)*vit1[i].theta-sin(pos1[i].phi)*vit1[i].phi;
     omg1_cart[i].y=cos(pos1[i].theta)*sin(pos1[i].phi)*vit1[i].theta+cos(pos1[i].phi)*vit1[i].phi;
     omg1_cart[i].z=-sin(pos1[i].theta)*vit1[i].theta;   
    // if(i>=run.N) printf("i=%d\t ||r,i||=%e\n", i, prod_scalaire(pos1_cart[i],pos1_cart[i]));
    }    
    printf("Distance initiale minimale/sigma entre 2 particules 1 %e\t\n",d_min_1/run.sigma);
    printf("Distance initiale minimale/sigma2 entre 2 particules 2 %e\t\n",d_min_2/run.sigma2);
    printf("Distance initiale minimale/run.sigma12 entre 1 particule 1 et 1 particule 2 %e\t\n",d_min_12/run.sigma12);
    FILE *fconf;
    char nconf[50];
    sprintf(nconf,"confiniNa%dNb%d.res",run.N,run.N2);
if ((fconf = fopen (nconf,"w")) == NULL)
    {
     fprintf(stderr,"Probleme d'ouverture de %s\n",nconf);
     exit(1);  
    }

  for(i=0;i<run.N+run.N2;i++)
  { 
   fprintf(fconf,"%e %e %e\n",pos1_cart[i].x,pos1_cart[i].y,pos1_cart[i].z);  
  }
  fclose(fconf);
 }
 else
 {
  for(i=0;i<run.N+run.N2;i++)
    {
     pos1[i].theta=acos(pos1_cart[i].z);
     if (pos1_cart[i].x!=0.0) {
     pos1[i].phi=atan(pos1_cart[i].y/pos1_cart[i].x);
     }
     else
     {
      if (pos1_cart[i].y>0) {pos1[i].phi=0.5*M_PI;} else  {pos1[i].phi=-0.5*M_PI;} 
     }
     omg1_cart[i].x=cos(pos1[i].theta)*cos(pos1[i].phi)*vit1[i].theta-sin(pos1[i].phi)*vit1[i].phi;
     omg1_cart[i].y=cos(pos1[i].theta)*sin(pos1[i].phi)*vit1[i].theta+cos(pos1[i].phi)*vit1[i].phi;
     omg1_cart[i].z=-sin(pos1[i].theta)*vit1[i].theta;      
    } 
 }
 
    
} // fin dyn_molec_init


//*******************************

/*! fonction qui reinitialise les vitesses suivant une loi gaussienne dÃ©pendant de T */

void reinit_vit(POS_cart *pos1_cart, POS_cart *omg1_cart,  POS_ang *pos1, POS_ang *vit1, double *arccos, RUN run, gsl_rng *rng)
{
  int i;
  double sumvphi=0; // somme des composantes des vitesses des part suivant phi, a retirer de chaque vphi,i
  double sumvtheta=0;
  
  
  for(i=0;i<run.N+run.N2;i++)
  {
    
    double u=gsl_rng_uniform(rng);
    double v=gsl_rng_uniform(rng); // u et v sont des variables uniformes sur [0,1]
    vit1[i].phi=sqrt(-1*run.temp*log(u))*cos(pi2*v);
    vit1[i].theta=sqrt(-1*run.temp*log(u))*sin(pi2*v); // vphi et vtheta sont des variables gaussiennes indep N(0,T)
    //vit[i].z=cos(vtheta); //vz=sin(theta)*vtheta
    sumvphi+=vit1[i].phi;
    sumvtheta+=vit1[i].theta;   
  }
  
  for(i=0;i<run.N+run.N2;i++)
    {
     vit1[i].phi-=1/((double) run.N+run.N2)*sumvphi;
     vit1[i].theta-=1/((double) run.N+run.N2)*sumvtheta;
     //printf("vphi=%e  vtheta=%e\n",vit1[i].phi,vit1[i].theta);
    }
    for(i=0;i<run.N+run.N2;i++)
    {
     
     //if(i>=run.N) printf("i%d\t pos1[i]%e\n",i,pos1_cart[i].z);
      pos1[i].theta=arccosinus(pos1_cart[i].z,arccos,run);
      if(pos1_cart[i].y>=0)
      {
       pos1[i].phi=arccosinus(pos1_cart[i].x/sqrt(pos1_cart[i].x*pos1_cart[i].x+pos1_cart[i].y*pos1_cart[i].y),arccos,run);
      }
      else
      {
       pos1[i].phi=pi2-arccosinus(pos1_cart[i].x/sqrt(pos1_cart[i].x*pos1_cart[i].x+pos1_cart[i].y*pos1_cart[i].y),arccos,run);
      }
      
   omg1_cart[i].x=cos(pos1[i].theta)*cos(pos1[i].phi)*vit1[i].theta-sin(pos1[i].phi)*vit1[i].phi;
   omg1_cart[i].y=cos(pos1[i].theta)*sin(pos1[i].phi)*vit1[i].theta+cos(pos1[i].phi)*vit1[i].phi;
   omg1_cart[i].z=-sin(pos1[i].theta)*vit1[i].theta;

    }
    
    
    POS_cart L_cin; L_cin.x=0; L_cin.y=0; L_cin.z=0;
    
    for(i=0;i<run.N+run.N2;i++)
    {
    // L_cin=somme_vect(L_cin,prod_vect(pos1_cart[i],omg1_cart[i]));
     L_cin=somme_vect(L_cin,omg1_cart[i]);
    }
    //printf("Lx=%e Ly=%e Lz=%e\n\n",L_cin.x,L_cin.y,L_cin.z);
    
    
  /*! annulation du moment cinÃ©tique */
    double a=0; double b=0; double c=0; double d=0; double e=0; double f=0; //coefficients de la matrice Ã  inverser
    
    for(i=0;i<run.N+run.N2;i++)
    {
     a+=pos1_cart[i].x*pos1_cart[i].x-1;
     b+=pos1_cart[i].x*pos1_cart[i].y;
     c+=pos1_cart[i].x*pos1_cart[i].z;
     d+=pos1_cart[i].y*pos1_cart[i].y-1;
     e+=pos1_cart[i].y*pos1_cart[i].z;
     f+=pos1_cart[i].z*pos1_cart[i].z-1;
    }
      
      //printf("a=%e\n",a); printf("b=%e\n",b); printf("c=%e\n",c); printf("d=%e\n",d); printf("e=%e\n",e);
      //printf("f=%e\n",f);
      
      
    POS_cart *omg=NULL; // vecteur intervenant dans le changement de coord qui sert Ã  annuler le moment cinÃ©tique
    omg=(POS_cart *) malloc(sizeof(POS_cart)*1); 
    //  omg[0].x=0; omg[0].y=0; omg[0].z=0;
      
    inversion_systeme(omg, a,b,c,d,e,f,L_cin.x,L_cin.y,L_cin.z);
  
    for(i=0;i<run.N+run.N2;i++)
    {
     omg1_cart[i].x+=omg[0].x-(pos1_cart[i].x*omg[0].x+pos1_cart[i].y*omg[0].y+pos1_cart[i].z*omg[0].z)*pos1_cart[i].x; 
     omg1_cart[i].y+=omg[0].y-(pos1_cart[i].x*omg[0].x+pos1_cart[i].y*omg[0].y+pos1_cart[i].z*omg[0].z)*pos1_cart[i].y; 
     omg1_cart[i].z+=omg[0].z-(pos1_cart[i].x*omg[0].x+pos1_cart[i].y*omg[0].y+pos1_cart[i].z*omg[0].z)*pos1_cart[i].z; 
    }
    
    //********************
  
  POS_cart L_cin_bis; L_cin_bis.x=0; L_cin_bis.y=0; L_cin_bis.z=0;
    
    for(i=0;i<run.N+run.N2;i++)
    {
    // L_cin=somme_vect(L_cin,prod_vect(pos1_cart[i],omg1_cart[i]));
    L_cin_bis=somme_vect(L_cin_bis,omg1_cart[i]);
    }
   // printf("Lx=%e Ly=%e Lz=%e\n\n",L_cin_bis.x,L_cin_bis.y,L_cin_bis.z);
  
  
  
} // fin reinit_vit

//************************************


void initial_pot(double *pot, double *pot2, double *pot12, RUN run)
{
    //double epsilon=1.0; // facteur dans le potentiel
    const double vpotshift=4*/* *epsilon**/(pow(2.5,-12)-pow(2.5,-6));
    int i;
#ifdef POT
FILE *fpot;
char npot[]="potentiel.res";
if ((fpot = fopen (npot,"w")) == NULL)
    {
     fprintf(stderr,"Probleme d'ouverture de %s\n",npot);
     exit (1);  
    }
#endif

//if(run.flagpara==1)
//{
#pragma omp parallel for default(shared), private(i)
//}
    for (i=0;i<run.maxpot;i++)
    {     
      //double r=(i+0.5)*2.5*run.sigma/((double) run.maxpot);
      //******* définition de v_11(r) (seul potentiel si monodisperse)
      double r=(i+0.5)*2.0*run.sigma/((double) run.maxpot)+0.5*run.sigma;
      double r2=r*r/(run.sigma*run.sigma);
      double r6=r2*r2*r2;
      double r12=r6*r6;
      pot[i]=4/* *epsilon*/*(1.0/r12-1.0/r6)-vpotshift;
      
      //****** v_22(r)
      double s=(i+0.5)*2.0*run.sigma2/((double) run.maxpot)+0.5*run.sigma2;
      double s2=s*s/(run.sigma2*run.sigma2);
      double s6=s2*s2*s2;
      double s12=s6*s6;
      pot2[i]=4*(1.0/s12-1.0/s6)-vpotshift;
      
      //****** v_12(r)
      double u=(i+0.5)*2.0*run.sigma12/((double) run.maxpot)+0.5*run.sigma12;
      double u2=u*u/(run.sigma12*run.sigma12);
      double u6=u2*u2*u2;
      double u12=u6*u6;
      pot12[i]=4*(1.0/u12-1.0/u6)-vpotshift;
      

      
#ifdef POT
      fprintf(fpot,"%e %e %e %e %e %e\n",r,pot[i],s,pot2[i],u,pot12[i]);
#endif       
    } // fin boucle sur i
#ifdef POT    
      fclose(fpot);
#endif 
} // fini initial_pot

void initial_force(double *force_A,RUN run) // tabule les A(r_ij)
  {
   // const double vpotshift=4*(pow(2.5,-12)-pow(2.5,-6));
   
   //double epsilon=1.0;
   //const double force_shift=(1/sin(2.5))*12*epsilon/run.sigma*(2*pow(2.5,-13)-pow(2.5,-7));
   
    int i;
#ifdef POT
FILE *f_force;
char nforce[]="force_A.res";
if ((f_force = fopen (nforce,"w")) == NULL)
    { fprintf(stderr,"Probleme d'ouverture de %s\n",nforce);
      exit(1);  
    }
#endif
//if(run.flagpara==1)
//{
#pragma omp parallel for default(shared), private(i)
//}
    for (i=0;i<run.maxpot;i++)
    {
      double r=(i+0.5)*2.5/* *run.sigma*//((double) run.maxpot);
      double r2=r*r/*/(run.sigma*run.sigma)*/;
      double r7=r2*r2*r2*r;
      double r13=r7*r2*r2*r2;
      force_A[i]=(12/* *epsilon*//run.sigma)*(2.0/r13-1.0/r7)*1/(sin((r*run.sigma)))/*-force_shift*/; //A(r) (faut-il shifter la force ?)
#ifdef POT
      fprintf(f_force,"%e %e\n",r,force_A[i]);
#endif       
    }
#ifdef POT    
      fclose(f_force);
#endif 
  }
  
  //*****************************************************

inline double terme_A(double x,double *force_A,RUN run)
{ 
   double upper;
   double invupper;//tds
   int i;
   upper=2.5*run.sigma;
   invupper=run.maxpot/upper;
  if(x>upper) {return (0);}
  else { i=(int) (x*invupper); 
#ifdef DEBUG
  if (i==0)  
  { printf("distance %e %e\n",x,x*invupper);
  }
  if (i>=run.maxpot) {printf("%e %d\n",x,i);}
#endif
    return(force_A[i]); 
  }
}

/*! Fonction inversant un système matriciel 3x3 */

void inversion_systeme(POS_cart *omg, double a, double b, double c, double d, double e, double f, double L_cin_x, double L_cin_y, double L_cin_z)
{
  double a_data[] = { a, b, c, 
                      b, d, e,
                      c, e, f};
		      
  double b_data[] = { L_cin_x, L_cin_y, L_cin_z };
  
  gsl_matrix_view m 
         = gsl_matrix_view_array (a_data, 3, 3);
     
       gsl_vector_view b_vect
         = gsl_vector_view_array (b_data, 3);
     
	// gsl_vector *x=NULL;
	// allocdouble(&x,3);
       gsl_vector *x = gsl_vector_alloc (3);
       
       int s;
     
       gsl_permutation *p = gsl_permutation_alloc (3);
     
       gsl_linalg_LU_decomp (&m.matrix, p, &s);
     
       gsl_linalg_LU_solve (&m.matrix, p, &b_vect.vector, x);
    
       
       omg[0].x=x->data[0];
       omg[0].y=x->data[1];
       omg[0].z=x->data[2];
       
       //double q=x->data[0];
       //gsl_vector_fprintf (stdout, x, "%g");
       
       
       gsl_permutation_free (p);
       gsl_vector_free (x);
  
}

