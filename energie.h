/*! Calcul potentiel d'interaction */

inline double poten(double x, double *pot, RUN run)
{ 
   double upper; 
   double invupper;
   long i;
   /*upper=2.5*run.sigma;
   invupper=run.maxpot/upper;
   */
   //***** pour tabulation entre 0.5*sigma et 2.5*sigma
   upper=2.0*run.sigma; 
   invupper=run.maxpot/upper;
  
  if(x>upper) {return 0.0;}
  else
  {
   
   i=(int)((x-0.5*run.sigma)*invupper); // pour tabulation entre 0.5*run.sigma et 2.5*run.sigma
  
#ifdef DEBUG
  if (i==0)  
  { printf("distance %e %e\n",x,x*invupper);
  }
  if (i>=run.maxpot) {printf("%e %d\n",x,i);}
#endif
  if(i<0)
  {
    return(0.0);
  }
  else
  {
    return(pot[i]);
  }
  }
}

/*! Calcul potentiel d'interaction v22(r) */

inline double poten2(double x, double *pot2, RUN run)
{ 
  double upper2;
  double invupper2;
    long i;
   upper2=2.0*run.sigma2;;
   invupper2=run.maxpot/upper2;
  
  if(x>upper2) {return 0.0;}
  else
  {  
   i=(int)((x-0.5*run.sigma2)*invupper2); // pour tabulation entre 0.5*run.sigma et 2.5*run.sigma
  
#ifdef DEBUG
  if (i==0)  
  {
   printf("distance %e %e\n",x,x*invupper2);
  }
  if (i>=run.maxpot) {printf("%e %d\n",x,i);}
#endif
  if(i<0)
  {
   return(0.0);
  }
  else
  {
   return(pot2[i]);
  }
  }
}

/*! Calcul potentiel d'interaction v12(r) */

inline double poten12(double x, double *pot12, RUN run)
{ 
  double upper12;
  double invupper12;
    long i;
   upper12=2.0*run.sigma12;;
   invupper12=run.maxpot/upper12;
  
  if(x>upper12) {return 0.0;}
  else
  {  
   i=(int)((x-0.5*run.sigma12)*invupper12); // pour tabulation entre 0.5*run.sigma et 2.5*run.sigma
  
#ifdef DEBUG
  if (i==0)  
  {
   printf("distance %e %e\n",x,x*invupper12);
  }
  if (i>=run.maxpot) {printf("%e %d\n",x,i);}
#endif
  if(i<0)
  {
   return(0.0);
  }
  else
  {
   return(pot12[i]);
  }
  }
}




/*! Calcul énergie potentielle */

inline double calcul_nrj_cart(POS_cart *pos1_cart, double *pot, double *pot2, double *pot12, double *arccos, RUN run)
{
  int i,j;
double Etot=0.0;
int nrun=run.N;
int nrun2=run.N2;
  const double prod_scal_max=cos(2.5*run.sigma);
  const double prod_scal_max2=cos(2.5*run.sigma2);
  const double prod_scal_max12=cos(2.5*run.sigma12);
  
#pragma omp parallel for default(shared), private(i,j), reduction(+:Etot)
for (i=0;i<nrun-1;i++)
  for(j=i+1;j<nrun;j++)
  { 
   // if(j!=i)
   // {
      double prod_scal=prod_scalaire(pos1_cart[i],pos1_cart[j]);
      if(prod_scal>prod_scal_max)
      {
    double x=acos(prod_scal);
   // if (x<0.7*run.sigma) printf("%e %d %d\n",x,i,j);
    Etot+=poten(x,pot,run);
      }
   // }
  }
  
 // if(nrun2>0)
 // {

  
#pragma omp parallel for default(shared), private(i,j), reduction(+:Etot)
for (i=nrun;i<nrun+nrun2-1;i++)
  for(j=i+1;j<nrun+nrun2;j++)
  { 
//     if(j!=i)
//     {
      double prod_scal=prod_scalaire(pos1_cart[i],pos1_cart[j]);
      if(prod_scal>prod_scal_max2)
      {
    double x=acos(prod_scal);
   // if (x<0.7*run.sigma) printf("%e %d %d\n",x,i,j);
    Etot+=poten2(x,pot2,run);
      }
  //  }
  }
  
#pragma omp parallel for default(shared), private(i,j), reduction(+:Etot)
for (i=0;i<nrun;i++)
  for(j=nrun;j<nrun+nrun2;j++)
  { 
   // if(j!=i)
   // {
      double prod_scal=prod_scalaire(pos1_cart[i],pos1_cart[j]);
      if(prod_scal>prod_scal_max12)
      {
    double x=acos(prod_scal);
   // if (x<0.7*run.sigma) printf("%e %d %d\n",x,i,j);
    Etot+=poten12(x,pot12,run);
      }
   // }
  }
  
return(Etot);
} // fin calcul_nrj_cart
 //************************************************
//fprintf(fspectre_energpot, "%e\t %e\n",1/((double) nfrac*run.dt*dt_energpot)*run.sigma*kf, spectre_nrj_pot[kf]); 





/*! Calcul énergie potentielle et son spectre */

void calcul_nrj_pot(int flag, int step, POS_cart *pos2_cart, double *pot, double *pot2, double *pot12, THERMO *thermo, double *arccos, RUN run)
{
 int i, j, k, t, ndefauts, f;
 static double *nrj_pot=NULL;
 static double *spectre_nrj_pot=NULL;
 int nstep=(int) (run.nstep);
 int nfrac=20000; // nombre de points pour frac6 et spectre_frac6
 
 if(flag==0)
 {
  
  allocdouble(&nrj_pot,nfrac);
  allocdouble(&spectre_nrj_pot,nfrac);
  t=0;
  //  printf("intialisation calcul_defauts OK\n");
 }
 
 if(flag==1)
 {
  t=step/run.moviestep;

  nrj_pot[t]=1/((double) (run.N+run.N2))*calcul_nrj_cart(pos2_cart,pot,pot2,pot12,arccos,run);
  
  
 //***************** test BUg 
//  if(step%1000==0)
//  {
//  int dt_energpot=nstep/nfrac;
//  printf("nfrac %d\n", nfrac);
//  printf("dt_energpot %d\n", dt_energpot);
//  printf("sigma %e\n", run.sigma);
//  printf("dt %e\n", run.dt);
//  printf("freq_max %e\n",1/(0.0004*1)*0.34);
//  printf("freq_max %e\n",run.sigma/(run.dt*dt_energpot));
//  }
  
  
 } // fin if flag == 1
 
 
 if(flag==2)
 {
   

   
 
/*! Séries temporelles des fractions de 6 et défauts */ 
 
 FILE *ffonction_energpot;
  char energpot[50];
 if(run.flagposanc==0)
  {

    sprintf(energpot,"energpotNa%dNb%deta%dTemp%dbstep%dstep%dtFs%d.res",run.N,run.N2,(int) (1000*run.eta),(int) (1000/run.beta),(int) run.bstep,(int) run.nstep,(int) (run.tFs));
  }
  else  //sinon le nouveau bstep est bstep+step_oldsimu
  {
    sprintf(energpot,"energpotNa%dNb%deta%dTemp%dbstep%dstep%dtFs%d.res",run.N,run.N2,(int) (1000*run.eta),(int) (1000/run.beta),(int) (run.bstep+run.bstep_old+run.step_old),(int) (run.nstep),(int) (run.tFs));
  }
 
    
  if ((ffonction_energpot = fopen (energpot,"w")) == NULL)
    {
     fprintf(stderr,"Probleme d'ouverture de %s\n",energpot);
     exit (1);  
    }
    
    
  int dt_energpot=nstep/nfrac;
  double delta_f=1/((double) nfrac*dt_energpot);
  
//   printf("delta_f %e\n", delta_f);
//   printf("dt_energpot %d\n", dt_energpot);
//   printf("sigma %e\n", run.sigma);
//   printf("dt %e\n", run.dt);
//   printf("freq_max %e\n",1/(0.0004*1)*0.34);
//   printf("freq_max %e\n",run.sigma/(run.dt*dt_energpot));
//   printf("freq_max %e\n\n", nfrac*(delta_f/run.dt)*run.sigma);
  
 for(t=0;t<nfrac;t++)
 {
  fprintf(ffonction_energpot, "%e\t %e\n",t*run.dt/run.sigma*dt_energpot, nrj_pot[t]); 
 }
 fclose(ffonction_energpot); 
 
 //********************* TFourier
 
// int delta_t=nfrac/run.N;
 //int delta_t=dt_energpot;
 
 
 
 
 
 complexe sum;
 complexe i0; i0.Re=0; i0.Im=1;
 int kt, kf;
 
 for(kf=0;kf<nfrac;kf++) 
 {
  //double f=pi2/(run.N*delta_t)*i;
 // sum.Re=0; sum.Im=0;
 // for(kt=0;kt<nfrac;kt++) sum=somme_complexe(sum,scalaire_complexe(1/((double) nfrac)*nrj_pot[kt],exp_complexe(scalaire_complexe(-pi2*kt*kf/nfrac,i0))));
 //  spectre_nrj_pot[kf]=dt_energpot*dt_energpot*module_complexe(sum)*module_complexe(sum);
 // spectre_nrj_pot[kf]=module_complexe(sum)*module_complexe(sum);
  double sum_cos, sum_sin;
  
  for(kt=0;kt<nfrac;kt++)
  {
   sum_cos+=nrj_pot[kt]*cos(pi2*kt*kf/nfrac);
   sum_sin+=nrj_pot[kt]*sin(pi2*kt*kf/nfrac);
  }
  spectre_nrj_pot[kf]=1/((double) nfrac*nfrac)*(sum_cos*sum_cos+sum_sin*sum_sin);
   
 }

 
  FILE *fspectre_energpot;
  char spectre_energpot[50];
 if(run.flagposanc==0)
  {

    sprintf(spectre_energpot,"spectre_energpotNa%dNb%deta%dTemp%dbstep%dstep%dtFs%d.res",run.N,run.N2,(int) (1000*run.eta),(int) (1000/run.beta),(int) run.bstep,(int) run.nstep,(int) (run.tFs));
  }
  else  //sinon le nouveau bstep est bstep+step_oldsimu
  {
    sprintf(spectre_energpot,"spectre_energpotNa%dNb%deta%dTemp%dbstep%dstep%dtFs%d.res",run.N,run.N2,(int) (1000*run.eta),(int) (1000/run.beta),(int) (run.bstep+run.bstep_old+run.step_old),(int) (run.nstep),(int) (run.tFs));
  }
 
    
  if ((fspectre_energpot = fopen (spectre_energpot,"w")) == NULL)
    {
     fprintf(stderr,"Probleme d'ouverture de %s\n",spectre_energpot);
     exit (1);  
    }
    
    nstep=(int) (run.nstep);
    nfrac=20000; 
    dt_energpot=nstep/nfrac;
    delta_f=1/((double) nfrac*dt_energpot);
  
  //****************** test bug
    printf("nfrac %d\n",nfrac);
    printf("delta_f %e\n", delta_f);
    printf("dt_energpot %d\n", dt_energpot);
    printf("sigma %e\n", run.sigma);
    printf("dt %e\n", run.dt);
    printf("freq_max %e\n", nfrac*(delta_f/run.dt)*run.sigma);
 // int dt_defauts=run.nstep/1000;
 for(kf=1;kf<nfrac;kf++)
 {
 // fprintf(fspectre_energpot, "%e\t %e\n",1/((double) nfrac*run.dt*dt_energpot)*run.sigma*kf, spectre_nrj_pot[kf]); 
//  fprintf(fspectre_energpot, "%e\t %e\n",(delta_f/run.dt)*run.sigma*kf, spectre_nrj_pot[kf]);
   fprintf(fspectre_energpot, "%e\t %e\n",run.sigma/(nfrac*run.dt)*kf, spectre_nrj_pot[kf]);
 }
 fclose(fspectre_energpot); 
 

 
 
 } // fin if flag == 2 
 
} // fin calcul_nrj_pot
 
 
 
/*! Calcul énergie potentielle et son spectre lissé */

void calcul_nrj_pot_lissee(int flag, int step, POS_cart *pos2_cart, double *pot, double *pot2, double *pot12, THERMO *thermo, double *arccos, RUN run)
{
 int i, j, k, t, ndefauts, f;
 static double *nrj_pot=NULL;
 static double *spectre_nrj_pot=NULL;
 int nstep=(int) (run.nstep);
 int nlissage=20; // nbre de points sur lesquels on lisse energ_pot
 int nfrac=20000/nlissage; // nombre de points pour energ_pot et spectre_energ_pot
 
 if(flag==0)
 {
  allocdouble(&nrj_pot,nfrac);
  allocdouble(&spectre_nrj_pot,nfrac);
  t=0;
 }
 
 if(flag==1)
 {
  t=step/run.moviestep;
  nrj_pot[t]=1/((double) (run.N+run.N2))*calcul_nrj_cart(pos2_cart,pot,pot2,pot12,arccos,run); 
 } // fin if flag == 1
 
 
if(flag==2)
{

 /*! Série temporelle de l'énergie potentielle */
 
FILE *ffonction_energpot;
char energpot[50];
if(run.flagposanc==0)
{
 sprintf(energpot,"energpot_lisseeNa%dNb%deta%dTemp%dbstep%dstep%dtFs%d.res",run.N,run.N2,(int) (1000*run.eta),(int) (1000/run.beta),(int) run.bstep,(int) run.nstep,(int) (run.tFs));
}
else  //sinon le nouveau bstep est bstep+step_oldsimu
{
 sprintf(energpot,"energpot_lisseeNa%dNb%deta%dTemp%dbstep%dstep%dtFs%d.res",run.N,run.N2,(int) (1000*run.eta),(int) (1000/run.beta),(int) (run.bstep+run.bstep_old+run.step_old),(int) (run.nstep),(int) (run.tFs));
}
 
    
if ((ffonction_energpot = fopen (energpot,"w")) == NULL)
{
 fprintf(stderr,"Probleme d'ouverture de %s\n",energpot);
 exit (1);  
}
    
    
int dt_energpot=nstep/nfrac;
double delta_f=1/((double) nfrac*dt_energpot);
for(t=0;t<nfrac;t++)
{
 fprintf(ffonction_energpot, "%e\t %e\n",t*run.dt/run.sigma*dt_energpot, nrj_pot[t]); 
}
fclose(ffonction_energpot); 
 
 //********************* TFourier

 int kt, kf;
 
for(kf=0;kf<nfrac;kf++) 
{
 double sum_cos, sum_sin;  
 for(kt=0;kt<nfrac;kt++)
 {
  sum_cos+=nrj_pot[kt]*cos(pi2*kt*kf/nfrac);
  sum_sin+=nrj_pot[kt]*sin(pi2*kt*kf/nfrac);
 }
 spectre_nrj_pot[kf]=1/((double) nfrac*nfrac)*(sum_cos*sum_cos+sum_sin*sum_sin);   
}
 
FILE *fspectre_energpot;
char spectre_energpot[50];
if(run.flagposanc==0)
{
 sprintf(spectre_energpot,"spectre_energpot_lisseeNa%dNb%deta%dTemp%dbstep%dstep%dtFs%d.res",run.N,run.N2,(int) (1000*run.eta),(int) (1000/run.beta),(int) run.bstep,(int) run.nstep,(int) (run.tFs));
}
else  //sinon le nouveau bstep est bstep+step_oldsimu
{
 sprintf(spectre_energpot,"spectre_energpot_lisseeNa%dNb%deta%dTemp%dbstep%dstep%dtFs%d.res",run.N,run.N2,(int) (1000*run.eta),(int) (1000/run.beta),(int) (run.bstep+run.bstep_old+run.step_old),(int) (run.nstep),(int) (run.tFs));
}
 
    
if ((fspectre_energpot = fopen (spectre_energpot,"w")) == NULL)
{
 fprintf(stderr,"Probleme d'ouverture de %s\n",spectre_energpot);
 exit (1);  
}
    
nstep=(int) (run.nstep);
nfrac=20000/nlissage; 
dt_energpot=nstep/nfrac;
delta_f=1/((double) nfrac*dt_energpot);
  
  //****************** test bug
    printf("nfrac %d\n",nfrac);
    printf("delta_f %e\n", delta_f);
    printf("dt_energpot %d\n", dt_energpot);
    printf("sigma %e\n", run.sigma);
    printf("dt %e\n", run.dt);
    printf("freq_max %e\n", nfrac*(delta_f/run.dt)*run.sigma);
 // int dt_defauts=run.nstep/1000;
for(kf=1;kf<nfrac;kf++)
{
 // fprintf(fspectre_energpot, "%e\t %e\n",(delta_f/run.dt)*run.sigma*kf, spectre_nrj_pot[kf]);
 fprintf(fspectre_energpot, "%e\t %e\n",run.sigma/(nfrac*run.dt)*kf, spectre_nrj_pot[kf]);
}
fclose(fspectre_energpot); 
 

 
 
} // fin if flag == 2 
 
} // fin calcul_nrj_pot_lissee 
 
 
 
 
 