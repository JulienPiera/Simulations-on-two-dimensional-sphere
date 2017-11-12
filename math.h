void legendre_tabu(double *legendre,RUN run)
 {
   //int maxacos=10000;
   int i; 
   int k_wave=(int)(run.kFs*pi/run.sigma);
   if(run.flagkwave==0)
   {
    k_wave=(int)(run.kFs*pi/run.sigma); // k~1/sigma
   }
   else
   {
    k_wave=(int)(run.kFs*pi); // k~1/R
   }
   printf("kwave %d\n\n",k_wave);
  // int k_wave=15; //test

   /*
FILE *flegendre;
char nlegendre[]="legendrecos.res";
if ((flegendre = fopen (nlegendre,"w")) == NULL)
    { fprintf(stderr,"Probleme d'ouverture de %s\n",nlegendre);
      exit (1);  
    }
*/
   double startwtime,endwtime;
//if(run.flagpara==1)
//{
#pragma omp parallel  
//}  
  startwtime =  omp_get_wtime();
  
//if(run.flagpara==1)
//{
#pragma omp parallel for  default(shared), private(i)
//}
    for (i=0;i<run.maxlegendre;i++)
    {
       double r=pi*(i+0.5)/((double) run.maxlegendre);
    legendre[i]=gsl_sf_legendre_Pl(k_wave,cos(r)); //test tabulation de legendre(cos())
  //  fprintf(flegendre, "%e %e\n",r,legendre[i]);
    }
  endwtime =  omp_get_wtime(); 
//if(run.flagpara==1)
{
#pragma omp barrier

if (omp_get_thread_num() == 0) printf("legendre dynam  = %e\n",(double)(endwtime-startwtime));
}   
//fclose(flegendre);
 }

 // Tabulation de arccosinus
 
 
 //****************************************
 
 
void arccos_tabu(double *arccos,RUN run)
 {
   //int maxacos=10000;
   int i;
   
double startwtime,endwtime;

//if(run.flagpara==1)
//{
#pragma omp parallel
if (omp_get_thread_num() == 0)    startwtime =  omp_get_wtime();
#pragma omp parallel for default(shared), private(i)
//}
    for (i=0;i<run.maxacos;i++)
    {
    double r=-1+2*(i+0.5)/((double) run.maxacos);
   // double r=min_acos+(max_acos-min_acos)*i/((double) run.maxacos); // pour tabulation entre cos(2.5*sigma) et cos(0.5*sigma)
    arccos[i]=acos(r);
    }
//if(run.flagpara==1)
//{
if (omp_get_thread_num() == 0)    endwtime =  omp_get_wtime(); 
if (omp_get_thread_num() == 0) printf("arccos dynam  = %e\n",(double)(endwtime-startwtime));
//}
 }
 
 //************************************************************
 
 inline double arccosinus(double x, double *arccos, RUN run)
 {
 
   int i;
   i=(int) ((x+1)*run.maxacos/2-0.5);
  
   if((i>run.maxacos)||(i<0))
   {
   printf("i=%d\n",i);
   }
   return(arccos[i]);
 }
 
 inline double calcul_legendre(double x, double *legendre, RUN run)
 {
   int i;
   i=(int) (x*run.maxlegendre/pi-0.5); //tabulation de legendre(cos())
   return(legendre[i]);
 }
  
 

inline double distance(POS posi,POS posj, double *arccos, RUN run) /*! distance en coordonnées cylindriques */
{
  return(arccosinus(posi.z*posj.z+sqrt((1-posi.z*posi.z)*(1-posj.z*posj.z))*cos(posi.phi-posj.phi),arccos,run));
}

inline double distance_ang(POS_ang posi,POS_ang posj, double *arccos, RUN run) /*! distance en coordonnées sphériques */
{
  return(arccosinus(cos(posi.theta)*cos(posj.theta)+sin(posi.theta)*sin(posj.theta)*cos(posi.phi-posj.phi),arccos,run));
}

 double distance_cart(POS_cart posi,POS_cart posj, double *arccos, RUN run) /*! distance en coordonnées cartésiennes */
{
  double alpha=posi.x*posj.x+posi.y*posj.y+posi.z*posj.z;
  
  double alpha_min;
  double alpha_max;
  
  if(run.flagposanc==0)
  {
   alpha_min=-1.0-1e-6;
   alpha_max=1.0+1e-6;
  }
  else
  {
   alpha_min=-1.0-1e-6;
   alpha_max=1.0+1e-6;
  }
  
 // printf("alpha %e\n",alpha);
  if((alpha>=-1.0)&&(alpha<=1.0))
  {
  return(arccosinus(alpha,arccos,run)); // distance geodesique en coordonnees cartesiennes
  }
  else // on traite a  part le cas ou l'argument sort de l'intervalle de definition de arccosinus
  {
    if((alpha<-1.0)&&(alpha>alpha_min)) // si -1-10^-5 < alpha < -1
    {
     return(pi);
    }
    else
    {
    if((alpha>1.0)&&(alpha<alpha_max)) // si 1 < alpha < 1+10^-5
    {
     return(0.0);
    }
    else
    {   
      printf("alpha (bug) %e \n",alpha);
      
      printf("posi.x %e posi.y %e posi.z %e posj.x %e posj.y %e pos.z %e\n",posi.x,posi.y,posi.z,posj.x,posj.y,posj.z);
      printf("Le programme s'arrÃªte Ã  cause du calcul d'un arccosinus.\n");
      exit(1);
      return(0);
    }
    }
  }
}

/*! calcul des vitesses a  t a  partir de r_i(t-dt) et r_i(t+dt) */

inline double vitesse(int i, POS_ang *pos1, POS_ang *pos2, POS_ang *pos3, double *arccos, RUN run)
{
// double theta_apres=distance(pos1[i],pos3[i],arccos,run); // theta de i Ã  t+dt dans le rÃ©fÃ©rentiel oÃ¹ i Ã  t-dt est au sommet
 //double vit=theta_apres/(2*run.dt); //composante de la vitesse
 
 double vit_theta=1/(2*run.dt)*(pos3[i].theta-pos1[i].theta);
 double vit_phi=1/(2*run.dt)*(pos3[i].phi-pos1[i].phi)*sin(pos2[i].theta); //v,phi=sin(theta)*d(phi,i)/dt
 double vit=sqrt(vit_phi*vit_phi+vit_theta*vit_theta);
 
 return(vit);
}

inline POS_cart scalaire(double lambda, POS_cart a) /*! multiplication d'un vecteur par un scalaire */
{
 POS_cart c;
 c.x=lambda*a.x;
 c.y=lambda*a.y;
 c.z=lambda*a.z;
 
 return(c);
}
    
inline double prod_scalaire(POS_cart a, POS_cart b) /*! produit scalaire */
{
 double alpha=a.x*b.x+a.y*b.y+a.z*b.z;
 return(alpha);
}
  
inline POS_cart prod_vect(POS_cart a, POS_cart b) /*! produit vectoriel */
{
 POS_cart c;
 c.x=a.y*b.z-a.z*b.y;
 c.y=a.z*b.x-a.x*b.z;
 c.z=a.x*b.y-a.y*b.x;
  return(c);
}

inline POS_cart somme_vect(POS_cart a, POS_cart b) /*! somme de 2 vecteurs */
{
 POS_cart c;
 c.x=a.x+b.x;
 c.y=a.y+b.y;
 c.z=a.z+b.z;
 
 return(c);
}

inline complexe somme_complexe(complexe a, complexe b) /*! somme de 2 complexes */
{
 complexe c;
 c.Re=a.Re+b.Re;
 c.Im=a.Im+b.Im;
 return(c); 
}

inline complexe scalaire_complexe(double lambda, complexe a) /*! multiplication d'un complexe par un scalaire */
{
 complexe c;
 c.Re=lambda*a.Re;
 c.Im=lambda*a.Im;
 return(c); 
}

inline complexe prod_complexe(complexe a, complexe b) /*! produit */
{
 complexe c;
 c.Re=a.Re*b.Re-a.Im*b.Im;
 c.Im=a.Re*b.Im+a.Im*b.Re;
 return(c); 
}

inline complexe conjug_complexe(complexe a) /*! conjugué complexe */
{
 complexe c;
 c.Re=a.Re;
 c.Im=-a.Im; 
 return(c);
}

inline double module_complexe(complexe a) /*! module */
{
 double c=sqrt(a.Re*a.Re+a.Im*a.Im);
 return(c); 
}

inline complexe exp_complexe(complexe a)
{
 complexe i;
 i.Re=0;
 i.Im=1;
 complexe c;
 complexe cos_im_a; // cos(Im(a)), défini comme complexe
 cos_im_a.Re=cos(a.Im); cos_im_a.Im=0;
 c=scalaire_complexe(exp(a.Re),somme_complexe(cos_im_a,scalaire_complexe(sin(a.Im),i)));
 return(c);
}

//********************* fonction calculant la TF d'un vecteur

// void fourier(double *a, double *tf_a, int npts, RUN run)
// {
//   double delta_t=nfrac/run.N;
//  double delta_f=1/(run.N*delta_t);
//  complexe sum;
//  complexe i0; i0.Re=0; i0.Im=1;
//  int kt, kf;
//  
//  for(kf=0;kf<nfrac;kf++) 
//  {
//   //double f=pi2/(run.N*delta_t)*i;
//   sum.Re=0; sum.Im=0;
//   for(kt=0;kt<nfrac;kt++) sum=somme_complexe(sum,scalaire_complexe(frac6_fonc_t[kt],exp_complexe(scalaire_complexe(-pi2*kt*kf/run.N,i0))));
//    spectre_frac6[kf]=delta_t*delta_t*module_complexe(sum)*module_complexe(sum);
//  }
// }

//********************* fonction comptant le nombre de lignes dans un fichier

int line_count(FILE *n)
{
int c; 
int lines = 0;

while ((c = fgetc(n)) != EOF)
{
if (c == '\n') ++lines;
}

return lines+1;
}

