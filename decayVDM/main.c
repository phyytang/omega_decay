/*====== Modules ===============
   Keys to switch on 
   various modules of micrOMEGAs  
================================*/
      
#define MASSES_INFO      
  /* Display information about mass spectrum  */
//#define OMEGA            
  /* Calculate relic density and display contribution of  individual channels */
#define INDIRECT_DETECTION  
  /* Compute spectra of gamma/positron/antiprotons/neutrinos for DM annihilation; 
     Calculate <sigma*v>;
     Integrate gamma signal over DM galactic squared density for given line 
     of sight; 
     Calculate galactic propagation of positrons and antiprotons.      
  */
      
/*#define RESET_FORMFACTORS*/
  /* Modify default nucleus form factors, 
    DM velocity distribution,
    A-dependence of Fermi-dencity
  */     
//#define CDM_NUCLEON     
  /* Calculate amplitudes and cross-sections for  CDM-mucleon collisions */  

/*#define CDM_NUCLEUS*/      
  /* Calculate number of events for 1kg*day and recoil energy distibution
      for various nuclei
  */
//#define NEUTRINO 

#define DECAYS
//#define CROSS_SECTIONS

  
/*===== end of Modules  ======*/

/*===== Options ========*/
/*#define SHOWPLOTS*/
     /* Display  graphical plots on the screen */ 

#define CLEAN  to clean intermediate files

/*===== End of DEFINE  settings ===== */


#include"../sources/micromegas.h"
#include"../sources/micromegas_aux.h"
#include"lib/pmodel.h"


int main(int argc,char** argv)
{  int err;
   char cdmName[10];
   int spin2, charge3,cdim;
   
  ForceUG=0;  /* to Force Unitary Gauge assign 1 */
  
  if(argc==1)
  { 
      printf(" Correct usage:  ./main  <file with parameters> \n");
      printf("Example: ./main data1.par\n");
      exit(1);
  }
                               
  err=readVar(argv[1]);
  
  if(err==-1)     { printf("Can not open the file\n"); exit(1);}
  else if(err>0)  { printf("Wrong file contents at line %d\n",err);exit(1);}
           
  

  err=sortOddParticles(cdmName);
  if(err) { printf("Can't calculate %s\n",cdmName); return 1;}
  
   qNumbers(cdmName, &spin2, &charge3, &cdim);
   printf("\nDark matter candidate is '%s' with spin=%d/2\n",
    cdmName,       spin2); 
   if(charge3) { printf("Dark Matter has electric charge %d*3\n",charge3); exit(1);}
   if(cdim!=1) { printf("Dark Matter ia a color particle\n"); exit(1);}
#ifdef MASSES_INFO
{
  printf("\n=== MASSES OF HIGG AND ODD PARTICLES: ===\n");
  printHiggs(stdout);
  printMasses(stdout,1);
}
#endif

#ifdef OMEGA
{ int fast=1;
  double Beps=1.E-5, cut=0.01;
  double Omega,Xf;   
  printf("\n==== Calculation of relic density =====\n");  
  Omega=darkOmega(&Xf,fast,Beps);
  printf("Xf=%.2e Omega=%.2e\n",Xf,Omega);
  printChannels(Xf,cut,Beps,1,stdout);   
}
#endif


#ifdef INDIRECT_DETECTION
{ 
  int err,i,ii;
  double Emin=1,/* Energy cut  in GeV   */  sigmaV;
  double vcs_gz,vcs_gg;
  char txt[100];
  double SpA[NZ],SpE[NZ],SpP[NZ];
  double FluxA[NZ],FluxE[NZ],FluxP[NZ];
  double * SpNe=NULL,*SpNm=NULL,*SpNl=NULL;
  double Etest=Mcdm/2;

  double lifetime=1.0E26;
  double scalefactor=1;

  printf("\n==== Indirect detection =======\n");  

  //  sigmaV=calcSpectrum(4,SpA,SpE,SpP,SpNe,SpNm,SpNl ,&err);
  sigmaV=calcSpectrum(6,SpA,SpE,SpP,SpNe,SpNm,SpNl ,&err);
    /* Returns sigma*v in cm^3/sec.     SpX - calculated spectra of annihilation.
       Use SpectdNdE(E, SpX) to calculate energy distribution in  1/GeV units.
       
       First parameter 1-includes W/Z polarization
                       2-includes gammas for 2->2+gamma
                       4-print cross sections             
    */

  /* *set an arbitary sigmaV, only affect the overall flux, but not the shape */
  /* *sigmaV=1.0E-24;*/
  DecayDM = 1;
  if(DecayDM==1)
    {
      scalefactor=1;
      lifetime=lifetime*scalefactor; // units [s]
      sigmaV=Mcdm/lifetime/rhoDM; // effective sigmaV
    }

  printf("lifetime=%.2E[s], sigmav=Mcdm/lifetime/rhoDM=%.2E[cm^3/s]\n",lifetime,sigmaV);  

  if(SpA)
  { 
    //     double fi=0.1,dfi=0.05; /* angle of sight and 1/2 of cone angle in [rad] */ 
     double fi=0.1,dfi=0.02; /* angle of sight and 1/2 of cone angle in [rad] */ 

     gammaFluxTab(fi,dfi, sigmaV, SpA,  FluxA);     
     printf("Photon flux  for angle of sight f=%.2f[rad]\n"
     "and spherical region described by cone with angle %.2f[rad]\n",fi,2*dfi);
#ifdef SHOWPLOTS
     sprintf(txt,"Photon flux[cm^2 s GeV]^{1} at f=%.2f[rad], cone angle %.2f[rad]",fi,2*dfi);
     displaySpectrum(FluxA,txt,Emin,Mcdm,1);
#endif
     printf("Photon flux = %.2E[cm^2 s GeV]^{-1} for E=%.1f[GeV]\n",SpectdNdE(Etest, FluxA), Etest);       

    if(DecayDM==1)
      {
	for(ii=0;ii<131;ii=ii+1)
	  {
	    //	Etest = Mcdm*pow(1.E-7,pow((double)(ii)/(double)(NZ),1.5));
	    //	Etest = (ii+1.0)*Mcdm/NZ;
	    Etest=1.0*pow(10.0,ii/40.0);
	    //	printf("%d\t %.2f\t %.2E\n",ii,Etest,SpectdNdE(Etest, FluxE)*pow(Etest,3.0));
	    //	fprintf(ofp,"%d\t %.4f\t %.2E\n",ii,Etest,SpectdNdE(Etest, FluxE)*pow(Etest,3.0));
	    //  fprintf(ofp,"%d\t  %.2E\n",ii,SpectdNdE(Etest, FluxA)*pow(Etest,2.0));
	    printf("%d\t %f\t  %.2E\n",ii,Etest,SpectdNdE(Etest, FluxA)*pow(Etest,2.0));
	  }
      }
  }

  char outputFilename[100] = "/home/ytang/test.dat";
  sprintf(outputFilename,"/home/ytang/decay_positron_%.0f.dat",2*Mcdm);

  FILE *ifp, *ofp;
  ofp = fopen(outputFilename, "w");

  if(SpE)
  { 
    double SMmev=550;// in MeV
    posiFluxTab(Emin, sigmaV, SpE,  FluxE);

    if(SMmev>0)  solarModulation(SMmev,0.0005,FluxE,FluxE);
#ifdef SHOWPLOTS     
    displaySpectrum(FluxE,"positron flux [cm^2 s sr GeV]^{-1}" ,Emin,Mcdm, 1);
#endif
    printf("Positron flux  =  %.2E[cm^2 sr s GeV]^{-1} for E=%.1f[GeV] \n",
    SpectdNdE(Etest, FluxE),  Etest);           

    if(DecayDM==1)
      {
	for(ii=0;ii<71;ii=ii+1)
	  {
	    //	Etest = Mcdm*pow(1.E-7,pow((double)(ii)/(double)(NZ),1.5));
	    //	Etest = (ii+1.0)*Mcdm/NZ;
	    Etest=1.0*pow(10.0,ii/20.0);
	    //	printf("%d\t %.2f\t %.2E\n",ii,Etest,SpectdNdE(Etest, FluxE)*pow(Etest,3.0));
	    //	fprintf(ofp,"%d\t %.4f\t %.2E\n",ii,Etest,SpectdNdE(Etest, FluxE)*pow(Etest,3.0));
	    fprintf(ofp,"%d\t  %.2E\n",ii,SpectdNdE(Etest, FluxE)*pow(Etest,3.0));
	  }
      }
  }
  
  fclose(ofp);

  if(SpP)
  { 
    double SMmev=550;
    pbarFluxTab(Emin, sigmaV, SpP,  FluxP  ); 

    if(SMmev>0)  solarModulation(SMmev,1,FluxP,FluxP);  
#ifdef SHOWPLOTS    
     displaySpectrum(FluxP,"antiproton flux [cm^2 s sr GeV]^{-1}" ,Emin, Mcdm,1);
#endif
    printf("Antiproton flux  =  %.2E[cm^2 sr s GeV]^{-1} for E=%.1f[GeV] \n",
    SpectdNdE(Etest, FluxP),  Etest);             
  }
}  
#endif

#ifdef RESET_FORMFACTORS
{
/* 
   The user has approach to form factors  which specifies quark contents 
   of  proton and nucleon via global parametes like
      <Type>FF<Nucleon><q>
   where <Type> can be "Scalar", "pVector", and "Sigma"; 
         <Nucleon>     "P" or "N" for proton and neutron
         <q>            "d", "u","s"

   calcScalarFF( Mu/Md, Ms/Md, sigmaPiN[MeV], sigma0[MeV])  
   calculates and rewrites Scalar form factors
*/
  printf("\n======== RESET_FORMFACTORS ======\n");
 
  printf("protonFF (default) d %.2E, u %.2E, s %.2E\n",ScalarFFPd, ScalarFFPu,ScalarFFPs);                               
  printf("neutronFF(default) d %.2E, u %.2E, s %.2E\n",ScalarFFNd, ScalarFFNu,ScalarFFNs);

  calcScalarFF(0.553,18.9,70.,35.);

  printf("protonFF (new)     d %.2E, u %.2E, s %.2E\n",ScalarFFPd, ScalarFFPu,ScalarFFPs);                               
  printf("neutronFF(new)     d %.2E, u %.2E, s %.2E\n",ScalarFFNd, ScalarFFNu,ScalarFFNs);

}
#endif

#ifdef CDM_NUCLEON
{ double pA0[2],pA5[2],nA0[2],nA5[2];
  double Nmass=0.939; /*nucleon mass*/
  double SCcoeff;        

printf("\n==== Calculation of CDM-nucleons amplitudes  =====\n");   

    nucleonAmplitudes(NULL, pA0,pA5,nA0,nA5);
    printf("CDM[antiCDM]-nucleon micrOMEGAs amplitudes:\n");
    printf("proton:  SI  %.3E [%.3E]  SD  %.3E [%.3E]\n",pA0[0], pA0[1],  pA5[0], pA5[1] );
    printf("neutron: SI  %.3E [%.3E]  SD  %.3E [%.3E]\n",nA0[0], nA0[1],  nA5[0], nA5[1] ); 

  SCcoeff=4/M_PI*3.8937966E8*pow(Nmass*Mcdm/(Nmass+ Mcdm),2.);
    printf("CDM[antiCDM]-nucleon cross sections[pb]:\n");
    printf(" proton  SI %.3E [%.3E] SD %.3E [%.3E]\n",
       SCcoeff*pA0[0]*pA0[0],SCcoeff*pA0[1]*pA0[1],3*SCcoeff*pA5[0]*pA5[0],3*SCcoeff*pA5[1]*pA5[1]);
    printf(" neutron SI %.3E [%.3E] SD %.3E [%.3E]\n",
       SCcoeff*nA0[0]*nA0[0],SCcoeff*nA0[1]*nA0[1],3*SCcoeff*nA5[0]*nA5[0],3*SCcoeff*nA5[1]*nA5[1]);

}
#endif
  
#ifdef CDM_NUCLEUS
{ double dNdE[300];
  double nEvents;

printf("\n======== Direct Detection ========\n");    

  nEvents=nucleusRecoil(Maxwell,73,Z_Ge,J_Ge73,SxxGe73,NULL,dNdE);

  printf("73Ge: Total number of events=%.2E /day/kg\n",nEvents);
  printf("Number of events in 10 - 50 KeV region=%.2E /day/kg\n",
                                   cutRecoilResult(dNdE,10,50));
                                                                                                         
#ifdef SHOWPLOTS
    displayRecoilPlot(dNdE,"Distribution of recoil energy of 73Ge",0,199);
#endif

  nEvents=nucleusRecoil(Maxwell,131,Z_Xe,J_Xe131,SxxXe131,NULL,dNdE);

  printf("131Xe: Total number of events=%.2E /day/kg\n",nEvents);
  printf("Number of events in 10 - 50 KeV region=%.2E /day/kg\n",
                                   cutRecoilResult(dNdE,10,50));                                   
#ifdef SHOWPLOTS
    displayRecoilPlot(dNdE,"Distribution of recoil energy of 131Xe",0,199);
#endif

  nEvents=nucleusRecoil(Maxwell,23,Z_Na,J_Na23,SxxNa23,NULL,dNdE);

  printf("23Na: Total number of events=%.2E /day/kg\n",nEvents);
  printf("Number of events in 10 - 50 KeV region=%.2E /day/kg\n",
                                   cutRecoilResult(dNdE,10,50));                                   
#ifdef SHOWPLOTS
    displayRecoilPlot(dNdE,"Distribution of recoil energy of 23Na",0,199);
#endif

  nEvents=nucleusRecoil(Maxwell,127,Z_I,J_I127,SxxI127,NULL,dNdE);

  printf("I127: Total number of events=%.2E /day/kg\n",nEvents);
  printf("Number of events in 10 - 50 KeV region=%.2E /day/kg\n",
                                   cutRecoilResult(dNdE,10,50));                                   
#ifdef SHOWPLOTS
    displayRecoilPlot(dNdE,"Distribution of recoil energy of 127I",0,199);
#endif
  
}
#endif 

#ifdef NEUTRINO
{ double nu[NZ], nu_bar[NZ],mu[NZ];
  double Ntot;
  int forSun=1;
  double Emin=1;
  
 printf("\n===============Neutrino Telescope=======  for  "); 
 if(forSun) printf("Sun\n"); else printf("Earth\n");  

  err=neutrinoFlux(Maxwell,forSun, nu,nu_bar);
#ifdef SHOWPLOTS
  displaySpectrum(nu,"nu flux from Sun [1/Year/km^2/GeV]",Emin,Mcdm,1);
  displaySpectrum(nu_bar,"nu-bar from Sun [1/Year/km^2/GeV]",Emin,Mcdm,1);
#endif
{ double Ntot;
  double Emin=1; //GeV
  spectrInfo(Emin/Mcdm,nu, &Ntot,NULL);
    printf(" E>%.1E GeV neutrino flux       %E [1/Year/km^2] \n",Emin,Ntot);
  spectrInfo(Emin/Mcdm,nu_bar, &Ntot,NULL);
    printf(" E>%.1E GeV anti-neutrino flux  %E [1/Year/km^2]\n",Emin,Ntot);  
} 
  
/* Upward events */
  
  muonUpward(nu,nu_bar, mu);
#ifdef SHOWPLOTS  
  displaySpectrum(mu,"Upward muons[1/Year/km^2/GeV]",1,Mcdm/2,1);
#endif
  { double Ntot;
    double Emin=1; //GeV
    spectrInfo(Emin/Mcdm,mu, &Ntot,NULL);
    printf(" E>%.1E GeV Upward muon flux    %E [1/Year/km^2]\n",Emin,Ntot);
  } 
  
/* Contained events */
  muonContained(nu,nu_bar,1., mu);
#ifdef SHOWPLOTS  
  displaySpectrum(mu,"Contained  muons[1/Year/km^3/GeV]",Emin,Mcdm,1); 
#endif
  { double Ntot;
    double Emin=1; //GeV
    spectrInfo(Emin/Mcdm,mu, &Ntot,NULL);
    printf(" E>%.1E GeV Contained muon flux %E [1/Year/km^3]\n",Emin,Ntot);
  }  
}        
#endif 


#ifdef DECAYS
{ 
  txtList LZ,Lh;
   double width,br;
   char * pname;
   
   pname=pdg2name(23);
   if(pname) 
   {  width=pWidth(pname,&LZ);
      printf("%s :   total width=%E \n and Branchings:\n",pname,width);
      printTxtList(LZ,stdout);
   } else printf("Z particle is detected in the model\n");

   pname=pdg2name(112);
   if(pname) 
   {  width=pWidth(pname,&LZ);
      printf("%s :   total width=%E \n and Branchings:\n",pname,width);
      printTxtList(LZ,stdout);
   } else printf("X particle is detected in the model\n");

}
#endif

#ifdef CROSS_SECTIONS
{
   char * pname1,*pname2;
   char process[30];
   double Pcm=250,cosmin=-0.99, cosmax=0.99;
   
   pname1= pdg2name(11);
   pname2= pdg2name(-11);
   if(pname1 && pname2)
   { numout*cc;
     sprintf(process,"%s,%s->2*x",pname1,pname2);
     cc=newProcess(process); 
     if(cc)
     {  int ntot,l;
        char * name[4];
        procInfo1(cc,&ntot,NULL,NULL);
        for(l=1;l<=ntot; l++)
        { int err;
          double cs;
          char txt[100];
          procInfo2(cc,l,name,NULL);
          sprintf(txt,"%3s,%3s -> %3s %3s  ",name[0],name[1],name[2],name[3]);
          cs= cs22(cc,l,Pcm,cosmin,cosmax,&err);
          if(err) printf("%-20.20s    Error\n",txt);
          else if(cs) printf("%-20.20s  %.2E [pb]\n",txt,cs); 
        }
     } 
   } 
}

#endif

#ifdef CLEAN

#endif

  killPlots();
  return 0;
}
