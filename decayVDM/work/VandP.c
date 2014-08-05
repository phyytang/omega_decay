#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "../../CalcHEP_src/include/extern.h"
#include "../../CalcHEP_src/include/VandP.h"
#include "autoprot.h"
extern int  FError;
/*  Special model functions  */

int nModelParticles=20;
static ModelPrtclsStr ModelPrtcls_[20]=
{
  {"G","G", 21, "0","0",2,8,0}
, {"A","A", 22, "0","0",2,1,0}
, {"Z","Z", 23, "MZ","wZ",2,1,0}
, {"W+","W-", 24, "MW","wW",2,1,3}
, {"h","h", 25, "Mh","wh",0,1,0}
, {"e","E", 11, "0","0",1,1,-3}
, {"ne","Ne", 12, "0","0",1,1,0}
, {"m","M", 13, "Mm","0",1,1,-3}
, {"nm","Nm", 14, "0","0",1,1,0}
, {"l","L", 15, "Ml","0",1,1,-3}
, {"nl","Nl", 16, "0","0",1,1,0}
, {"d'","D'", 81, "0","0",1,3,-1}
, {"u","U", 2, "0","0",1,3,2}
, {"s'","S'", 83, "0","0",1,3,-1}
, {"c","C", 4, "Mc","0",1,3,2}
, {"b","B", 5, "Mb","0",1,3,-1}
, {"t","T", 6, "Mtp","wt",1,3,2}
, {"~dk","~Dk", 111, "MDk","0",1,1,0}
, {"Xm","Xm", 112, "MXm","0",2,1,0}
, {"phi","phi", 113, "Mphi","0",0,1,0}
};
ModelPrtclsStr *ModelPrtcls=ModelPrtcls_; 
int nModelVars=21;
int nModelFunc=6;
static char*varNames_[27]={
 "EE","alfSMZ","SW","Mm","Ml","Q","McMc","MbMb","Mtp","MZ"
,"Mh","wt","wZ","wW","MDk","MXm","Mphi","ye","ym","yl"
,"yA","CW","MW","qcdOk","Mb","Mt","Mc"};
char**varNames=varNames_;
static REAL varValues_[27]={
   3.122300E-01,  1.172000E-01,  4.810000E-01,  1.057000E-01,  1.777000E+00,  1.000000E+02,  1.200000E+00,  4.230000E+00,  1.730000E+02,  9.118840E+01
,  2.000000E+02,  1.590000E+00,  2.494440E+00,  2.088950E+00,  5.005000E+02,  1.000000E+03,  1.000000E+00,  1.000000E+00,  1.000000E+00,  1.000000E+00
,  0.000000E+00};
REAL*varValues=varValues_;
int calcMainFunc(void)
{
   int i;
   static REAL * VV=NULL;
   static int iQ=-1;
   static int cErr=1;
   REAL *V=varValues;
   FError=0;
   if(VV && cErr==0)
   { for(i=0;i<nModelVars;i++) if(i!=iQ && VV[i]!=V[i]) break;
     if(i==nModelVars)      {if(iQ>=0 && VV[iQ]!=V[iQ]) goto FirstQ; else return 0;} 
   }
  cErr=1;
   V[21]=sqrt(1-pow(V[2],2));
   if(!finite(V[21]) || FError) return 21;
   V[22]=V[9]*V[21];
   if(!finite(V[22]) || FError) return 22;
   V[23]=initQCD(V[1],V[6],V[7],V[8]);
   if(!finite(V[23]) || FError) return 23;
 FirstQ:
 cErr=1;
   V[24]=MbEff(V[5])*1;
   if(!finite(V[24]) || FError) return 24;
   V[25]=MtEff(V[5])*1;
   if(!finite(V[25]) || FError) return 25;
   V[26]=McEff(V[5])*1;
   if(!finite(V[26]) || FError) return 26;
   if(VV==NULL) 
   {  VV=malloc(sizeof(REAL)*nModelVars);
      for(i=0;i<nModelVars;i++) if(strcmp(varNames[i],"Q")==0) iQ=i;
   }
   for(i=0;i<nModelVars;i++) VV[i]=V[i];
   cErr=0;
   return 0;
}
