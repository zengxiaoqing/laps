#include <stdio.h>
#include <string.h>
#include <math.h>
#include <stdlib.h>

#include "CtesMathMeteo.h"
#include "Define.h"
#include "Structures.h"

#define OUI         1

/*!
==================================================================================
\brief Determiner la probabilite de temps violents a chaque point
\date 22 janvier 2010
\version 0.1
\author Anna-Belle Filion

param[in]  ni     : Dimension en x
param[in]  nj     : Dimension en y 
param[in]  region : Region de prevision (AB => Alberta, ON => Ontario, QC => Quebec 
param[in]  Indice : Tableau contenant les indice utilent a la prevision d'orage 
param[in]  VarTV  : Tableau contenant les ingredients/signatures necessaires et 
                    particulier aux orages violents
param[out] ProbTV : Vecteur contenant l'indice d'intensite d'orage violent a
                    chaque point de grille

-----------------------------------------------------------------------------------
 */
int calculprobtempsviolent_( int DEBUG, int ni, int nj, CFG *cfg, int nbr_var_2D,
									 int nbr_var_TV, float Champs2D[nbr_var_2D][ni*nj], 
									 float VarTV[nbr_var_TV][ni*nj], float Prob_TV[ni*nj])
{
   if ( DEBUG == OUI )
      fprintf(stderr,"FCN ==> CalculProbTempsViolent\n");

   int i = 0, j = 0, k = 0, m = 0, p = 0;
   int Surplomb = 0, Updraft = 0, Meso = 0, WER = 0;
   
   float sweat = 0, vmf = 0, lwc = 0, lwc700 = 0;
   float vif = 0, prob = 0, shear3 = 0, shear6 = 0;  
   
   int *check = (int*)NULL;
   int *liste = (int*)NULL;  
   
   float *Prob = (float*)NULL;
   float **Champs = (float**)NULL;
	
	char *region = (char*)NULL;
	
	region=cfg->region;

   //--------------------------------
   // Allocation de la memoire
   //--------------------------------
   check = (int*)malloc(ni*nj*sizeof(int));
   Prob  = (float*)malloc(ni*nj*sizeof(float));  
   Champs = (float**)malloc(nbr_var_2D*sizeof(float*));
	
   //-----------------------------
   // Initialiser les variables
   //-----------------------------
   for ( i = 0; i < (ni*nj); i++)
   {
      check[i] = 0;
      Prob[i]  = 0;
   }

   for (i=0;i<nbr_var_2D;i++)
	{
		Champs[i] = (float*)malloc(ni*nj*sizeof(float));
		
		for (j=0;j<(ni*nj);j++)
		{
			Champs[i][j] =  Champs2D[i][j];
		}
	}
   
   //----------------------------------------------
   // Identifie les elements/ingredients presents
   // et necessaires aux temps violents
   //----------------------------------------------
   for ( i = 0; i < (ni*nj); i++)
   {
      if ( (VarTV[Hauteur_MaxR_aboveT0][i] > 40. && check[i] == 0) || ( (VarTV[zoneSurplomb][i] == 1 || VarTV[Severe_Updraft][i] == 1 || VarTV[Mesocyclone][i] == 1) && check[i] == 0) )
      {
         if ( VarTV[zoneSurplomb][i] == 1 || VarTV[Severe_Updraft][i] == 1 || VarTV[Mesocyclone][i] == 1 )
         {
            m = 1;
            p = 0;
            
            liste = (int*)malloc(m*sizeof(int));
            liste[p] = i;
            check[i] = 1;
               
            while ( p < m )
            {
               if ( (liste[p]+1) < (ni*nj) && ( VarTV[zoneSurplomb][liste[p]+1] == 1 || VarTV[Severe_Updraft][liste[p]+1] == 1 || VarTV[Mesocyclone][liste[p]+1] == 1 )  && check[liste[p]+1] == 0)
               {
                  liste = (int*)realloc(liste,(m+1)*sizeof(int));
                  liste[m] = liste[p]+1;
                  check[liste[p]+1] = 1;
                  m++;
               }
               if ( (liste[p]-1) > 0 && ( VarTV[zoneSurplomb][liste[p]-1] == 1 || VarTV[Severe_Updraft][liste[p]-1] == 1 || VarTV[Mesocyclone][liste[p]-1] == 1 ) && check[liste[p]-1] == 0)
               {
                  liste = (int*)realloc(liste,(m+1)*sizeof(int));
                  liste[m] = liste[p]-1;
                  check[liste[p]-1] = 1;
                  m++;
               }
               if ( (liste[p]+ni) < (ni*nj) && ( VarTV[zoneSurplomb][liste[p]+ni] == 1 || VarTV[Severe_Updraft][liste[p]+ni] == 1 || VarTV[Mesocyclone][liste[p]+ni] == 1 ) && check[liste[p]+ni] == 0)
               {
                  liste = (int*)realloc(liste,(m+1)*sizeof(int));
                  liste[m] = liste[p]+ni;
                  check[liste[p]+ni] = 1;
                  m++;
               }
               if ( (liste[p]+ni+1) < (ni*nj) && ( VarTV[zoneSurplomb][liste[p]+ni+1] == 1 || VarTV[Severe_Updraft][liste[p]+ni+1] == 1 || VarTV[Mesocyclone][liste[p]+ni+1] == 1 ) && check[liste[p]+ni+1] == 0)
               {
                  liste = (int*)realloc(liste,(m+1)*sizeof(int));
                  liste[m] = liste[p]+ni+1;
                  check[liste[p]+ni+1] = 1;
                  m++;
               }
               if ( (liste[p]+ni-1) < (ni*nj) && ( VarTV[zoneSurplomb][liste[p]+ni-1] == 1 || VarTV[Severe_Updraft][liste[p]+ni-1] == 1 || VarTV[Mesocyclone][liste[p]+ni-1] == 1 ) && check[liste[p]+ni-1] == 0)
               {
                  liste = (int*)realloc(liste,(m+1)*sizeof(int));
                  liste[m] = liste[p]+ni-1;
                  check[liste[p]+ni-1] = 1;
                  m++;
               }
               if ( (liste[p]-ni) > 0 && ( VarTV[zoneSurplomb][liste[p]-ni] == 1 || VarTV[Severe_Updraft][liste[p]-ni] == 1 || VarTV[Mesocyclone][liste[p]-ni] == 1 ) && check[liste[p]-ni] == 0)
               {
                  liste = (int*)realloc(liste,(m+1)*sizeof(int));
                  liste[m] = liste[p]-ni;
                  check[liste[p]-ni] = 1;
                  m++;
               }
               if ( (liste[p]-ni+1) > 0 && ( VarTV[zoneSurplomb][liste[p]-ni+1] == 1 || VarTV[Severe_Updraft][liste[p]-ni+1] == 1 || VarTV[Mesocyclone][liste[p]-ni+1] == 1 ) && check[liste[p]-ni+1] == 0)
               {
                  liste = (int*)realloc(liste,(m+1)*sizeof(int));
                  liste[m] = liste[p]-ni+1;
                  check[liste[p]-ni+1] = 1;
                  m++;
               }
               if ( (liste[p]-ni-1) > 0 && ( VarTV[zoneSurplomb][liste[p]-ni-1] == 1 || VarTV[Severe_Updraft][liste[p]-ni-1] == 1 || VarTV[Mesocyclone][liste[p]-ni-1] == 1 ) && check[liste[p]-ni-1] == 0)
               {
                  liste = (int*)realloc(liste,(m+1)*sizeof(int));
                  liste[m] = liste[p]-ni-1;
                  check[liste[p]-ni-1] = 1;
                  m++;
               }
               p++;
            }
             
            Surplomb = 0;
            Updraft  = 0;
            Meso     = 0;
            sweat    = 0;
            vmf      = 0;
            lwc      = 0;
            lwc700   = 0;
            WER      = 0;
            shear3   = 0;
            shear6   = 0;
            vif      = 0;
               
            for ( k = 0; k < m; k++ )
            {
               if ( VarTV[zoneSurplomb][liste[k]] == 1 )
                  Surplomb = 1;
               if ( VarTV[Severe_Updraft][liste[k]] == 1 )
                  Updraft = 1;
               if ( VarTV[Mesocyclone][liste[k]] == 1 )
                  Meso = 1;
               if ( VarTV[VEF][liste[k]] == 1 )
                  WER = 1;
               if ( Champs[SWEAT][liste[k]] > sweat )
                  sweat = Champs[SWEAT][liste[k]];
               if ( VarTV[Vertical_Moisture_Flux][liste[k]] > vmf )
                  vmf = VarTV[Vertical_Moisture_Flux][liste[k]];
               if ( (Champs[LWC][liste[k]]*1000) > lwc )
                  lwc = Champs[LWC][liste[k]]  * 1000.;
               if ( (Champs[LWC_700][liste[k]]*1000)> lwc700 )
                  lwc700 = Champs[LWC_700][liste[k]]  * 1000.;
               if ( VarTV[Cisaillement_3km][liste[k]] > 60 )
                  shear3 = 1;
               if ( VarTV[Cisaillement_6km][liste[k]] > 90 )
                  shear6 = 1;
               if ( QIT1 != VALEUR_BIDON && QJT1 != VALEUR_BIDON && QHT1 != VALEUR_BIDON )
               {
                  if ( VarTV[Vertical_Ice_Flux][liste[k]] > vif )
                     vif = VarTV[Vertical_Ice_Flux][liste[k]];
               }
            }
            
            prob = 0;
            if ( strcmp(region,"AB") == 0 )
            {
	       if ( Surplomb == 1 && Updraft == 1 && Meso == 1 && WER == 1 )
		  prob = 100;
               else if ( Surplomb == 1 && Updraft == 1 && Meso == 1 )
                  prob = 90;
	       else if ( Surplomb == 1 && Updraft == 1 && Meso == 0 && WER == 1 && lwc > 25000 && lwc700 > 10000 && vmf > 4000 )
                  prob = 80;		  
	       else if ( Surplomb == 1 && Updraft == 0 && Meso == 1 && WER == 1 && lwc > 20000 && lwc700 > 10000 && vmf > 3500 )
                  prob = 70;
	       else if ( Surplomb == 0 && Updraft == 1 && Meso == 0 && WER == 0 && lwc > 30000 && lwc700 > 12000 && vmf > 5000 && vif > 100 && (shear3 == 1 || shear6 == 1) )
                  prob = 70;
	       else if ( Surplomb == 1 && Updraft == 1 && Meso == 0 && WER == 1 && lwc > 20000 && lwc700 > 10000 && vmf > 3500 )
                  prob = 60;
               else if ( Surplomb == 1 && Updraft == 1 && Meso == 0 && WER == 0 && lwc > 25000 && lwc700 > 10000 && vmf > 4000 )
                  prob = 60;
	       else if ( Surplomb == 1 && Updraft == 1 && Meso == 0 && WER == 0 && lwc > 20000 && lwc700 > 10000 && vmf > 3500 )
                  prob = 50;	
	       else if ( Surplomb == 0 && Updraft == 1 && Meso == 0 && WER == 0 && lwc > 25000 && lwc700 > 11000 && vmf > 4000 && vif > 150 )
                  prob = 50;	
	       else if ( Surplomb == 1 && Updraft == 0 && Meso == 1 && WER == 0 && lwc > 20000 && lwc700 > 10000 && vmf > 3500 )
                  prob = 40; 
	       else if ( Surplomb == 1 && Updraft == 0 && Meso == 0 && WER == 0 && lwc > 25000 && lwc700 > 10000 && vmf > 4000 )
                  prob = 40;
	       else if ( Surplomb == 0 && Updraft == 1 && Meso == 0 && WER == 0 && lwc > 25000 && lwc700 > 10000 && vmf > 4000 && vif > 100 && (shear3 == 1 || shear6 ==1) )
		  prob = 40;
	       else if ( Surplomb == 1 && Updraft == 0 && Meso == 0 && WER == 1 && lwc > 20000 && lwc700 > 10000 && vmf > 3500 )
                  prob = 30; 
	       else if ( Surplomb == 1 && Updraft == 0 && Meso == 0 && WER == 0 && lwc > 20000 && lwc700 > 10000 && vmf > 3500 )
                  prob = 20;
               else if ( Surplomb == 1 && Updraft == 0 && Meso == 0 && WER == 0 && lwc > 25000 && lwc700 > 10000 && vmf < 4000 )
                  prob = 20;
               else if ( Surplomb == 1 && Updraft == 0 && Meso == 0 && WER == 0 && lwc > 20000 && lwc700 > 10000 && vmf < 4000 )
                  prob = 10;
               else if ( Surplomb == 0 && Updraft == 1 && Meso == 0 && WER == 0 && lwc > 25000 && lwc700 > 10000 && vmf > 4000 && shear3 == 1 )
                  prob = 10;
            }
            else if ( strcmp(region,"ON") == 0 )
            {
	       if ( Surplomb ==1 && Updraft == 1 && Meso == 1 && WER == 1 )
                  prob = 100;
	       else if ( Surplomb ==1 && Updraft == 1 && Meso == 1 )
		  prob = 90;
               else if ( Surplomb == 1 && Updraft == 1 && Meso == 0 && WER == 1 && lwc > 40000 && lwc700 > 10000 && vmf > 6000 && sweat > 200 )
                  prob = 80;
	       else if ( Surplomb == 1 && Updraft == 1 && Meso == 0 && WER == 0 && lwc > 35000 && lwc700 > 10000 && vmf > 6000 && sweat > 200 )
                  prob = 70;
	       else if ( Surplomb == 1 && Updraft == 1 && Meso == 0 && WER == 0 && lwc > 30000 && lwc700 > 9000 && vmf > 5000 && sweat > 300 )
                  prob = 70;
	       else if ( Surplomb == 1 && Updraft == 0 && Meso == 0 && WER == 0 && lwc > 40000 && lwc700 > 10000 && vmf > 6000 && sweat > 300 )
                  prob = 60;
	       else if ( Surplomb == 1 && Updraft == 0 && Meso == 0 && WER == 1 && lwc > 35000 && lwc700 > 10000 && vmf > 5000 && sweat > 200 )
                  prob = 60;
	       else if ( Surplomb == 1 && Updraft == 1 && Meso == 0 && WER == 1 && lwc > 35000 && lwc700 > 9000 && vmf > 5000 && sweat > 200 )
                  prob = 50;
	       else if ( Surplomb == 0 && Updraft == 1 && Meso == 0 && WER == 0 && lwc > 40000 && lwc700 > 12000 && vmf > 6000 && sweat > 200 )
                  prob = 50;
	       else if ( Surplomb == 1 && Updraft == 0 && Meso == 1 && WER == 1 && lwc > 35000 && lwc700 > 9000 && vmf > 5000 && sweat > 200 )
                  prob = 40;
               else if ( Surplomb == 1 && Updraft == 1 && Meso == 0 && WER == 0 && lwc > 30000 && lwc700 > 9000 && vmf > 5000 && sweat > 200 )
                  prob = 40;
	       else if ( Surplomb == 1 && Updraft == 1 && Meso == 0 && WER == 0 && lwc > 30000 && lwc700 > 9000 && vmf < 5000 && sweat > 200 )
                  prob = 30;
               else if ( Surplomb == 1 && Updraft == 1 && Meso == 0 && WER == 0 && lwc > 30000 && lwc700 < 9000 && vmf > 5000 && sweat > 200 )
                  prob = 30;
               else if ( Surplomb == 1 && Updraft == 1 && Meso == 0 && WER == 0 && lwc < 30000 && lwc700 > 9000 && vmf > 5000 && sweat > 200 )
                  prob = 30;
               else if ( Surplomb == 1 && Updraft == 0 && Meso == 1 && WER == 0 && lwc > 40000 && lwc700 > 10000 && vmf > 6000 && sweat > 200 )
                  prob = 20;
               else if ( Surplomb == 1 && Updraft == 1 && Meso == 0 && WER == 0 && lwc > 30000 && lwc700 > 9000 && vmf > 5000 && sweat < 200 )
                  prob = 20;
               else if ( Surplomb == 1 && Updraft == 0 && Meso == 1 && WER == 0 && (lwc < 40000 || lwc700 < 10000 || vmf < 6000) && sweat > 200 )
                  prob = 10;
               if ( sweat < 200. )
                  prob = 0;
            }
            else if ( strcmp(region,"QC") == 0)
            {
   	       if ( Surplomb ==1 && Updraft == 1 && Meso == 1 && WER == 1 )
		  prob = 100;
	       else if ( Surplomb ==1 && Updraft == 1 && Meso == 1 )
                  prob = 90;
	       else if ( Surplomb == 1 && Updraft == 1 && Meso == 0 && WER == 1 && lwc > 40000 && lwc700 > 10000 && vmf > 6000 && sweat > 200 )
                  prob = 80;
	       else if ( Surplomb == 1 && Updraft == 1 && Meso == 0 && WER == 1 && lwc > 35000 && lwc700 > 9000 && vmf > 5000 && sweat > 200 )
                  prob = 70;
	       else if ( Surplomb == 1 && Updraft == 0 && Meso == 1 && WER == 1 && lwc > 35000 && lwc700 > 9000 && vmf > 5000 && sweat > 200 )
                  prob = 70;
	       else if ( Surplomb == 0 && Updraft == 1 && Meso == 0 && WER == 0 && lwc > 40000 && lwc700 > 10000 && vmf > 5000 && sweat > 200 && shear3 == 1 && shear6 == 1 )
                  prob = 60;
	       else if ( Surplomb == 1 && Updraft == 1 && Meso == 0 && WER == 1 && lwc > 30000 && lwc700 > 8000 && vmf > 5000 && sweat > 200 )
                  prob = 50;
	       else if ( Surplomb == 1 && Updraft == 1 && Meso == 0 && WER == 0 && lwc > 35000 && lwc700 > 10000 && vmf > 5000 && sweat > 200 )
                  prob = 50;	
	       else if ( Surplomb == 0 && Updraft == 1 && Meso == 0 && WER == 0 && lwc > 30000 && lwc700 > 9000 && vmf > 5000 && sweat > 200 && shear3 == 1 && shear6 == 1 )
                  prob = 40;
	       else if ( Surplomb == 1 && Updraft == 0 && Meso == 0 && WER == 1 && lwc > 35000 && lwc700 > 10000 && vmf > 5000 && sweat > 200 )
                  prob = 30;
	       else if ( Surplomb == 1 && Updraft == 0 && Meso == 0 && WER == 0 && lwc > 40000 && lwc700 > 10000 && vmf > 5000 && sweat > 300 )
                  prob = 30;	
	       else if ( Surplomb == 1 && Updraft == 0 && Meso == 1 && WER == 0 && lwc > 35000 && lwc700 > 9000 && vmf > 5000 && sweat > 200 )
                  prob = 20;	
	       else if ( Surplomb == 0 && Updraft == 1 && Meso == 0 && WER == 0 && lwc > 25000 && lwc700 > 8000 && vmf > 5000 && sweat > 200 && shear3 == 1 && shear6 == 1 )
                  prob = 10;
	       if ( sweat < 200. )
               prob = 0;
            } 
            else if ( strcmp(region,"EAST") == 0 )
            {
               if ( Surplomb ==1 && Updraft == 1 && Meso == 1 && WER == 1 )
                  prob = 100;
               else if ( Surplomb ==1 && Updraft == 1 && Meso == 1 )
                  prob = 90;
               /*else if ( Surplomb == 1 && Updraft == 1 && Meso == 0 && WER == 1 && lwc > 50000 && lwc700 > 20000 && vmf > 12000 && sweat > 200 )
                  prob = 80;
               else if ( Surplomb == 1 && Updraft == 1 && Meso == 0 && WER == 1 && lwc > 45000 && lwc700 > 20000 && vmf > 12000 && sweat > 200 )
                  prob = 70;
               else if ( Surplomb == 1 && Updraft == 0 && Meso == 1 && WER == 1 && lwc > 45000 && lwc700 > 20000 && vmf > 10000 && sweat > 200 )
                  prob = 70;*/
               else if ( Surplomb == 1 && Updraft == 1 && Meso == 0 && WER == 0 && lwc > 50000 && lwc700 > 20000 && vmf > 10000 && sweat > 200 )
                  prob = 70;
               else if ( Surplomb == 1 && Updraft == 1 && Meso == 0 && WER == 0 && lwc > 45000 && lwc700 > 19000 && vmf > 10000 && sweat > 200 )
                  prob = 60;
               else if ( Surplomb == 0 && Updraft == 1 && Meso == 0 && WER == 0 && lwc > 50000 && lwc700 > 20000 && vmf > 12000 && sweat > 200 && shear3 == 1 && shear6 == 1 )
                  prob = 60;
               /*else if ( Surplomb == 1 && Updraft == 1 && Meso == 0 && WER == 1 && lwc > 40000 && lwc700 > 18000 && vmf > 10000 && sweat > 200 )
                  prob = 60;*/
               /*else if ( Surplomb == 1 && Updraft == 0 && Meso == 0 && WER == 1 && lwc > 40000 && lwc700 > 18000 && vmf > 10000 && sweat > 200 )
                  prob = 50;*/
               else if ( Surplomb == 0 && Updraft == 1 && Meso == 0 && WER == 0 && lwc > 50000 && lwc700 > 20000 && vmf > 10000 && sweat > 200 && shear3 == 1 && shear6 == 1 )
                  prob = 50;
               else if ( Surplomb == 1 && Updraft == 1 && Meso == 0 && WER == 0 && lwc > 40000 && lwc700 > 18000 && vmf > 10000 && sweat > 200 )
                  prob = 50;             
               else if ( Surplomb == 1 && Updraft == 0 && Meso == 0 && WER == 0 && lwc > 50000 && lwc700 > 20000 && vmf > 10000 && sweat > 300 )
                  prob = 40;
               /*else if ( Surplomb == 1 && Updraft == 0 && Meso == 1 && WER == 1 && lwc > 40000 && lwc700 > 18000 && vmf > 10000 && sweat > 200 )
                  prob = 40;*/
               else if ( Surplomb == 0 && Updraft == 1 && Meso == 0 && WER == 0 && lwc > 45000 && lwc700 > 19000 && vmf > 10000 && sweat > 200 && shear3 == 1 && shear6 == 1 )
                  prob = 40;
               else if ( Surplomb == 1 && Updraft == 1 && Meso == 0 && WER == 0 && lwc > 40000 && lwc700 > 18000 && vmf < 10000 && sweat > 200 )
                  prob = 30;
               else if ( Surplomb == 1 && Updraft == 1 && Meso == 0 && WER == 0 && lwc > 40000 && lwc700 < 18000 && vmf > 10000 && sweat > 200 )
                  prob = 30;
               else if ( Surplomb == 1 && Updraft == 1 && Meso == 0 && WER == 0 && lwc < 40000 && lwc700 > 18000 && vmf > 10000 && sweat > 200 )
                  prob = 30;
               else if ( Surplomb == 1 && Updraft == 0 && Meso == 0 && WER == 0 && lwc > 45000 && lwc700 > 19000 && vmf > 1000 && sweat > 200 )
                  prob = 20;
               else if ( Surplomb == 1 && Updraft == 1 && Meso == 0 && WER == 0 && lwc > 40000 && lwc700 > 18000 && vmf > 10000 && sweat < 200 )
                  prob = 20;
               else if ( Surplomb == 1 && Updraft == 0 && Meso == 1 && WER == 0 && lwc > 45000 && lwc700 > 18000 && vmf > 7000 && sweat > 200 )
                  prob = 10; 
               else if ( Surplomb == 0 && Updraft == 1 && Meso == 0 && WER == 0 && lwc > 40000 && lwc700 > 18000 && vmf > 10000 && sweat > 200 && shear3 == 1 && shear6 == 1 )
                  prob = 10;
               if ( sweat < 200. )
               prob = 0;
            }  

            for ( k = 0; k < m; k++ )
            {
               Prob[liste[k]] = prob;
            }
            free(liste);
            liste = NULL;
         }
      }
   }
   
   //----------------------------
   // Retourner les variables
   //----------------------------   
	for (j=0;j<(ni*nj);j++)
	{
		Prob_TV[j] = Prob[j];
	}
	free(Prob);
   
   if ( DEBUG == OUI )
      fprintf(stderr,"FCN ==> CalculProbTempsViolent completed with success\n\n");
   return SUCCESS;
}

