#include <stdio.h>
#include <string.h>
#include <stdlib.h>

#include "TV_Algorithms.h"
#include "Structures.h"
#include "Define.h"

#define     OUI                             1
#define     RADIAN_A_DEGRE                  57.295779513082320876798155
#define	   METRES_PAR_SEC_EN_NOEUDS        1.94254
#define     NOEUDS_A_METRES_PAR_SEC         1./METRES_PAR_SEC_EN_NOEUDS
#define     T_ZERO_C	                       273.15

/*!
==================================================================================
\brief Trouver les ingredients necessaire aux orages violent, detecter la 
       structure particuliere aux orages violents, et evaluer le type du 
       temps violent prevu
       
\date 22 janvier 2010
\version 0.1
\author Anna-Belle Filion

param[in]  fichiers_prev   : Nom des fichiers de prevision 
param[out] Grele           : Vecteur de la grosseur des grelons prevu

-----------------------------------------------------------------------------------
 */
int potentielle_temps_violent_( int DEBUG, int nbr_niveau, int ni, int nj, CFG *cfg, 
										  int nbr_var_3D, float Champs3D[nbr_var_3D][nbr_niveau][ni*nj], 
										  int nbr_var_TV, float VarTV[nbr_var_TV][ni*nj] )
{
   if ( DEBUG == OUI )
      fprintf(stderr,"FCN ==> PotentielleTempsViolent\n");  

   int *IndiceT0 = (int*)NULL;
   int *Indice_MaxR_aboveT0 = (int*)NULL;
   
   int k = 0, i, j, ier = 0;
   int nbr_pts = 0;

	float ***champs3D = (float***)NULL;
	float **Variable_TV = (float**)NULL;

	//--------------------------------
   // Allocation de la memoire
   //--------------------------------
   nbr_pts = ni*nj;
   Variable_TV = (float**)malloc(nbr_variableTV*sizeof(float*));
	champs3D = (float***)malloc(nbr_var_3D*sizeof(float**));
	
	//--------------------------------
	// Read the input arrays
	//--------------------------------
	for (i=0;i<nbr_var_3D;i++)
	{
		champs3D[i] = (float**)malloc(nbr_niveau*sizeof(float*));
		
		for (k=0;k<nbr_niveau;k++)
		{
			champs3D[i][k] = (float*)malloc(ni*nj*sizeof(float));
			
			for (j=0;j<(ni*nj);j++)
			{
				champs3D[i][k][j] =  Champs3D[i][k][j];
			}
		}
	}

   //-------------------------------------------------------------
   // Calculer le cisaillement entre 0 et 3km, 0 et 6km
   //------------------------------------------------------------- // Verifier //vent en kts
   if ( DEBUG == OUI )
      fprintf(stderr,"          FCN ==> CalculerCisaillement_3kmET6km\n");
   ier = CalculerCisaillement_3kmET6km(champs3D[GZ], champs3D[WD], champs3D[UV], ni, nj, nbr_niveau, &Variable_TV[Cisaillement_3km], &Variable_TV[Cisaillement_6km]);
   if ( ier != SUCCESS )
   {
      fprintf(stderr,"ERREUR dans la fonction ==> TrouverCisaillement_3kmET6km\n");
      return ERROR;
   }
   
   //----------------------------------
   // Liberer la memoire
   //---------------------------------- 
   for ( k = 0; k < nbr_niveau; k++)
   {
      free(champs3D[UV][k]); 
      free(champs3D[WD][k]);
   }  
   free(champs3D[UV]);
   free(champs3D[WD][k]);
   
   //-------------------------------------
   // Calculer le vertical ice flux
   //-------------------------------------
   if ( DEBUG == OUI )    
   fprintf(stderr,"          FCN ==> CalculerVerticalIceFlux\n");
   ier = CalculerVerticalIceFlux(champs3D[TT], champs3D[WW], champs3D[QJT1], champs3D[QIT1], champs3D[QHT1], ni, nj, nbr_niveau, &Variable_TV[Vertical_Ice_Flux]);
   if ( ier != SUCCESS )
   {
      fprintf(stderr,"ERREUR dans la fonction ==>VerticalIceFlux\n");
      return ERROR;
   }
   
   //----------------------------------
   // Liberer la memoire
   //---------------------------------- 
   for ( k = 0; k < nbr_niveau; k++)
   {
      free(champs3D[QIT1][k]); 
      free(champs3D[QJT1][k]);
      free(champs3D[QHT1][k]);
   }  
   free(champs3D[QIT1]); 
   free(champs3D[QJT1]);
   free(champs3D[QHT1]);
   
   //-----------------------------------------------------
   // Trouver le niveau de congelation pour chaque point
   //----------------------------------------------------- //verifier
   if ( DEBUG == OUI ) 
      fprintf(stderr,"          FCN ==> TrouverNiveauCongelation\n");
   ier = TrouverNiveauCongelation(champs3D[TT], ni, nj, nbr_niveau, &IndiceT0);
   if ( ier != SUCCESS)
   {
      fprintf(stderr,"ERREUR dans la fonction ==> TrouverNiveauCongelation\n");
      return ERROR;
   }
   
   //----------------------------------
   // Liberer la memoire
   //---------------------------------- 
   for ( k = 0; k < nbr_niveau; k++)
   {
      free(champs3D[TT][k]); 
   }  
   free(champs3D[TT]);
   
   //----------------------------------------
   // Trouver la quantite d'eau surfondue 
   //---------------------------------------- //Verifier
   if ( SLW != VALEUR_BIDON )
   {
      if ( DEBUG == OUI )
         fprintf(stderr,"          FCN ==> CalculerQteEauSurfondue\n");
      ier = CalculerQteEauSurfondue( champs3D[SLW], IndiceT0, ni, nj, nbr_niveau, &Variable_TV[Supercooled_Liquid_Water]);
      if ( ier != SUCCESS)
      {
         fprintf(stderr,"ERREUR dans la fonction ==> CalculerQteEauSurfondue\n");
         return ERROR;
      }
      
      //----------------------------------
      // Liberer la memoire
      //---------------------------------- 
      for ( k = 0; k < nbr_niveau; k++)
      {
			free(champs3D[SLW][k]); 
      }  
      free(champs3D[SLW]);
   }
   else
      Variable_TV[Supercooled_Liquid_Water] = (float*)malloc(nbr_pts*sizeof(float));
   
   //--------------------------------------------
   // Calculer le transport vertical d'humidite
   //--------------------------------------------
   if ( DEBUG == OUI )   
      fprintf(stderr,"          FCN ==> CalculerVerticalMoistureFlux\n");
   ier = CalculerVerticalMoistureFlux(champs3D[WW], champs3D[HU], ni, nj, nbr_niveau, &Variable_TV[Vertical_Moisture_Flux],&Variable_TV[Vertical_Moisture_Flux_NEW]);
   if ( ier != SUCCESS )
   {
      fprintf(stderr,"ERREUR dans la fonction ==> VerticalMoistureFlux\n");
      return ERROR;
   }
   
   //----------------------------------
   // Liberer la memoire
   //---------------------------------- 
   for ( k = 0; k < nbr_niveau; k++)
   {
      free(champs3D[HU][k]); 
   }  
   free(champs3D[HU]);
   
   //--------------------------------
   // Trouver les updraft forte
   //-------------------------------- //Verifier
   if ( DEBUG == OUI )
      fprintf(stderr,"          FCN ==> DetecterUpdraft\n");
   ier =  DetecterUpdraft( champs3D[WW], champs3D[GZ], ni, nj, nbr_niveau, (*cfg).thres_severe_updraft, &Variable_TV[Severe_Updraft]);
   if ( ier != SUCCESS)
   {
      fprintf(stderr,"ERREUR dans la fonction ==> DetecterUpdraft\n");
      return ERROR;
   }

   //----------------------------------
   // Liberer la memoire
   //---------------------------------- 
   for ( k = 0; k < nbr_niveau; k++)
   {
      free(champs3D[WW][k]); 
   }  
   free(champs3D[WW]);
   
   //--------------------------------------------------------------------
   // Trouver la reflectivite maximum en haut du point de congelation
   //-------------------------------------------------------------------- //Verifier
   if ( DEBUG == OUI ) 
      fprintf(stderr,"          FCN ==> TrouverMaxRAboveT0\n");
   ier = TrouverMaxRAboveT0(champs3D[ZET], champs3D[GZ], IndiceT0, ni, nj, nbr_niveau, 0., &Variable_TV[MaxR_aboveT0], &Indice_MaxR_aboveT0, &Variable_TV[Hauteur_MaxR_aboveT0]);
   if ( ier != SUCCESS)
   {
      fprintf(stderr,"ERREUR dans la fonction ==> TrouverMaxRAboveT0\n");
      return ERROR;
   }
   
   //-----------------------------------------
   // Trouver la grele prevu par le modele
   //-----------------------------------------
   if ( DMH != VALEUR_BIDON )
   {
      if ( DEBUG == OUI )   
         fprintf(stderr,"          FCN ==> TrouverGrele\n");
      ier = TrouverGrele(champs3D[DMH], Variable_TV[MaxR_aboveT0], ni, nj, nbr_niveau, &Variable_TV[Grosseur_Grele_Max], &Variable_TV[Grele_sfc]);
      if ( ier != SUCCESS )
      {
         fprintf(stderr,"ERREUR dans la fonction ==> TrouverGrele\n");
         return ERROR;
      }
      
      //----------------------------------
      // Liberer la memoire
      //---------------------------------- 
      for ( k = 0; k < nbr_niveau; k++)
      {
			free(champs3D[DMH][k]); 
      }  
      free(champs3D[DMH]);
   }
   else
   {
      Variable_TV[Grosseur_Grele_Max] = (float*)malloc(nbr_pts*sizeof(float));
      Variable_TV[Grele_sfc] = (float*)malloc(nbr_pts*sizeof(float));
   }
     
   //---------------------------
   // Liberer la memoire
   //---------------------------
   free(IndiceT0);
   free(Indice_MaxR_aboveT0);
   IndiceT0 = NULL;
   Indice_MaxR_aboveT0 = NULL;
     
   //---------------------------------------
   // Trouver les couplets de trourbillon
   //---------------------------------------
   if ( DEBUG == OUI )
      fprintf(stderr,"          FCN ==> DetecterCoupletTourbillion\n");
   ier =  DetecterCoupletTourbillion( champs3D[QR], Variable_TV[MaxR_aboveT0], champs3D[GZ], ni, nj, nbr_niveau, (*cfg).thres_couplet_tourbillon, \
   &Variable_TV[Couplet_Tourbillon], &Variable_TV[Couplet_Positif], &Variable_TV[Couplet_Negatif]);
   if ( ier != SUCCESS)
   {
      fprintf(stderr,"ERREUR dans la fonction ==> DetecterCoupletTourbillion\n");
      return ERROR;
   } 
   
   //------------------------------------------------
   // Verifier si il y a des tornades
   //------------------------------------------------
   if ( (*cfg).resolution <= 1.0 )
   {
      if ( DEBUG == OUI )   
         fprintf(stderr,"          FCN ==> DetecterTornade\n");
      ier = DetecterTornade(champs3D[QR], champs3D[GZ], ni, nj, nbr_niveau, (*cfg).thres_tornade, &Variable_TV[Tornade]);
      if ( ier != SUCCESS)
      {
         fprintf(stderr,"ERREUR dans la fonction ==> DetecterTornade");
         return ERROR;
      }
   }
   else
      Variable_TV[Tornade] = (float*)malloc(nbr_pts*sizeof(float));
     
   //-----------------------------------
   // Trouver les zones de surplomb
   //-----------------------------------
   if ( DEBUG == OUI )    
      fprintf(stderr,"          FCN ==> repererSurplomb\n");
   ier = repererSurplomb(champs3D[ZET], champs3D[GZ], Variable_TV[MaxR_aboveT0], ni, nj, nbr_niveau, (*cfg).resolution, &Variable_TV[zoneSurplomb], &Variable_TV[VEF]);
   if ( ier != SUCCESS)
   {
      fprintf(stderr,"ERREUR dans la fonction ==> repererSurplomb\n");
      return ERROR;
   }
   
   //--------------------------------------
   // Verifier si il y a des mesocyclones
   //--------------------------------------
   if ( DEBUG == OUI )   
      fprintf(stderr,"          FCN ==> DetecterMesocyclone\n");
   ier = DetecterMesocyclone(champs3D[QR], champs3D[GZ], Variable_TV[zoneSurplomb], Variable_TV[MaxR_aboveT0], ni, nj, nbr_niveau, (*cfg).thres_meso, \
	&Variable_TV[Mesocyclone], &Variable_TV[Meso_base]);
   if ( ier != SUCCESS )
   {
      fprintf(stderr,"ERREUR dans la fonction ==>DetecterMesocyclone\n");
      return ERROR;
   }
   
   //----------------------------------
   // Liberer la memoire
   //---------------------------------- 
   for ( k = 0; k < nbr_niveau; k++)
   {
      free(champs3D[QR][k]); 
   }  
   free(champs3D[QR]);
   
   //----------------------------------
   // Liberer la memoire
   //---------------------------------- 
   for ( k = 0; k < nbr_niveau; k++)
   {
      free(champs3D[GZ][k]); 
      free(champs3D[ZET][k]); 
   }  
   free(champs3D[GZ]); 
   free(champs3D[ZET]); 
    
	//--------------------------
   // Retourner les resultats
   //-------------------------- 
	for (i=0;i<19;i++)
	{
		for (j=0;j<(ni*nj);j++)
		{
			VarTV[i][j] = Variable_TV[i][j];
		}
		free(Variable_TV[i]);
	}
	free(Variable_TV);
   
   if ( DEBUG == OUI )	
      fprintf(stderr,"FCN ==> PotentielleTempsViolent completed with success\n\n");
   return SUCCESS;
}
