#include <stdio.h>
#include <string.h>
#include <math.h>
#include <stdlib.h>

#include "CtesMathMeteo.h"

/*!
=============================================================================
\brief Trouver si le modele a prevu de la grele, si oui trouver la grosseur
 maximum de la grele dans la colonne et la grosseur de la grele a la surface 


\param[in]  DMH         : Tableau de diametre moyen de la grele (m)
\param[in]  ni          : Dimension en x
\param[in]  nj          : Dimension en y
\param[in]  nk          : Dimension en z
\param[out] HailSizeSfc : Vecteur de la grosseur moyenne de la grele a la
                             surface (mm)
\param[out] HailSizeMax : Vecteur de grosseur maximum de la grele dans
                             la colonne (mm)
\param[out] Hail        : Vecteur de presence de grele dans la colonne
                              0 => Pas present 1 => Present 


\date 2 sept 2010
\version 0.1
\author Anna-Belle Filion

-----------------------------------------------------------------------------
 */
int TrouverGrele( float **DMH, float *MaxR, int ni, int nj, int nk, 
                  float **HailMax, float **HailSfc)
{
   int i, k;
   
   float max;
   
   float *GreleMax = (float*)NULL;
   float *GreleSfc = (float*)NULL;
   
   //---------------------------
   // Allocation de memoire
   //---------------------------
   GreleMax = (float*)malloc(ni*nj*sizeof(float));
   GreleSfc = (float*)malloc(ni*nj*sizeof(float));
   
   for ( i = 0; i < (ni*nj); i++ )
   {
      if ( MaxR[i] > 40. )
      {
         max = 0;
         for ( k = 0; k < nk; k++ )
         {
            if ( DMH[k][i] > max )
               max = DMH[k][i];
         }
         if ( max != 0 )
         {
            GreleMax[i] = max * 1000.;
         }
         else
         {
            GreleMax[i] = 0;
         }
         GreleSfc[i] = DMH[1][i] * 1000.;
      }
      else
      {
         GreleMax[i] = 0.;
         GreleSfc[i] = 0.;
      }
   }
   
   //--------------------------
   // Retourner les resultats
   //--------------------------
   *HailMax = GreleMax;
   *HailSfc = GreleSfc;
         
   return SUCCESS;
}


/*!
=============================================================================
\brief Verifier si le modele a prevu des rafales severes (selon Brasseur)
       les rafales sont considere severe si la rafale estimee > 70km/h


\param[in]  WGE            : Vecteur de rafales estimees (m/s)
\param[in]  WGX            : Vecteur de rafale maximum   (m/s)
\param[in]  ni             : Dimension en x
\param[in]  nj             : Dimension en y
\param[out] SevereGust     : Vecteur de rafale severe prevu (km/h)
\param[out] SevereGustMax  : Vecteur de rafale severe maximale prevu (km/h)


\date 2 sept 2010
\version 0.1
\author Anna-Belle Filion

-----------------------------------------------------------------------------
 */
int BrasseurSevereGust( float *WGE, float *WGX, int ni, int nj, float **SevereGust,
                        float **SevereGustMax)
{
   int i;
   
   float *Gust    = (float*)NULL;
   float *GustMax = (float*)NULL;
    
   //---------------------------
   // Allocation de memoire
   //---------------------------
   Gust    = (float*)malloc(ni*nj*sizeof(float));
   GustMax = (float*)malloc(ni*nj*sizeof(float));
   
   for ( i = 0; i < (ni*nj); i++ )
   {
      if ( (WGE[i] * METRES_PAR_SEC_A_KMH) >= 70. )
      {
         Gust[i]    = WGE[i] * METRES_PAR_SEC_A_KMH;
         GustMax[i] = WGX[i] * METRES_PAR_SEC_A_KMH;
      }
      else
      {
         Gust[i]    = 0;
         GustMax[i] = 0;
      }
   }
   
   //--------------------------
   // Retourner les resultats
   //--------------------------
   *SevereGust    = Gust;
   *SevereGustMax = GustMax;
   
   return SUCCESS;
}




/*!
=============================================================================
\brief Calculer la quantite d'eau surfondue en haut


\param[in]  SLW            : Tableau d'eau surfondue (kg/m^3)
\param[in]  indiceT0       : Vecteur d'indices du point de congelation
\param[in]  ni             : Dimension en x
\param[in]  nj             : Dimension en y
\param[in]  nk             : Dimension en k
\param[out] EauSurfondue   : Vecteur de quantite d'eau surfondue (g/m^3)


\date 2 sept 2010
\version 0.1
\author Anna-Belle Filion

-----------------------------------------------------------------------------
 */
int CalculerQteEauSurfondue( float **SLW, int *indiceT0, int ni, int nj, int nk,
                             float **EauSurfondue)
{
   int i, k;
   
   float *QteEauSurf = (float*)NULL;
   
   //----------------------------
   // Allocation de memoire
   //----------------------------
   QteEauSurf = (float*)malloc(ni*nj*sizeof(float));

   for ( i = 0; i < (ni*nj); i++)
   {
      QteEauSurf[i] = 0.;
      for ( k = (indiceT0[i]+1); k < nk; k++)
      {
         QteEauSurf[i] += (SLW[k][i] * 1000.);
      }
   }
   
   //--------------------------
   // Retourner les variables
   //--------------------------
   *EauSurfondue = QteEauSurf;
   
   return SUCCESS;
}


/*!
=============================================================================
\brief Detection des mesocyclones


\param[in]  QR               : Tableau du tourbillon (s^(-1))
\param[in]  GZ               : Tableau de la hauteur geopotentielle (Dam)
\param[in]  MaxR             : Vecteur de reflectivite maximum dans une colonne (dBz)
\param[in]  ni               : dimension en x
\param[in]  nj               : dimension en y
\param[in]  nk               : dimension en k
\param[out] Mesocyclone      : Vecteur des mesocyclones detectes ( 1 => Oui / 0 => Non )
\param[out] BaseMesoyclone   : Vecteur de la base des mesocyclones detectes (m)

\date 24 aout 2010
\version 0.1
\author Anna-Belle Filion

-----------------------------------------------------------------------------
 */
int DetecterMesocyclone(float **QR, float **GZ, float *Surplomb, float *MaxR,
                        int ni, int nj, int nk, float seuil, float **Mesocyclone, 
                        float **BaseMesoyclone)
{
   int i, k, m, n, o, p;
   
   int indice_base, top, base, start;
   
   int *size  = (int*)NULL;
   int *size1 = (int*)NULL;
   int *size2 = (int*)NULL;
   int **check  = (int**)NULL;
   int **liste  = (int**)NULL;
   int **liste1 = (int**)NULL;
   int **liste2 = (int**)NULL;
   
   float *mesocyclone = (float*)NULL;
   float *base_mesocyclone = (float*)NULL;
   
   //-----------------------------
   // Allocation de memoire
   //-----------------------------
   check = (int**)malloc(nk*sizeof(int*));
   for ( k = 0; k < nk; k++ )
   {
      check[k] = (int*)malloc(ni*nj*sizeof(int));
   }
   mesocyclone = (float*)malloc(ni*nj*sizeof(float));
   base_mesocyclone = (float*)malloc(ni*nj*sizeof(float));
   
   //-------------------------------
   // Initialiser les variables
   //-------------------------------
   for ( k = 0; k < nk; k++ )
   {
      for ( i = 0; i < (ni*nj); i++)
      {
         check[k][i] = 0;
         if( k == 0 )
         {
            mesocyclone[i] = 0.;
            base_mesocyclone[i] = 0.;
         }
      }
   }

   //----------------------------------
   // Trouver les mesocyclones
   //----------------------------------
   for ( i = 0; i < (ni*nj); i++)
   {
      if ( MaxR[i] > 40. && Surplomb[i] == 1 )
      {
         for ( k = 0; k < nk; k++)
         {
            if ( (GZ[k][i] - GZ[0][i]) < 1200. && (GZ[k][i] - GZ[0][i]) > 200. )
            {
               m = 0;
               n = 0;
               if ( QR[k][i] >= seuil && check[k][i] == 0 )
               {
		  top=k;
		  base=k;
                  indice_base = k;
                  start = k;
                  size = (int*)malloc(1*sizeof(int));
                  liste = (int**)malloc(1*sizeof(int*));
                  liste[n] = (int*)(malloc)(1*sizeof(int));
                  liste[n][m] = i;
                  check[k][i] = 1;
                  m++;
                  
                  if ( (i+1) < (ni*nj) && ((i+1)%ni) != 0 && QR[k][i+1] > seuil && check[k][i+1] == 0 )
                  {
                     liste[n] = (int*)realloc(liste[n],(m+1)*sizeof(int));
                     liste[n][m] = i+1;
                     check[k][i+1] = 1;
                     m++;
                  }
                  if ( (i-1) >= 0 && (i%ni) != 0 && QR[k][i-1] > seuil && check[k][i-1] == 0 )
                  {
                     liste[n] = (int*)realloc(liste[n],(m+1)*sizeof(int));
                     liste[n][m] = i-1;
                     check[k][i-1] = 1;
                     m++;           
                  }
                  if ( (i+ni-1) < (ni*nj) && (i%ni) != 0 && QR[k][i+ni-1] > seuil && check[k][i+ni-1] == 0 )
                  {
                     liste[n] = (int*)realloc(liste[n],(m+1)*sizeof(int));
                     liste[n][m] = i+ni-1;
                     check[k][i+ni-1] = 1;
                     m++;        
                  }
                  if ( (i+ni) < (ni*nj) && QR[k][i+ni] > seuil && check[k][i+ni] == 0 )
                  {
                     liste[n] = (int*)realloc(liste[n],(m+1)*sizeof(int));
                     liste[n][m] = i+ni;
                     check[k][i+ni] = 1;
                     m++;
                  }
                  if ( (i+ni+1) < (ni*nj) && ((i+1)%ni) != 0 && QR[k][i+ni+1] > seuil && check[k][i+ni+1] == 0 )
                  {
                     liste[n] = (int*)realloc(liste[n],(m+1)*sizeof(int));
                     liste[n][m] = i+ni+1;
                     check[k][i+ni+1] = 1;
                     m++;           
                  }
                  if ( (i-ni-1) >= 0 && (i%ni) != 0 && QR[k][i-ni-1] > seuil && check[k][i-ni-1] == 0 )
                  {
                     liste[n] = (int*)realloc(liste[n],(m+1)*sizeof(int));
                     liste[n][m] = i-ni-1;
                     check[k][i-ni-1] = 1;
                     m++;           
                  }
                  if ( (i-ni) >= 0 && QR[k][i-ni] > seuil && check[k][i-ni] == 0 )
                  {
                     liste[n] = (int*)realloc(liste[n],(m+1)*sizeof(int));
                     liste[n][m] = i-ni;
                     check[k][i-ni] = 1;
                     m++;            
                  }
                  if ( (i-ni+1) >= 0 && ((i+1)%ni) != 0 && QR[k][i-ni+1] > seuil && check[k][i-ni+1] == 0 ) 
                  {
                     liste[n] = (int*)realloc(liste[n],(m+1)*sizeof(int));
                     liste[n][m] = i-ni+1;
                     check[k][i-ni+1] = 1;
                     m++;           
                  }

                  for ( o = 1; o < m; o++)
                  {
                     if ( (liste[n][o]+1) < (ni*nj) && ((liste[n][o]+1)%ni) != 0 && QR[k][liste[n][o]+1] > seuil && check[k][liste[n][o]+1] != 1 )
                     {
                        liste[n] = (int*)realloc(liste[n],(m+1)*sizeof(int));
                        liste[n][m] = liste[n][o]+1;
                        check[k][liste[n][o]+1] = 1;
                        m++;
                     }
                     if ( (liste[n][o]-1) >= 0 && (liste[n][o]%ni) != 0 && QR[k][liste[n][o]-1] > seuil && check[k][liste[n][o]-1] != 1 )
                     {
                        liste[n] = (int*)realloc(liste[n],(m+1)*sizeof(int));
                        liste[n][m] = liste[n][o]-1;
                        check[k][liste[n][o]-1] = 1;
                        m++;
                     }
                     if ( (liste[n][o]+ni-1) < (ni*nj) && (liste[n][o]%ni) != 0 && QR[k][liste[n][o]+ni-1] > seuil && check[k][liste[n][o]+ni-1] != 1 )
                     {
                        liste[n] = (int*)realloc(liste[n],(m+1)*sizeof(int));
                        liste[n][m] = liste[n][o]+ni-1;
                        check[k][liste[n][o]+ni-1] = 1;
                        m++;
                     }
                     if ( (liste[n][o]+ni) < (ni*nj) && QR[k][liste[n][o]+ni] > seuil && check[k][liste[n][o]+ni] != 1 )
                     {
                        liste[n] = (int*)realloc(liste[n],(m+1)*sizeof(int));
                        liste[n][m] = liste[n][o]+ni;
                        check[k][liste[n][o]+ni] = 1;
                        m++;
                     }
                     if ( (liste[n][o]+ni+1) < (ni*nj) && ((liste[n][o]+1)%ni) != 0 && QR[k][liste[n][o]+ni+1] > seuil && check[k][liste[n][o]+ni+1] != 1 )
                     {
                        liste[n] = (int*)realloc(liste[n],(m+1)*sizeof(int));
                        liste[n][m] = liste[n][o]+ni+1;
                        check[k][liste[n][o]+ni+1] = 1;
                        m++;
                     }
                     if ( (liste[n][o]-ni-1) >= 0 && (liste[n][o]%ni) != 0 && QR[k][liste[n][o]-ni-1] > seuil && check[k][liste[n][o]-ni-1] != 1 )
                     {  
                        liste[n] = (int*)realloc(liste[n],(m+1)*sizeof(int));
                        liste[n][m] = liste[n][o]-ni-1;
                        check[k][liste[n][o]-ni-1] = 1;
                        m++;
                     }
                     if ( (liste[n][o]-ni) >= 0 && QR[k][liste[n][o]-ni] > seuil && check[k][liste[n][o]-ni] != 1 )
                     {
                        liste[n] = (int*)realloc(liste[n],(m+1)*sizeof(int));
                        liste[n][m] = liste[n][o]-ni;
                        check[k][liste[n][o]-ni] = 1;
                        m++;
                     }
                     if ( (liste[n][o]-ni+1) > 0 && ((liste[n][o]+1)%ni) != 0 && QR[k][liste[n][o]-ni+1] > seuil && check[k][liste[n][o]-ni+1] != 1 ) 
                     {
                        liste[n] = (int*)realloc(liste[n],(m+1)*sizeof(int));
                        liste[n][m] = liste[n][o]-ni+1;
                        check[k][liste[n][o]-ni+1] = 1;
                        m++;
                     }
                  }

                  size[n] = m;
                  while (  (GZ[k][i] - GZ[0][i]) < 1200. && m != 0 && k < (nk-1) )
                  {
                     n++;
                     k++;
                     liste = (int**)realloc(liste,(n+1)*sizeof(int*));
                     liste[n] = (int*)(malloc)(1*sizeof(int));
                              
                     p = 0;
                     for ( o = 0; o < m; o++)
                     {
                        if ( QR[k][liste[n-1][o]] > seuil && check[k][liste[n-1][o]] == 0 )
                        {
                           liste[n] = (int*)realloc(liste[n],(p+1)*sizeof(int));
                           liste[n][p] = liste[n-1][o];
                           check[k][liste[n][p]] = 1;
                           p++;
                        }
                        if ( (liste[n-1][o]+1) < (ni*nj) && ((liste[n-1][o]+1)%ni) != 0 && QR[k][liste[n-1][o]+1] > seuil && check[k][liste[n-1][o]+1] == 0 )
                        {
                           liste[n] = (int*)realloc(liste[n],(p+1)*sizeof(int));
                           liste[n][p] = liste[n-1][o]+1;
                           check[k][liste[n][p]] = 1;
                           p++;
                        }
                        if ( (liste[n-1][o]-1) > 0 && (liste[n-1][o]%ni) != 0 && QR[k][liste[n-1][o]-1] > seuil && check[k][liste[n-1][o]-1] == 0 )
                        {
                           liste[n] = (int*)realloc(liste[n],(p+1)*sizeof(int));
                           liste[n][p] = liste[n-1][o]-1;
                           check[k][liste[n][p]] = 1;
                           p++;
                        }
                        if ( (liste[n-1][o]+ni-1) < (ni*nj) && (liste[n-1][o]%ni) != 0 && QR[k][liste[n-1][o]+ni-1] > seuil && check[k][liste[n-1][o]+ni-1] == 0 )
                        {
                           liste[n] = (int*)realloc(liste[n],(p+1)*sizeof(int));
                           liste[n][p] = liste[n-1][o]+ni-1;
                           check[k][liste[n][p]] = 1;
                           p++;
                        }
                        if ( (liste[n-1][o]+ni) < (ni*nj) && QR[k][liste[n-1][o]+ni] > seuil && check[k][liste[n-1][o]+ni] == 0 )
                        {
                           liste[n] = (int*)realloc(liste[n],(p+1)*sizeof(int));
                           liste[n][p] = liste[n-1][o]+ni;
                           check[k][liste[n][p]] = 1;
                           p++;
                        }
                        if ( (liste[n-1][o]+ni+1) < (ni*nj) && ((liste[n-1][o]+1)%ni) != 0 && QR[k][liste[n-1][o]+ni+1] > seuil && check[k][liste[n-1][o]+ni+1] == 0 )
                        {
                           liste[n] = (int*)realloc(liste[n],(p+1)*sizeof(int));
                           liste[n][p] = liste[n-1][o]+ni+1;
                           check[k][liste[n][p]] = 1;
                           p++;
                        }
                        if ( (liste[n-1][o]-ni-1) > 0 && (liste[n-1][o]%ni) != 0 && QR[k][liste[n-1][o]-ni-1] > seuil && check[k][liste[n-1][o]-ni-1] == 0 )
                        {
                           liste[n] = (int*)realloc(liste[n],(p+1)*sizeof(int));
                           liste[n][p] = liste[n-1][o]-ni-1;
                           check[k][liste[n][p]] = 1;
                           p++;
                        }
                        if ( (liste[n-1][o]-ni) > 0 && QR[k][liste[n-1][o]-ni] > seuil && check[k][liste[n-1][o]-ni] == 0 )
                        {
                           liste[n] = (int*)realloc(liste[n],(p+1)*sizeof(int));
                           liste[n][p] = liste[n-1][o]-ni;
                           check[k][liste[n][p]] = 1;
                           p++;
                        }
                        if ( (liste[n-1][o]-ni+1) > 0 && ((liste[n-1][o]+1)%ni) != 0 && QR[k][liste[n-1][o]-ni+1] > seuil && check[k][liste[n-1][o]-ni+1] == 0 )
                        {
                           liste[n] = (int*)realloc(liste[n],(p+1)*sizeof(int));
                           liste[n][p] = liste[n-1][o]-ni+1;
                           check[k][liste[n][p]] = 1;
                           p++;
                        }
                     }
                     for ( o = 0; o < p; o++)
                     {
                        if ( (liste[n][o]+1) < (ni*nj) && ((liste[n][o]+1)%ni) != 0 && QR[k][liste[n][o]+1] > seuil && check[k][liste[n][o]+1] == 0 )
                        {
                           liste[n] = (int*)realloc(liste[n],(p+1)*sizeof(int));
                           liste[n][p] = liste[n][o]+1;
                           check[k][liste[n][p]] = 1;
                           p++;
                        }
                        if ( (liste[n][o]-1) > 0 && (liste[n][o]%ni) != 0 && QR[k][liste[n][o]-1] > seuil && check[k][liste[n][o]-1] == 0 )
                        {
                           liste[n] = (int*)realloc(liste[n],(p+1)*sizeof(int));
                           liste[n][p] = liste[n][o]-1;
                           check[k][liste[n][p]] = 1;
                           p++;
                        }
                        if ( (liste[n][o]+ni-1) < (ni*nj) && (liste[n][o]%ni) != 0 && QR[k][liste[n][o]+ni-1] > seuil && check[k][liste[n][o]+ni-1] == 0 )
                        {  
                           liste[n] = (int*)realloc(liste[n],(p+1)*sizeof(int));
                           liste[n][p] = liste[n][o]+ni-1;
                           check[k][liste[n][p]] = 1;
                           p++;
                        }
                        if ( (liste[n][o]+ni) < (ni*nj) && QR[k][liste[n][o]+ni] > seuil && check[k][liste[n][o]+ni] == 0 )
                        {
                           liste[n] = (int*)realloc(liste[n],(p+1)*sizeof(int));
                           liste[n][p] = liste[n][o]+ni;
                           check[k][liste[n][p]] = 1;
                           p++;
                        }
                        if ( (liste[n][o]+ni+1) < (ni*nj) && ((liste[n][o]+1)%ni) != 0 && QR[k][liste[n][o]+ni+1] > seuil && check[k][liste[n][o]+ni+1] == 0 )
                        {
                           liste[n] = (int*)realloc(liste[n],(p+1)*sizeof(int));
                           liste[n][p] = liste[n][o]+ni+1;
                           check[k][liste[n][p]] = 1;
                           p++;
                        }
                        if ( (liste[n][o]-ni-1) > 0 && (liste[n][o]%ni) != 0 && QR[k][liste[n][o]-ni-1] > seuil && check[k][liste[n][o]-ni-1] == 0 )
                        {  
                           liste[n] = (int*)realloc(liste[n],(p+1)*sizeof(int));
                           liste[n][p] = liste[n][o]-ni-1;
                           check[k][liste[n][p]] = 1;
                           p++;
                        }
                        if ( (liste[n][o]-ni) > 0 && QR[k][liste[n][o]-ni] > seuil && check[k][liste[n][o]-ni] == 0 )
                        {
                           liste[n] = (int*)realloc(liste[n],(p+1)*sizeof(int));
                           liste[n][p] = liste[n][o]-ni;
                           check[k][liste[n][p]] = 1;
                           p++;
                        }
                        if ( (liste[n][o]-ni+1) > 0 && ((liste[n][o]+1)%ni) != 0 && QR[k][liste[n][o]-ni+1] > seuil && check[k][liste[n][o]-ni+1] == 0 ) 
                        {
                           liste[n] = (int*)realloc(liste[n],(p+1)*sizeof(int));
                           liste[n][p] = liste[n][o]-ni+1;
                           check[k][liste[n][p]] = 1;
                           p++;
                        }
                     }
                     size = (int*)realloc(size,(n+1)*sizeof(int));
                     size[n] = p;
                     m = p;
                  }                
                  top = k;
                  size1 = (int*)malloc(1*sizeof(int));
                  liste1 = (int**)malloc(1*sizeof(int*));
                  liste1[0] = (int*)malloc(size[0]*sizeof(int));
         
                  size1[0] = size[0];
                  for ( o = 0; o < size1[0]; o++)
                  {
                     liste1[0][o] = liste[0][o];
                  }
                  
                  n = 0;
                  k = start;
                  m = size1[0];
                  while ( m != 0 && k > 0 )
                  {
                     n++;
                     k--;
                     liste1 = (int**)realloc(liste1,(n+1)*sizeof(int*));
                     liste1[n] = (int*)(malloc)(1*sizeof(int));

                     p = 0;
                     for ( o = 0; o < m; o++)
                     {
                        if ( QR[k][liste1[n-1][o]] > seuil && check[k][liste1[n-1][o]] == 0 )
                        {
                           liste1[n] = (int*)realloc(liste1[n],(p+1)*sizeof(int));
                           liste1[n][p] = liste1[n-1][o];
                           check[k][liste1[n][p]] = 1;
                           p++;
                        }
                        if ( (liste1[n-1][o]+1) < (ni*nj) && ((liste1[n-1][o]+1)%ni) != 0 && QR[k][liste1[n-1][o]+1] > seuil && check[k][liste1[n-1][o]+1] == 0 )
                        {
                           liste1[n] = (int*)realloc(liste1[n],(p+1)*sizeof(int));
                           liste1[n][p] = liste1[n-1][o]+1;
                           check[k][liste1[n][p]] = 1;
                           p++;
                        }
                        if ( (liste1[n-1][o]-1) > 0 && (liste1[n-1][o]%ni) != 0 && QR[k][liste1[n-1][o]-1] > seuil && check[k][liste1[n-1][o]-1] == 0 )
                        {        
                           liste1[n] = (int*)realloc(liste1[n],(p+1)*sizeof(int));
                           liste1[n][p] = liste1[n-1][o]-1;
                           check[k][liste1[n][p]] = 1;
                           p++;
                        }
                        if ( (liste1[n-1][o]+ni-1) < (ni*nj) && (liste1[n-1][o]%ni) != 0 && QR[k][liste1[n-1][o]+ni-1] > seuil && check[k][liste1[n-1][o]+ni-1] == 0 )
                        {
                           liste1[n] = (int*)realloc(liste1[n],(p+1)*sizeof(int));
                           liste1[n][p] = liste1[n-1][o]+ni-1;
                           check[k][liste1[n][p]] = 1;
                           p++;
                        }
                        if ( (liste1[n-1][o]+ni) < (ni*nj) && QR[k][liste1[n-1][o]+ni] > seuil && check[k][liste1[n-1][o]+ni] == 0 )
                        {
                           liste1[n] = (int*)realloc(liste1[n],(p+1)*sizeof(int));
                           liste1[n][p] = liste1[n-1][o]+ni;
                           check[k][liste1[n][p]] = 1;
                           p++;
                        }
                        if ( (liste1[n-1][o]+ni+1) < (ni*nj) && ((liste1[n-1][o]+1)%ni) != 0 && QR[k][liste1[n-1][o]+ni+1] > seuil && check[k][liste1[n-1][o]+ni+1] == 0 )
                        {
                           liste1[n] = (int*)realloc(liste1[n],(p+1)*sizeof(int));
                           liste1[n][p] = liste1[n-1][o]+ni+1;
                           check[k][liste1[n][p]] = 1;
                           p++;
                        }
                        if ( (liste1[n-1][o]-ni-1) > 0 && (liste1[n-1][o]%ni) != 0 && QR[k][liste1[n-1][o]-ni-1] > seuil && check[k][liste1[n-1][o]-ni-1] == 0 )
                        {
                           liste1[n] = (int*)realloc(liste1[n],(p+1)*sizeof(int));
                           liste1[n][p] = liste1[n-1][o]-ni-1;
                           check[k][liste1[n][p]] = 1;
                           p++;
                        }
                        if ( (liste1[n-1][o]-ni) > 0 && QR[k][liste1[n-1][o]-ni] > seuil && check[k][liste1[n-1][o]-ni] == 0 )
                        {
                           liste1[n] = (int*)realloc(liste1[n],(p+1)*sizeof(int));
                           liste1[n][p] = liste1[n-1][o]-ni;
                           check[k][liste1[n][p]] = 1;
                           p++;
                        }
                        if ( (liste1[n-1][o]-ni+1) > 0 && ((liste1[n-1][o]+1)%ni) != 0 && QR[k][liste1[n-1][o]-ni+1] > seuil && check[k][liste1[n-1][o]-ni+1] == 0 )
                        {
                           liste1[n] = (int*)realloc(liste1[n],(p+1)*sizeof(int));
                           liste1[n][p] = liste1[n-1][o]-ni+1;
                           check[k][liste1[n][p]] = 1;
                           p++;
                        }
                     }
                     for ( o = 0; o < p; o++)
                     {
                        if ( (liste1[n][o]+1) < (ni*nj) && ((liste1[n][o]+1)%ni) != 0 && QR[k][liste1[n][o]+1] > seuil && check[k][liste1[n][o]+1] == 0 )
                        {
                           liste1[n] = (int*)realloc(liste1[n],(p+1)*sizeof(int));
                           liste1[n][p] = liste1[n][o]+1;
                           check[k][liste1[n][p]] = 1;
                           p++;
                        }
                        if ( (liste1[n][o]-1) > 0 && (liste1[n][o]%ni) != 0 && QR[k][liste1[n][o]-1] > seuil && check[k][liste1[n][o]-1] == 0 )
                        {
                           liste1[n] = (int*)realloc(liste1[n],(p+1)*sizeof(int));
                           liste1[n][p] = liste1[n][o]-1;
                           check[k][liste1[n][p]] = 1;
                           p++;
                        }
                        if ( (liste1[n][o]+ni-1) < (ni*nj) && (liste1[n][o]%ni) != 0 && QR[k][liste1[n][o]+ni-1] > seuil && check[k][liste1[n][o]+ni-1] == 0 )
                        {  
                           liste1[n] = (int*)realloc(liste1[n],(p+1)*sizeof(int));
                           liste1[n][p] = liste1[n][o]+ni-1;
                           check[k][liste1[n][p]] = 1;
                           p++;
                        }
                        if ( (liste1[n][o]+ni) < (ni*nj) && QR[k][liste1[n][o]+ni] > seuil && check[k][liste1[n][o]+ni] == 0 )
                        {
                           liste1[n] = (int*)realloc(liste1[n],(p+1)*sizeof(int));
                           liste1[n][p] = liste1[n][o]+ni;
                           check[k][liste1[n][p]] = 1;
                           p++;
                        }
                        if ( (liste1[n][o]+ni+1) < (ni*nj) && ((liste1[n][o]+1)%ni) != 0 && QR[k][liste1[n][o]+ni+1] > seuil && check[k][liste1[n][o]+ni+1] == 0 )
                        {
                           liste1[n] = (int*)realloc(liste1[n],(p+1)*sizeof(int));
                           liste1[n][p] = liste1[n][o]+ni+1;
                           check[k][liste1[n][p]] = 1;
                           p++;
                        }        
                        if ( (liste1[n][o]-ni-1) > 0 && (liste1[n][o]%ni) != 0 && QR[k][liste1[n][o]-ni-1] > seuil && check[k][liste1[n][o]-ni-1] == 0 )
                        {         
                           liste1[n] = (int*)realloc(liste1[n],(p+1)*sizeof(int));
                           liste1[n][p] = liste1[n][o]-ni-1;
                           check[k][liste1[n][p]] = 1;
                           p++;
                        }
                        if ( (liste1[n][o]-ni) > 0 && QR[k][liste1[n][o]-ni] > seuil && check[k][liste1[n][o]-ni] == 0 )
                        {
                           liste1[n] = (int*)realloc(liste1[n],(p+1)*sizeof(int));
                           liste1[n][p] = liste1[n][o]-ni;
                           check[k][liste1[n][p]] = 1;
                           p++;
                        }
                        if ( (liste1[n][o]-ni+1) > 0 && ((liste1[n][o]+1)%ni) != 0 && QR[k][liste1[n][o]-ni+1] > seuil && check[k][liste1[n][o]-ni+1] == 0 ) 
                        {
                           liste1[n] = (int*)realloc(liste1[n],(p+1)*sizeof(int));
                           liste1[n][p] = liste1[n][o]-ni+1;
                           check[k][liste1[n][p]] = 1;
                           p++;
                        }
                     }   
                     size1 = (int*)realloc(size1,(n+1)*sizeof(int));
                     size1[n] = p;
                     m = p;
                  }
                  base = k;
                  size2 = (int*)malloc((top-base-1)*sizeof(int));
                  liste2 = (int**)malloc((top-base-1)*sizeof(int*));
                  for ( o = (start-base-1); o > 0; o--)
                  {
                     liste2[start-base-1-o] = (int*)malloc(size1[o]*sizeof(int));
                     size2[start-base-1-o] = size1[o];
                     for ( p = 0; p < size1[o]; p++ )
                     {
                        liste2[start-base-1-o][p] = liste1[o][p];
                     }
                  }
                  for ( o = 0; o < (top-start); o++)
                  {
                     liste2[start-base-1+o] = (int*)malloc(size[o]*sizeof(int));
                     size2[start-base-1+o] = size[o];

                     for ( p = 0; p < size[o]; p++ )
                     {
                        liste2[start-base-1+o][p] = liste[o][p];
                     }
                  }

                  for ( o = 0; o < (top-start+1); o++)
                  {
                     free(liste[o]);
                  }
                  for ( o = 0; o < (start-base+1); o++)
                  {
                     free(liste1[o]);
                  }
                  free(size);
                  free(size1);
                  free(liste);
                  free(liste1);
                  size = NULL;
                  size1 = NULL;
                  liste = NULL;
                  liste1 = NULL;

                  if ( (GZ[top-1][i] - GZ[base+1][i]) > 200. )
                  {
                     p = 0;
                     for ( o = 0; o < (top-base-1); o++ )
                     {
                        for ( p = 0; p < size2[o]; p++ )
                        {
                           mesocyclone[liste2[o][p]] =  1;
                           base_mesocyclone[liste2[o][p]] = GZ[base+1][i] - GZ[0][i];
                        }
                     }
                  }
                  
                  for ( o = 0; o < (top-base-1); o++)
                  {
                     free(liste2[o]);
                  }
                  free(size2);
                  free(liste2);
                  size2 = NULL;
                  liste2 = NULL;
               }
            }
         }
      } 
   } 
   
   //-------------------------
   // Liberer la memoire
   //-------------------------
   for ( k = 0; k < nk; k++ )
   {
      free(check[k]);
   }
   free(check);
   check = NULL;
   
   //----------------------------
   // Retourner les variables
   //----------------------------
   *Mesocyclone = mesocyclone;
   *BaseMesoyclone = base_mesocyclone;
   
   return SUCCESS;
}
                       
                       

/*!
=============================================================================
\brief Calcul les composantes du vent


\param[in]  UV         : Tableau du module du vent (noeuds)
\param[in]  WD         : Tableau de la direction du vent (deg meteo)
\param[in]  ni         : Dimension en x
\param[in]  nj         : Dimension en y
\param[in]  nbr_niveau : Dimension en z
\param[out] UU         : Tableau de la composante en x du vent (noeuds)
\param[out] VV         : Tableau de la composante en y du vent (noeuds)

\date 15 juin 2010
\version 0.1
\author Anna-Belle Filion

-----------------------------------------------------------------------------
*/
int CalculerUU_VV( float **UV, float **WD, int ni, int nj, int nbr_niveau,
                   float **UU, float **VV)
{
   int i, k;
   
   double wd_math;
   
   //-----------------------------------
   // Calculer les composantes du vent
   //-----------------------------------
   for ( k = 0; k < nbr_niveau; k++)
   {
      for ( i = 0; i < (ni*nj); i++)
      {
         if ( WD[k][i] == 360.)
            WD[k][i] = 0.;
         
         wd_math = 270. - (double)WD[k][i];
         
         if ( wd_math < 0. )
            wd_math = wd_math + 360.;
         
         UU[k][i] = (float)(UV[k][i] * cos(wd_math * DEGRE_A_RADIAN));
         VV[k][i] = (float)(UV[k][i] * sin(wd_math * DEGRE_A_RADIAN));
      }
   }
   
   //--------------------------
   // Retourner les resultats
   //--------------------------
   
   return SUCCESS;
}


/*!
=============================================================================
\brief Trouver la convergence d'humidite maximum dans une colonne


\param[in]  HU                     : Tableau d'humidite specifique (kg/kg)
\param[in]  DD                     : Tableau de la divergence (s^(-1))
\param[in]  GZ                     : Tableau des hauteurs geopotentielle (DAM)
\param[in]  ni                     : Dimension en x
\param[in]  nj                     : Dimension en y
\param[in]  nbr_niveau             : Dimension en z
\param[out] MoistureConvergenceX   : Vecteur de convergence maximum d'humidite 
                                     (g/kg *s-1)
\param[out] HauteurMoistureConvX   : Vecteur de la hauteur de la congervence 
                                     maximum d'humidite (m) 

\date 05 mai 2010
\version 0.1
\author Anna-Belle Filion

-----------------------------------------------------------------------------
*/
int ConvergenceHumidite( float **HU, float **DD, float **GZ, int ni, int nj, 
                         int nbr_niveau, float **MoistureConvergenceX, 
                         float **HauteurMoistureConvX)
{
   int i, k;
   
   float tmp;
   
   float *Conv_Moist = (float*)NULL;
   float *Haut_Conv_Moist = (float*)NULL;
   
   Conv_Moist = (float*)malloc(ni*nj*sizeof(float));
   Haut_Conv_Moist = (float*)malloc(ni*nj*sizeof(float));
   
   for ( i = 0; i < (ni*nj); i++)
   {
      Conv_Moist[i] = 0;
      for ( k = 0; k < nbr_niveau; k++)
      {
         if ( DD[k][i] < 0 )
            tmp = (float)( HU[k][i] * 1000. * fabs((double)(DD[k][i])));
         else 
            tmp = 0;

         if ( tmp > Conv_Moist[i] )
         {
            Conv_Moist[i] = tmp;
            Haut_Conv_Moist[i] = (GZ[k][i] - GZ[0][i]) * 10.;
         }
      }
   }

   //--------------------------
   // Retourner les resultats
   //--------------------------
   *MoistureConvergenceX = Conv_Moist;
   *HauteurMoistureConvX = Haut_Conv_Moist;
   
   return SUCCESS;
}




/*!
==============================================================================
\brief Detecter les tornades


\param[in]  QR          : Tableau du tourbillon relatif (s^-1)
\param[in]  GZ          : Tableau des hauteurs geopotentielle (Dam)
\param[in]  ni          : Dimension en x
\param[in]  nj          : Dimension en y
\param[in]  nk          : Dimension en z
\param[out] Tornade     : Vecteur de tornades possible ( 1 => Oui / 0 => Non )

\date 05 mai 2010
\version 0.1
\author Anna-Belle Filion

------------------------------------------------------------------------------
*/
int DetecterTornade(float **QR, float **GZ, int ni, int nj, int nk, float seuil,
                    float **tornade)
{
   int i, j, k;
   
   float *Tornade = (float*)NULL;
   
   //------------------------------
   // Allocation de memoire
   //------------------------------
   Tornade = (float*)malloc(ni*nj*sizeof(float));
   
   for ( i = 0; i < (ni*nj); i++)
   {
      k = 0;
      j = 0;
      Tornade[i] = 0;
      while ( (GZ[k][i] - GZ[0][i]) <= 200. )
      {
         if ( QR[k][i] >= seuil )
            j++;
         else 
            j = 0;
         
         if ( j == 3 )
            Tornade[i] = 1;
         k++;
      }
   }
   
   //--------------------------
   // Retourner les resultats
   //--------------------------
   *tornade = Tornade;
   
   return SUCCESS;
}


/*!
===============================================================================
\brief Detecter les couplets de tourbillon

\param[in]  QR      : Tableau de tourbillon relatif (s^-1)
\param[in]  GZ      : Tableau des hauteurs geopotentielle (DAM)
\param[in]  ni      : Dimension en x
\param[in]  nj      : Dimension en y
\param[in]  nk      : Dimension en z
\param[out] couplet : Vecteur de couplet de tourbillon


\date 17 main 2010
\version 0.0
\author Anna-Belle Filion

--------------------------------------------------------------------------------
*/
int DetecterCoupletTourbillion(float **QR, float *MaxR, float **GZ, int ni, int nj,
                               int nk, float seuil[2][2], float **couplet_tourbillon,
                               float **couplet_positif, float **couplet_negatif)
{
   int i, j, k, l, m, n, o, p, q;
   
   int CoupletPotentiel, base;
   int hauteurMinNV, hauteurMaxNV, hauteurMinPV, hauteurMaxPV;
   
   int *size = (int*)NULL;
   int **liste = (int**)NULL;
   int **checkNV = (int**)NULL;
   int **checkPV = (int**)NULL;
   int *checkNV_ = (int*)NULL;
   int *checkPV_ = (int*)NULL;
   
   float seuil1 = 0, seuil3 = 0, seuil4 = 0, seuil5 = 0;
   
   float *Negative_Vorticity = (float*)NULL;
   float *Positive_Vorticity = (float*)NULL;
   float *Vorticity_couplet  = (float*)NULL;

   //---------------------------------------------
   // Attribution des seuils selon la resolution
   //---------------------------------------------
   seuil1 = seuil[1][0];
   seuil3 = seuil[0][0];
   seuil4 = seuil[1][1];
   seuil5 = seuil[0][1];

   //---------------------------
   // Allocation de memoire
   //---------------------------
   checkNV = (int**)malloc(nk*sizeof(int*));
   checkPV = (int**)malloc(nk*sizeof(int*));
   
   for ( i = 0; i < nk; i++ )
   {
      checkNV[i] = (int*)malloc(ni*nj*sizeof(int));
      checkPV[i] = (int*)malloc(ni*nj*sizeof(int));
   }
   
   checkNV_ = (int*)malloc(ni*nj*sizeof(int));
   checkPV_ = (int*)malloc(ni*nj*sizeof(int));
   
   Negative_Vorticity = (float*)malloc(ni*nj*sizeof(float));
   Positive_Vorticity = (float*)malloc(ni*nj*sizeof(float));
   Vorticity_couplet = (float*)malloc(ni*nj*sizeof(float));
   
   //--------------------------------
   // Initialisation des variables
   //--------------------------------
   for ( k = 0; k < nk; k++ )
   {
      for ( i = 0; i < (ni*nj); i++)
      {
         checkNV[k][i] = 0;
         checkPV[k][i] = 0;
         if ( k == 0 )
         {
            Negative_Vorticity[i] = 0;
            Positive_Vorticity[i] = 0;
            Vorticity_couplet[i]  = 0;
         }
      }
   }
   
   //--------------------------------
   // Trouver les vortex positif
   //--------------------------------
   for ( i = 0; i < (ni*nj); i++)
   {
      l = 0;
      j = 0;
      for ( k = 0; k < nk; k++)
      {
         if ( ((GZ[k][i] - GZ[0][i]) < 1200.) && ((GZ[k][i] - GZ[0][i]) > 100.) )
         {
            m = 0;
            n = 0;
            if ( QR[k][i] > seuil1 && MaxR[i] > 40. && checkPV[k][i] == 0 )
            {                
               base = k;
               for ( k = base; k > 0; k-- )
               {
                  if ( QR[k][i] > seuil4 && checkPV[k][i] == 0 )
                     break;
               }

               size = (int*)malloc(1*sizeof(int));
               liste = (int**)malloc(1*sizeof(int*));
               liste[n] = (int*)(malloc)(1*sizeof(int));
               liste[n][m] = i;
               checkPV[k][i] = 1;
	       base = k;
               m++;

               if ( (i+1) < (ni*nj) && ((i+1)%ni) != 0 && QR[k][i+1] > seuil4 && checkPV[k][i+1] == 0 )
               {
                  liste[n] = (int*)realloc(liste[n],(m+1)*sizeof(int));
                  liste[n][m] = i+1;
                  checkPV[k][i+1] = 1;
                  m++;
               }
               if ( (i-1) > 0 && (i%ni) != 0 && QR[k][i-1] > seuil4 && checkPV[k][i-1] == 0 )
               {
                  liste[n] = (int*)realloc(liste[n],(m+1)*sizeof(int));
                  liste[n][m] = i-1;
                  checkPV[k][i-1] = 1;
                  m++;   
               }
               if ( (i+ni-1) < (ni*nj) && (i%ni) != 0 && QR[k][i+ni-1] > seuil4 && checkPV[k][i+ni-1] == 0 )
               {
                  liste[n] = (int*)realloc(liste[n],(m+1)*sizeof(int));
                  liste[n][m] = i+ni-1;
                  checkPV[k][i+ni-1] = 1;
                  m++;
               }
               if ( (i+ni) < (ni*nj) && QR[k][i+ni] > seuil4 && checkPV[k][i+ni] == 0 )
               {
                  liste[n] = (int*)realloc(liste[n],(m+1)*sizeof(int));
                  liste[n][m] = i+ni;
                  checkPV[k][i+ni] = 1;
                  m++;   
               }
               if ( (i+ni+1) < (ni*nj) && ((i+1)%ni) != 0 && QR[k][i+ni+1] > seuil4 && checkPV[k][i+ni+1] == 0 )
               {
                  liste[n] = (int*)realloc(liste[n],(m+1)*sizeof(int));
                  liste[n][m] = i+ni+1;
                  checkPV[k][i+ni+1] = 1;
                  m++;   
               }
               if ( (i-ni-1) > 0 && (i%ni) != 0 && QR[k][i-ni-1] > seuil4 && checkPV[k][i-ni-1] == 0 )
               {
                  liste[n] = (int*)realloc(liste[n],(m+1)*sizeof(int));
                  liste[n][m] = i-ni-1;
                  checkPV[k][i-ni-1] = 1;
                  m++;   
               }
               if ( (i-ni) > 0 && QR[k][i-ni] > seuil4 && checkPV[k][i-ni] == 0 )
               {
                  liste[n] = (int*)realloc(liste[n],(m+1)*sizeof(int));
                  liste[n][m] = i-ni;
                  checkPV[k][i-ni] = 1;
                  m++;  
               }
               if ( (i-ni+1) > 0 && ((i+1)%ni) != 0 && QR[k][i-ni+1] > seuil4 && checkPV[k][i-ni+1] == 0 ) 
               {
                  liste[n] = (int*)realloc(liste[n],(m+1)*sizeof(int));
                  liste[n][m] = i-ni+1;
                  checkPV[k][i-ni+1] = 1;
                  m++;   
               }
               for ( o = 1; o < m; o++)
               {
		  if ( (liste[n][o]+1) < (ni*nj) && ((liste[n][o]+1)%ni) != 0 && QR[k][liste[n][o]+1] > seuil4 && checkPV[k][liste[n][o]+1] == 0 )
                  {
                     liste[n] = (int*)realloc(liste[n],(m+1)*sizeof(int));
                     liste[n][m] = liste[n][o]+1;
                     checkPV[k][liste[n][m]] = 1;
                     m++;
                  }
                  if ( (liste[n][o]-1) > 0 && (liste[n][o]%ni) != 0 && QR[k][liste[n][o]-1] > seuil4 && checkPV[k][liste[n][o]-1] == 0 )
                  {
                     liste[n] = (int*)realloc(liste[n],(m+1)*sizeof(int));
                     liste[n][m] = liste[n][o]-1;
                     checkPV[k][liste[n][m]] = 1;
                     m++;
                  }
                  if ( (liste[n][o]+ni-1) < (ni*nj) && (liste[n][o]%ni) != 0 && QR[k][liste[n][o]+ni-1] > seuil4 && checkPV[k][liste[n][o]+ni-1] == 0 )
		  {
                     liste[n] = (int*)realloc(liste[n],(m+1)*sizeof(int));
                     liste[n][m] = liste[n][o]+ni-1;
                     checkPV[k][liste[n][m]] = 1;
                     m++;
                  }
                  if ( (liste[n][o]+ni) < (ni*nj) && QR[k][liste[n][o]+ni] > seuil4 && checkPV[k][liste[n][o]+ni] == 0 )
                  {
                     liste[n] = (int*)realloc(liste[n],(m+1)*sizeof(int));
                     liste[n][m] = liste[n][o]+ni;
                     checkPV[k][liste[n][m]] = 1;
                     m++;
                  }
                  if ( (liste[n][o]+ni+1) < (ni*nj) && ((liste[n][o]+1)%ni) != 0 && QR[k][liste[n][o]+ni+1] > seuil4 && checkPV[k][liste[n][o]+ni+1] == 0 )
                  {
                     liste[n] = (int*)realloc(liste[n],(m+1)*sizeof(int));
                     liste[n][m] = liste[n][o]+ni+1;
                     checkPV[k][liste[n][m]] = 1;
                     m++;
                  }
                  if ( (liste[n][o]-ni-1) > 0 && (liste[n][o]%ni) != 0 && QR[k][liste[n][o]-ni-1] > seuil4 && checkPV[k][liste[n][o]-ni-1] == 0 )
                  {  
                     liste[n] = (int*)realloc(liste[n],(m+1)*sizeof(int));
                     liste[n][m] = liste[n][o]-ni-1;
                     checkPV[k][liste[n][m]] = 1;
                     m++;
                  }
                  if ( (liste[n][o]-ni) > 0 && QR[k][liste[n][o]-ni] > seuil4 && checkPV[k][liste[n][o]-ni] == 0 )
                  {
                     liste[n] = (int*)realloc(liste[n],(m+1)*sizeof(int));
                     liste[n][m] = liste[n][o]-ni;
                     checkPV[k][liste[n][m]] = 1;
                     m++;
                  }
                  if ( (liste[n][o]-ni+1) > 0 && ((liste[n][o]+1)%ni) != 0 && QR[k][liste[n][o]-ni+1] > seuil4 && checkPV[k][liste[n][o]-ni+1] == 0 ) 
                  {
                     liste[n] = (int*)realloc(liste[n],(m+1)*sizeof(int));
                     liste[n][m] = liste[n][o]-ni+1;
                     checkPV[k][liste[n][m]] = 1;
                     m++;
                  }                 
               }
                           
               size[n] = m;
               while ( (GZ[k][i]  - GZ[0][i]) < 1200. && m != 0 && k < (nk-1) )
               {
                  n++;
		  k++;
                  liste = (int**)realloc(liste,(n+1)*sizeof(int*));
                  liste[n] = (int*)(malloc)(1*sizeof(int));
                         
                  p = 0;
                  for ( o = 0; o < m; o++)
                  {
                     if ( QR[k][liste[n-1][o]] > seuil4 && checkPV[k][liste[n-1][o]] == 0 )
                     {
                        liste[n] = (int*)realloc(liste[n],(p+1)*sizeof(int));
                        liste[n][p] = liste[n-1][o];
                        checkPV[k][liste[n][p]] = 1;
                        p++;
                     }
                     if ( (liste[n-1][o]+1) < (ni*nj) && ((liste[n-1][o]+1)%ni) != 0 && QR[k][liste[n-1][o]+1] > seuil4 && checkPV[k][liste[n-1][o]+1] == 0 )
                     {
                        liste[n] = (int*)realloc(liste[n],(p+1)*sizeof(int));
                        liste[n][p] = liste[n-1][o]+1;
                        checkPV[k][liste[n][p]] = 1;
                        p++;
                     }
                     if ( (liste[n-1][o]-1) > 0 && (liste[n-1][o]%ni) != 0 && QR[k][liste[n-1][o]-1] > seuil4 && checkPV[k][liste[n-1][o]-1] == 0 )
                     {
                        liste[n] = (int*)realloc(liste[n],(p+1)*sizeof(int));
                        liste[n][p] = liste[n-1][o]-1;
                        checkPV[k][liste[n][p]] = 1;
                        p++;
                     }
                     if ( (liste[n-1][o]+ni-1) < (ni*nj) && (liste[n-1][o]%ni) != 0 && QR[k][liste[n-1][o]+ni-1] > seuil4 && checkPV[k][liste[n-1][o]+ni-1] == 0 )
                     {
                        liste[n] = (int*)realloc(liste[n],(p+1)*sizeof(int));
                        liste[n][p] = liste[n-1][o]+ni-1;
                        checkPV[k][liste[n][p]] = 1;
                        p++;
                     }
                     if ( (liste[n-1][o]+ni) < (ni*nj) && QR[k][liste[n-1][o]+ni] > seuil4 && checkPV[k][liste[n-1][o]+ni] == 0 )
                     {
                        liste[n] = (int*)realloc(liste[n],(p+1)*sizeof(int));
                        liste[n][p] = liste[n-1][o]+ni;
                        checkPV[k][liste[n][p]] = 1;
                        p++;
                     }
                     if ( (liste[n-1][o]+ni+1) < (ni*nj) && ((liste[n-1][o]+1)%ni) != 0 && QR[k][liste[n-1][o]+ni+1] > seuil4 && checkPV[k][liste[n-1][o]+ni+1] == 0 )
                     {
                        liste[n] = (int*)realloc(liste[n],(p+1)*sizeof(int));
                        liste[n][p] = liste[n-1][o]+ni+1;
                        checkPV[k][liste[n][p]] = 1;
                        p++;
                     }
                     if ( (liste[n-1][o]-ni-1) > 0 && (liste[n-1][o]%ni) != 0 && QR[k][liste[n-1][o]-ni-1] > seuil4 && checkPV[k][liste[n-1][o]-ni-1] == 0 )
                     {
                        liste[n] = (int*)realloc(liste[n],(p+1)*sizeof(int));
                        liste[n][p] = liste[n-1][o]-ni-1;
                        checkPV[k][liste[n][p]] = 1;
                        p++;
                     }
                     if ( (liste[n-1][o]-ni) > 0 && QR[k][liste[n-1][o]-ni] > seuil4 && checkPV[k][liste[n-1][o]-ni] == 0 )
                     {
                        liste[n] = (int*)realloc(liste[n],(p+1)*sizeof(int));
                        liste[n][p] = liste[n-1][o]-ni;
                        checkPV[k][liste[n][p]] = 1;
                        p++;
                     }
                     if ( (liste[n-1][o]-ni+1) > 0 && ((liste[n-1][o]+1)%ni) != 0 && QR[k][liste[n-1][o]-ni+1] > seuil4 && checkPV[k][liste[n-1][o]-ni+1] == 0 )
                     {
                        liste[n] = (int*)realloc(liste[n],(p+1)*sizeof(int));
                        liste[n][p] = liste[n-1][o]-ni+1;
                        checkPV[k][liste[n][p]] = 1;
                        p++;
                     }
                  }
                  for ( o = 0; o < p; o++)
                  {
                     if ( (liste[n][o]+1) < (ni*nj) && ((liste[n][o]+1)%ni) != 0 && QR[k][liste[n][o]+1] > seuil4 && checkPV[k][liste[n][o]+1] == 0 )
                     {
                        liste[n] = (int*)realloc(liste[n],(p+1)*sizeof(int));
                        liste[n][p] = liste[n][o]+1;
                        checkPV[k][liste[n][p]] = 1;
                        p++;
                     }
                     if ( (liste[n][o]-1) > 0 && (liste[n][o]%ni) != 0 && QR[k][liste[n][o]-1] > seuil4 && checkPV[k][liste[n][o]-1] == 0 )
                     {
                        liste[n] = (int*)realloc(liste[n],(p+1)*sizeof(int));
                        liste[n][p] = liste[n][o]-1;
                        checkPV[k][liste[n][p]] = 1;
                        p++;
                     }
                     if ( (liste[n][o]+ni-1) < (ni*nj) && (liste[n][o]%ni) != 0 && QR[k][liste[n][o]+ni-1] > seuil4 && checkPV[k][liste[n][o]+ni-1] == 0 )
                     {  
                        liste[n] = (int*)realloc(liste[n],(p+1)*sizeof(int));
                        liste[n][p] = liste[n][o]+ni-1;
                        checkPV[k][liste[n][p]] = 1;
                        p++;
                     }
                     if ( (liste[n][o]+ni) < (ni*nj) && QR[k][liste[n][o]+ni] > seuil4 && checkPV[k][liste[n][o]+ni] == 0 )
                     {
                        liste[n] = (int*)realloc(liste[n],(p+1)*sizeof(int));
                        liste[n][p] = liste[n][o]+ni;
                        checkPV[k][liste[n][p]] = 1;
                        p++;
                     }
                     if ( (liste[n][o]+ni+1) < (ni*nj) && ((liste[n][o]+1)%ni) != 0 && QR[k][liste[n][o]+ni+1] > seuil4 && checkPV[k][liste[n][o]+ni+1] == 0 )
                     {
                        liste[n] = (int*)realloc(liste[n],(p+1)*sizeof(int));
                        liste[n][p] = liste[n][o]+ni+1;
                        checkPV[k][liste[n][p]] = 1;
                        p++;
                     }
                     if ( (liste[n][o]-ni-1) > 0 && (liste[n][o]%ni) != 0 && QR[k][liste[n][o]-ni-1] > seuil4 && checkPV[k][liste[n][o]-ni-1] == 0 )
                     {  
                        liste[n] = (int*)realloc(liste[n],(p+1)*sizeof(int));
                        liste[n][p] = liste[n][o]-ni-1;
                        checkPV[k][liste[n][p]] = 1;
                        p++;
                     }
                     if ( (liste[n][o]-ni) > 0 && QR[k][liste[n][o]-ni] > seuil4 && checkPV[k][liste[n][o]-ni] == 0 )
                     {
                        liste[n] = (int*)realloc(liste[n],(p+1)*sizeof(int));
                        liste[n][p] = liste[n][o]-ni;
                        checkPV[k][liste[n][p]] = 1;
                        p++;
                     }
                     if ( (liste[n][o]-ni+1) > 0  && ((liste[n][o]+1)%ni) != 0 && QR[k][liste[n][o]-ni+1] > seuil4 && checkPV[k][liste[n][o]-ni+1] == 0 ) 
                     {
                        liste[n] = (int*)realloc(liste[n],(p+1)*sizeof(int));
                        liste[n][p] = liste[n][o]-ni+1;
                        checkPV[k][liste[n][p]] = 1;
                        p++;
                     }
                  }
                  size = (int*)realloc(size,(n+1)*sizeof(int));
                  size[n] = p;
                  m = p;
                  //k++;
               }
               //n--;
               if ( (GZ[k][i]  - GZ[base][i]) > 200. )
               {
                  for ( o = 0; o < n; o++)
                  {
                     for ( q = 0; q < size[o]; q++)
                     {
                        Positive_Vorticity[liste[o][q]] = 1;
                     }
                  }
               }
               for ( o = 0; o < (n+1); o++)
               {
                  free (liste[o]);
               }
               free(size);
               free(liste);
               size = NULL;
               liste = NULL;
            }
            //break;
         }
      }
   } 
   
   for ( i = 0; i < (ni*nj); i++)
   {
      l = 0;
      j = 0;
      for ( k = 0; k < nk; k++)
      {
	 if ( ((GZ[k][i] - GZ[0][i]) < 1200.) && ((GZ[k][i] - GZ[0][i]) > 100.) )
         {
            m = 0;
            n = 0;
            if ( QR[k][i] < seuil3 && checkNV[k][i] == 0 )
            {
               for ( k = 0; k < nk; k++ )
               {
                  if ( QR[k][i] < seuil5 && checkNV[k][i] == 0 )
                     break;
               }
               
               size = (int*)malloc(1*sizeof(int));
               liste = (int**)malloc(1*sizeof(int*));
               liste[n] = (int*)(malloc)(1*sizeof(int));
               liste[n][m] = i;
               checkNV[k][i] = 1;
	       base = k;
               m++;

               if ( (i+1) < (ni*nj) && ((i+1)%ni) != 0 && QR[k][i+1] < seuil5 && checkNV[k][i+1] == 0 )
               {
                  liste[n] = (int*)realloc(liste[n],(m+1)*sizeof(int));
                  liste[n][m] = i+1;
                  checkNV[k][i+1] = 1;
                  m++;
               }
               if ( (i-1) > 0 && (i%ni) != 0 && QR[k][i-1] < seuil5 && checkNV[k][i-1] == 0 )
               {
                  liste[n] = (int*)realloc(liste[n],(m+1)*sizeof(int));
                  liste[n][m] = i-1;
                  checkNV[k][i-1] = 1;
                  m++;   
               }
               if ( (i+ni-1) < (ni*nj) && (i%ni) != 0 && QR[k][i+ni-1] < seuil5 && checkNV[k][i+ni-1] == 0 )
               {
                  liste[n] = (int*)realloc(liste[n],(m+1)*sizeof(int));
                  liste[n][m] = i+ni-1;
                  checkNV[k][i+ni-1] = 1;
                  m++;
               }
               if ( (i+ni) < (ni*nj) && QR[k][i+ni] < seuil5 && checkNV[k][i+ni] == 0 )
               {
                  liste[n] = (int*)realloc(liste[n],(m+1)*sizeof(int));
                  liste[n][m] = i+ni;
                  checkNV[k][i+ni] = 1;
                  m++;   
               }
               if ( (i+ni+1) < (ni*nj) && ((i+1)%ni) != 0 && QR[k][i+ni+1] < seuil5 && checkNV[k][i+ni+1] == 0 )
               {
                  liste[n] = (int*)realloc(liste[n],(m+1)*sizeof(int));
                  liste[n][m] = i+ni+1;
                  checkNV[k][i+ni+1] = 1;
                  m++;   
               }
               if ( (i-ni-1) > 0 && (i%ni) != 0 && QR[k][i-ni-1] < seuil5 && checkNV[k][i-ni-1] == 0 )
               {
                  liste[n] = (int*)realloc(liste[n],(m+1)*sizeof(int));
                  liste[n][m] = i-ni-1;
                  checkNV[k][i-ni-1] = 1;
                  m++;   
               }
               if ( (i-ni) > 0 && QR[k][i-ni] < seuil5 && checkNV[k][i-ni] == 0 )
               {
                  liste[n] = (int*)realloc(liste[n],(m+1)*sizeof(int));
                  liste[n][m] = i-ni;
                  checkNV[k][i-ni] = 1;
                  m++;  
               }
               if ( (i-ni+1) > 0 && ((i+1)%ni) != 0 && QR[k][i-ni+1] < seuil5 && checkNV[k][i-ni+1] == 0 ) 
               {
                  liste[n] = (int*)realloc(liste[n],(m+1)*sizeof(int));
                  liste[n][m] = i-ni+1;
                  checkNV[k][i-ni+1] = 1;
                  m++;   
               }
               for ( o = 1; o < m; o++)
               {
                  if ( (liste[n][o]+1) < (ni*nj) && ((liste[n][o]+1)%ni) != 0 && QR[k][liste[n][o]+1] < seuil5 && checkNV[k][liste[n][o]+1] == 0 )
                  {
                     liste[n] = (int*)realloc(liste[n],(m+1)*sizeof(int));
                     liste[n][m] = liste[n][o]+1;
                     checkNV[k][liste[n][m]] = 1;
                     m++;
                  }
                  if ( (liste[n][o]-1) > 0 && (liste[n][o]%ni) != 0 && QR[k][liste[n][o]-1] < seuil5 && checkNV[k][liste[n][o]-1] == 0 )
                  {
                     liste[n] = (int*)realloc(liste[n],(m+1)*sizeof(int));
                     liste[n][m] = liste[n][o]-1;
                     checkNV[k][liste[n][m]] = 1;
                     m++;
                  }
                  if ( (liste[n][o]+ni-1) < (ni*nj) && (liste[n][o]%ni) != 0 && QR[k][liste[n][o]+ni-1] < seuil5 && checkNV[k][liste[n][o]+ni-1] == 0 )
                  {
                     liste[n] = (int*)realloc(liste[n],(m+1)*sizeof(int));
                     liste[n][m] = liste[n][o]+ni-1;
                     checkNV[k][liste[n][m]] = 1;
                     m++;
                  }
                  if ( (liste[n][o]+ni) < (ni*nj) && QR[k][liste[n][o]+ni] < seuil5 && checkNV[k][liste[n][o]+ni] == 0 )
                  {
                     liste[n] = (int*)realloc(liste[n],(m+1)*sizeof(int));
                     liste[n][m] = liste[n][o]+ni;
                     checkNV[k][liste[n][m]] = 1;
                     m++;
                  }
                  if ( (liste[n][o]+ni+1) < (ni*nj) && ((liste[n][o]+1)%ni) != 0 && QR[k][liste[n][o]+ni+1] < seuil5 && checkNV[k][liste[n][o]+ni+1] == 0 )
                  {
                     liste[n] = (int*)realloc(liste[n],(m+1)*sizeof(int));
                     liste[n][m] = liste[n][o]+ni+1;
                     checkNV[k][liste[n][m]] = 1;
                     m++;
                  }
                  if ( (liste[n][o]-ni-1) > 0  && (liste[n][o]%ni) != 0 && QR[k][liste[n][o]-ni-1] < seuil5 && checkNV[k][liste[n][o]-ni-1] == 0 )
                  {  
                     liste[n] = (int*)realloc(liste[n],(m+1)*sizeof(int));
                     liste[n][m] = liste[n][o]-ni-1;
                     checkNV[k][liste[n][m]] = 1;
                     m++;
                  }
                  if ( (liste[n][o]-ni) > 0 && QR[k][liste[n][o]-ni] < seuil5 && checkNV[k][liste[n][o]-ni] == 0 )
                  {
                     liste[n] = (int*)realloc(liste[n],(m+1)*sizeof(int));
                     liste[n][m] = liste[n][o]-ni;
                     checkNV[k][liste[n][m]] = 1;
                     m++;
                  }
                  if ( (liste[n][o]-ni+1) > 0 && ((liste[n][o]+1)%ni) != 0 && QR[k][liste[n][o]-ni+1] < seuil5 && checkNV[k][liste[n][o]-ni+1] == 0 ) 
                  {
                     liste[n] = (int*)realloc(liste[n],(m+1)*sizeof(int));
                     liste[n][m] = liste[n][o]-ni+1;
                     checkNV[k][liste[n][m]] = 1;
                     m++;
                  }
               }
                           
               size[n] = m;
               while ( ((GZ[k][i]  - GZ[0][i]) < 1200.) && m != 0 && k < (nk-1) )
               {
                  k++;
                  n++;

                  liste = (int**)realloc(liste,(n+1)*sizeof(int*));
                  liste[n] = (int*)(malloc)(1*sizeof(int));
                              
                  p = 0;
                  for ( o = 0; o < m; o++)
                  {
                     if ( QR[k][liste[n-1][o]] < seuil5 && checkNV[k][liste[n-1][o]] == 0 )
                     {
                        liste[n] = (int*)realloc(liste[n],(p+1)*sizeof(int));
                        liste[n][p] = liste[n-1][o];
                        checkNV[k][liste[n][p]] = 1;
                        p++;
                     }
                     if ( (liste[n-1][o]+1) < (ni*nj) && ((liste[n-1][o]+1)%ni) != 0 && QR[k][liste[n-1][o]+1] < seuil5 && checkNV[k][liste[n-1][o]+1] == 0 )
                     {
                        liste[n] = (int*)realloc(liste[n],(p+1)*sizeof(int));
                        liste[n][p] = liste[n-1][o]+1;
                        checkNV[k][liste[n][p]] = 1;
                        p++;
                     }
                     if ( (liste[n-1][o]-1) > 0 && (liste[n-1][o]%ni) != 0 && QR[k][liste[n-1][o]-1] < seuil5 && checkNV[k][liste[n-1][o]-1] == 0 )
                     {
                        liste[n] = (int*)realloc(liste[n],(p+1)*sizeof(int));
                        liste[n][p] = liste[n-1][o]-1;
                        checkNV[k][liste[n][p]] = 1;
                        p++;
                     }
                     if ( (liste[n-1][o]+ni-1) < (ni*nj) && (liste[n-1][o]%ni) != 0 && QR[k][liste[n-1][o]+ni-1] < seuil5 && checkNV[k][liste[n-1][o]+ni-1] == 0 )
                     {
                        liste[n] = (int*)realloc(liste[n],(p+1)*sizeof(int));
                        liste[n][p] = liste[n-1][o]+ni-1;
                        checkNV[k][liste[n][p]] = 1;
                        p++;
                     }
                     if ( (liste[n-1][o]+ni) < (ni*nj) && QR[k][liste[n-1][o]+ni] < seuil5 && checkNV[k][liste[n-1][o]+ni] == 0 )
                     {
                        liste[n] = (int*)realloc(liste[n],(p+1)*sizeof(int));
                        liste[n][p] = liste[n-1][o]+ni;
                        checkNV[k][liste[n][p]] = 1;
                        p++;
                     }
                     if ( (liste[n-1][o]+ni+1) < (ni*nj) && ((liste[n-1][o]+1)%ni) != 0 && QR[k][liste[n-1][o]+ni+1] < seuil5 && checkNV[k][liste[n-1][o]+ni+1] == 0 )
                     {
                        liste[n] = (int*)realloc(liste[n],(p+1)*sizeof(int));
                        liste[n][p] = liste[n-1][o]+ni+1;
                        checkNV[k][liste[n][p]] = 1;
                        p++;
                     }
                     if ( (liste[n-1][o]-ni-1) > 0 && (liste[n-1][o]%ni) != 0 && QR[k][liste[n-1][o]-ni-1] < seuil5 && checkNV[k][liste[n-1][o]-ni-1] == 0 )
                     {
                        liste[n] = (int*)realloc(liste[n],(p+1)*sizeof(int));
                        liste[n][p] = liste[n-1][o]-ni-1;
                        checkNV[k][liste[n][p]] = 1;
                        p++;
                     }
                     if ( (liste[n-1][o]-ni) > 0 && QR[k][liste[n-1][o]-ni] < seuil5 && checkNV[k][liste[n-1][o]-ni] == 0 )
                     {
                        liste[n] = (int*)realloc(liste[n],(p+1)*sizeof(int));
                        liste[n][p] = liste[n-1][o]-ni;
                        checkNV[k][liste[n][p]] = 1;
                        p++;
                     }
                     if ( (liste[n-1][o]-ni+1) > 0 && ((liste[n-1][o]+1)%ni) != 0 && QR[k][liste[n-1][o]-ni+1] < seuil5 && checkNV[k][liste[n-1][o]-ni+1] == 0 )
                     {
                        liste[n] = (int*)realloc(liste[n],(p+1)*sizeof(int));
                        liste[n][p] = liste[n-1][o]-ni+1;
                        checkNV[k][liste[n][p]] = 1;
                        p++;
                     }
                  }
                  for ( o = 0; o < p; o++)
                  {
                     if ( (liste[n][o]+1) < (ni*nj) && ((liste[n][o]+1)%ni) != 0 && QR[k][liste[n][o]+1] < seuil5 && checkNV[k][liste[n][o]+1] == 0 )
                     {
                        liste[n] = (int*)realloc(liste[n],(p+1)*sizeof(int));
                        liste[n][p] = liste[n][o]+1;
                        checkNV[k][liste[n][p]] = 1;
                        p++;
                     }
                     if ( (liste[n][o]-1) > 0 && (liste[n][o]%ni) != 0 && QR[k][liste[n][o]-1] < seuil5 && checkNV[k][liste[n][o]-1] == 0 )
                     {
                        liste[n] = (int*)realloc(liste[n],(p+1)*sizeof(int));
                        liste[n][p] = liste[n][o]-1;
                        checkNV[k][liste[n][p]] = 1;
                        p++;
                     }
                     if ( (liste[n][o]+ni-1) < (ni*nj) && (liste[n][o]%ni) != 0 && QR[k][liste[n][o]+ni-1] < seuil5 && checkNV[k][liste[n][o]+ni-1] == 0 )
                     {  
                        liste[n] = (int*)realloc(liste[n],(p+1)*sizeof(int));
                        liste[n][p] = liste[n][o]+ni-1;
                        checkNV[k][liste[n][p]] = 1;
                        p++;
                     }
                     if ( (liste[n][o]+ni) < (ni*nj) && QR[k][liste[n][o]+ni] < seuil5 && checkNV[k][liste[n][o]+ni] == 0 )
                     {
                        liste[n] = (int*)realloc(liste[n],(p+1)*sizeof(int));
                        liste[n][p] = liste[n][o]+ni;
                        checkNV[k][liste[n][p]] = 1;
                        p++;
                     }
                     if ( (liste[n][o]+ni+1) < (ni*nj) && ((liste[n][o]+1)%ni) != 0 && QR[k][liste[n][o]+ni+1] < seuil5 && checkNV[k][liste[n][o]+ni+1] == 0 )
                     {
                        liste[n] = (int*)realloc(liste[n],(p+1)*sizeof(int));
                        liste[n][p] = liste[n][o]+ni+1;
                        checkNV[k][liste[n][p]] = 1;
                        p++;
                     }
                     if ( (liste[n][o]-ni-1) > 0 && (liste[n][o]%ni) != 0 && QR[k][liste[n][o]-ni-1] < seuil5 && checkNV[k][liste[n][o]-ni-1] == 0 )
                     {  
                        liste[n] = (int*)realloc(liste[n],(p+1)*sizeof(int));
                        liste[n][p] = liste[n][o]-ni-1;
                        checkNV[k][liste[n][p]] = 1;
                        p++;
                     }
                     if ( (liste[n][o]-ni) > 0 && QR[k][liste[n][o]-ni] < seuil5 && checkNV[k][liste[n][o]-ni] == 0 )
                     {
                        liste[n] = (int*)realloc(liste[n],(p+1)*sizeof(int));
                        liste[n][p] = liste[n][o]-ni;
                        checkNV[k][liste[n][p]] = 1;
                        p++;
                     }
                     if ( (liste[n][o]-ni+1) > 0 && ((liste[n][o]+1)%ni) != 0 && QR[k][liste[n][o]-ni+1] < seuil5 && checkNV[k][liste[n][o]-ni+1] == 0 ) 
                     {
                        liste[n] = (int*)realloc(liste[n],(p+1)*sizeof(int));
                        liste[n][p] = liste[n][o]-ni+1;
                        checkNV[k][liste[n][p]] = 1;
                        p++;
                     }
                  }
                  size = (int*)realloc(size,(n+1)*sizeof(int));
                  size[n] = p;
                  m = p;
               }
               //n--;
               if ( (GZ[k][i]  - GZ[base][i]) > 200. )
               {
                  for ( o = 0; o < n; o++)
                  {
                     for ( q = 0; q < size[o]; q++)
                     {
                        Negative_Vorticity[liste[o][q]] = 1;
                     }
                  }
               }  
               for ( o = 0; o < (n+1); o++)
               {
                  free (liste[o]);
               }
               free(size);
               free(liste);
               size = NULL;
               liste = NULL;
            }
            //break;
         }
      }   
   } 

   
   
   //--------------------------------
   // Initialisation des variables
   //--------------------------------
   for ( i = 0; i < (ni*nj); i++)
   {
         checkNV_[i] = 0;
         checkPV_[i] = 0;
   }
   
   liste = (int**)malloc(1*sizeof(int*));
   n = 0;
   for ( i = 0; i < (ni*nj); i++)
   {
      CoupletPotentiel = 0;
      m = 0;
      if ( Positive_Vorticity[i] == 1 && checkPV_[i] == 0 )
      {
         liste[n] = (int*)(malloc)(1*sizeof(int));
         liste[n][m] = i;
         checkPV_[i] = 1;
         m++;
         
         for ( o = 0; o < m; o++)
         {
            if ( (liste[n][o]+1) < (ni*nj) && ((liste[n][o]+1)%ni) != 0 && Positive_Vorticity[liste[n][o]+1] == 1 && checkPV_[liste[n][o]+1] == 0 )
            {
               liste[n] = (int*)realloc(liste[n],(m+1)*sizeof(int));
               liste[n][m] = liste[n][o]+1;
               checkPV_[liste[n][o]+1] = 1;
               m++;
            }
            if ( (liste[n][o]-1) > 0 && (liste[n][o]%ni) != 0 && Positive_Vorticity[liste[n][o]-1] == 1 && checkPV_[liste[n][o]-1] == 0 )
            {
               liste[n] = (int*)realloc(liste[n],(m+1)*sizeof(int));
               liste[n][m] = liste[n][o]-1;
               checkPV_[liste[n][o]-1] = 1;
               m++;
            }
            if ( (liste[n][o]+ni-1) < (ni*nj) && (liste[n][o]%ni) != 0 && Positive_Vorticity[liste[n][o]+ni-1] == 1 && checkPV_[liste[n][o]+ni-1] == 0 )
            {
               liste[n] = (int*)realloc(liste[n],(m+1)*sizeof(int));
               liste[n][m] = liste[n][o]+ni-1;
               checkPV_[liste[n][o]+ni-1] = 1;
               m++;
            }
            if ( (liste[n][o]+ni) < (ni*nj) && Positive_Vorticity[liste[n][o]+ni] == 1 && checkPV_[liste[n][o]+ni] == 0 )
            {
               liste[n] = (int*)realloc(liste[n],(m+1)*sizeof(int));
               liste[n][m] = liste[n][o]+ni;
               checkPV_[liste[n][o]+ni] = 1;
               m++;
            }
            if ( (liste[n][o]+ni+1) < (ni*nj) && ((liste[n][o]+1)%ni) != 0 && Positive_Vorticity[liste[n][o]+ni+1] == 1 && checkPV_[liste[n][o]+ni+1] == 0 )
            {
               liste[n] = (int*)realloc(liste[n],(m+1)*sizeof(int));
               liste[n][m] = liste[n][o]+ni+1;
               checkPV_[liste[n][o]+ni+1] = 1;
               m++;
            }
            if ( (liste[n][o]-ni-1) > 0 && (liste[n][o]%ni) != 0 && Positive_Vorticity[liste[n][o]-ni-1] == 1 && checkPV_[liste[n][o]-ni-1] == 0 )
            {
               liste[n] = (int*)realloc(liste[n],(m+1)*sizeof(int));
               liste[n][m] = liste[n][o]-ni-1;
               checkPV_[liste[n][o]-ni-1] = 1;
               m++;
            }
            if ( (liste[n][o]-ni) > 0 && Positive_Vorticity[liste[n][o]-ni] == 1 && checkPV_[liste[n][o]-ni] == 0 )
            {
               liste[n] = (int*)realloc(liste[n],(m+1)*sizeof(int));
               liste[n][m] = liste[n][o]-ni;
               checkPV_[liste[n][o]-ni] = 1;
               m++;
            }
            if ( (liste[n][o]-ni+1) > 0 && ((liste[n][o]+1)%ni) != 0 && Positive_Vorticity[liste[n][o]-ni+1] == 1 && checkPV_[liste[n][o]-ni+1] == 0 )
            {
               liste[n] = (int*)realloc(liste[n],(m+1)*sizeof(int));
               liste[n][m] = liste[n][o]-ni+1;
               checkPV_[liste[n][o]-ni+1] = 1;
               m++;
            }
         }
         for ( o = 0; o < m; o++)
         {
            if ( Negative_Vorticity[liste[n][o]] == 1 )
               CoupletPotentiel = 1;
         }
         if ( CoupletPotentiel == 1 )
         {
            for ( o = 0; o < m; o++)
            {
               if ( (liste[n][o]+1) < (ni*nj) && ((liste[n][o]+1)%ni) != 0 && Negative_Vorticity[liste[n][o]+1] == 1 && checkNV_[liste[n][o]+1] == 0 )
               {
                  liste[n] = (int*)realloc(liste[n],(m+1)*sizeof(int));
                  liste[n][m] = liste[n][o]+1;
                  checkNV_[liste[n][o]+1] = 1;
                  m++;
               }
               if ( (liste[n][o]-1) > 0 && (liste[n][o]%ni) != 0 && Negative_Vorticity[liste[n][o]-1] == 1 && checkNV_[liste[n][o]-1] == 0 )
               {
                  liste[n] = (int*)realloc(liste[n],(m+1)*sizeof(int));
                  liste[n][m] = liste[n][o]-1;
                  checkNV_[liste[n][o]-1] = 1;
                  m++;
               }
               if ( (liste[n][o]+ni-1) < (ni*nj) && (liste[n][o]%ni) != 0 && Negative_Vorticity[liste[n][o]+ni-1] == 1 && checkNV_[liste[n][o]+ni-1] == 0 )
               {
                  liste[n] = (int*)realloc(liste[n],(m+1)*sizeof(int));
                  liste[n][m] = liste[n][o]+ni-1;
                  checkNV_[liste[n][o]+ni-1] = 1;
                  m++;
               }
               if ( (liste[n][o]+ni) < (ni*nj) && Negative_Vorticity[liste[n][o]+ni] == 1 && checkNV_[liste[n][o]+ni] == 0 )
               {
                  liste[n] = (int*)realloc(liste[n],(m+1)*sizeof(int));
                  liste[n][m] = liste[n][o]+ni;
                  checkNV_[liste[n][o]+ni] = 1;
                  m++;
               }
               if ( (liste[n][o]+ni+1) < (ni*nj) && ((liste[n][o]+1)%ni) != 0 && Negative_Vorticity[liste[n][o]+ni+1] == 1 && checkNV_[liste[n][o]+ni+1] == 0 )
               {
                  liste[n] = (int*)realloc(liste[n],(m+1)*sizeof(int));
                  liste[n][m] = liste[n][o]+ni+1;
                  checkNV_[liste[n][o]+ni+1] = 1;
                  m++;
               }
               if ( (liste[n][o]-ni-1) > 0 && (liste[n][o]%ni) != 0 && Negative_Vorticity[liste[n][o]-ni-1] == 1 && checkNV_[liste[n][o]-ni-1] == 0 )
               {
                  liste[n] = (int*)realloc(liste[n],(m+1)*sizeof(int));
                  liste[n][m] = liste[n][o]-ni-1;
                  checkNV_[liste[n][o]-ni-1] = 1;
                  m++;
               }
               if ( (liste[n][o]-ni) > 0 && Negative_Vorticity[liste[n][o]-ni] == 1 && checkNV_[liste[n][o]-ni] == 0 )
               {
                  liste[n] = (int*)realloc(liste[n],(m+1)*sizeof(int));
                  liste[n][m] = liste[n][o]-ni;
                  checkNV_[liste[n][o]-ni] = 1;
                  m++;
               }
               if ( (liste[n][o]-ni+1) > 0 && ((liste[n][o]+1)%ni) != 0 && Negative_Vorticity[liste[n][o]-ni+1] == 1 && checkNV_[liste[n][o]-ni+1] == 0 )
               {
                  liste[n] = (int*)realloc(liste[n],(m+1)*sizeof(int));
                  liste[n][m] = liste[n][o]-ni+1;
                  checkNV_[liste[n][o]-ni+1] = 1;
                  m++;
               }
            }

            hauteurMinNV = nk;
            hauteurMaxNV = 0;
            hauteurMinPV = nk;
            hauteurMaxPV = 0;
            for ( o = 0; o < m; o++)
            {
               for ( k = 0; k < (nk-1); k++)
               {
                  if ( checkPV[k][liste[n][o]] == 0 && checkPV[k+1][liste[n][o]] == 1 )
                  {
                     if ( (k+1) < hauteurMinPV )
                        hauteurMinPV = k + 1;
                  }
                  if ( checkNV[k][liste[n][o]] == 0 && checkNV[k+1][liste[n][o]] == 1 )
                  {
                     if ( (k+1) < hauteurMinNV )
                        hauteurMinNV = k + 1;
                  }
                  if ( checkPV[k][liste[n][o]] == 1 && checkPV[k+1][liste[n][o]] == 0 )
                  {
                     if ( k > hauteurMaxPV )
                        hauteurMaxPV = k;
                  }
                  if ( checkNV[k][liste[n][o]] == 1 && checkNV[k+1][liste[n][o]] == 0 )
                  {
                     if ( k > hauteurMaxNV )
                        hauteurMaxNV = k;
                  }
               } 
            }
            if ( (hauteurMinNV < hauteurMaxPV) && (hauteurMaxNV > hauteurMinPV) )
            {
               for ( o = 0; o < m; o++)
               {
                  Vorticity_couplet[liste[n][o]] = 1;
               }
            }
         }
         free(liste[n]);
      }
   }
   
   free(liste);
   liste = NULL;
   
   //---------------------
   // Liberer la memoire
   //---------------------
   for ( k = 0; k < nk; k++ )
   {
      free(checkNV[k]);
      free(checkPV[k]);
   }
   free(checkNV);
   free(checkPV);
   free(checkNV_);
   free(checkPV_);
   checkNV = NULL;
   checkPV = NULL;
   checkNV_ = NULL;
   checkPV_ = NULL;
   
   //--------------------------
   // Retourner les resultats
   //--------------------------
   *couplet_negatif = Negative_Vorticity;
   *couplet_positif = Positive_Vorticity;
   *couplet_tourbillon = Vorticity_couplet;
   
   return SUCCESS;
}


/*!
===============================================================================
\brief Detecter les mouvements verticaux forts



\param[in]  WW        : Tableau du mouvement vertical (Pa/s)
\param[in]  GZ        : Tableau de la hauteur geopotentielle (DAM)
\param[in]  ni        : Dimension en x
\param[in]  nj        : Dimension en y
\param[in]  nk        : Dimension en z
\param[out] Updraft   : Vecteur contenant les mouvements vertical forts detecter
                        ( 0 => Pas de mouvement vertical fort detecter
                          1 => Mouvement vertical fort detecter )


\date 17 main 2010
\version 0.0
\author Anna-Belle Filion

--------------------------------------------------------------------------------
 */
int DetecterUpdraft(float **WW, float **GZ, int ni, int nj, int nk,
		    float seuil[2], float **Updraft)
{
   int i, j, k, l, m, n, o, p;  
   int start, top, base;
   
   int *size  = (int*)NULL;
   int *size1 = (int*)NULL; 
   int *size2 = (int*)NULL; 
   int **liste  = (int**)NULL;
   int **liste1 = (int**)NULL;
   int **liste2 = (int**)NULL;
   int **check  = (int**)NULL;
   
   float seuil1 = 0, seuil4 = 0;
   
   float *Severe_Updraft = (float*)NULL;
   
   //-------------------------
   // Attribution des seuils 
   //-------------------------
   seuil1 = seuil[0];
   seuil4 = seuil[1];

   //---------------------
   // Allouer la memoire
   //---------------------
   check = (int**)malloc(nk*sizeof(int*));
   for ( k = 0; k < nk; k++ )
   { 
      check[k] = (int*)malloc(ni*nj*sizeof(int));
   }
   Severe_Updraft = (float*)malloc(ni*nj*sizeof(float));
   
   //----------------------------
   // Initialiser les variables
   //----------------------------
   for ( k = 0; k < nk; k++ )
   {
      for ( i = 0; i < (ni*nj); i++ )
      {
         check[k][i] = 0;
         if ( k == 0 )
            Severe_Updraft[i] = 0.;
      }
   }
   
   //-----------------------------------------
   // Trouver les mouvements verticaux forts
   //-----------------------------------------
   for ( i = 0; i < (ni*nj); i++)
   {
      l = 0;
      j = 0;

      for ( k = 0; k < nk; k++)
      {
         if ( (GZ[k][i] - GZ[0][i]) < 1200.)
         {
            m = 0;
            n = 0;
            if ( WW[k][i] < seuil1 && check[k][i] == 0 )
            {
               size = (int*)malloc(1*sizeof(int));
               liste = (int**)malloc(1*sizeof(int*));
               liste[n] = (int*)(malloc)(1*sizeof(int));
               liste[n][m] = i;
               check[k][i] = 1;
               start = k;
               m++;

               if ( (i+1) < (ni*nj) && ((i+1)%ni) != 0 && WW[k][i+1] < seuil4 && Severe_Updraft[i+1] == 0 && check[k][i+1] == 0 )
               {
                  liste[n] = (int*)realloc(liste[n],(m+1)*sizeof(int));
                  liste[n][m] = i+1;
                  check[k][i+1] = 1;
                  m++;
               }
               if ( (i-1) > 0 && (i%ni) != 0 && WW[k][i-1] < seuil4 && Severe_Updraft[i-1] == 0 && check[k][i-1] == 0 )
               {
                  liste[n] = (int*)realloc(liste[n],(m+1)*sizeof(int));
                  liste[n][m] = i-1;
                  check[k][i-1] = 1;
                  m++;  
               }
               if ( (i+ni-1) < (ni*nj) && (i%ni) != 0 && WW[k][i+ni-1] < seuil4 && Severe_Updraft[i+ni-1] == 0  && check[k][i+ni-1] == 0 )
               {
                  liste[n] = (int*)realloc(liste[n],(m+1)*sizeof(int));
                  liste[n][m] = i+ni-1;
                  check[k][i+ni-1] = 1;
                  m++;
               }
               if ( (i+ni) < (ni*nj) && WW[k][i+ni] < seuil4 && Severe_Updraft[i+ni] == 0 && check[k][i+ni] == 0 )
               {
                  liste[n] = (int*)realloc(liste[n],(m+1)*sizeof(int));
                  liste[n][m] = i+ni;
                  check[k][i+ni] = 1;
                  m++;      
               }
               if ( (i+ni+1) < (ni*nj) && ((i+1)%ni) != 0 && WW[k][i+ni+1] < seuil4 && Severe_Updraft[i+ni+1] == 0 && check[k][i+ni+1] == 0 )
               {
                  liste[n] = (int*)realloc(liste[n],(m+1)*sizeof(int));
                  liste[n][m] = i+ni+1;
                  check[k][i+ni+1] = 1;
                  m++;     
               }
               if ( (i-ni-1) > 0 && (i%ni) != 0 && WW[k][i-ni-1] < seuil4 && Severe_Updraft[i-ni-1] == 0 && check[k][i-ni-1] == 0 )
               {
                  liste[n] = (int*)realloc(liste[n],(m+1)*sizeof(int));
                  liste[n][m] = i-ni-1;
                  check[k][i-ni-1] = 1;
                  m++;      
               }
               if ( (i-ni) > 0 && WW[k][i-ni] < seuil4 && Severe_Updraft[i-ni] == 0 && check[k][i-ni] == 0 )
               {
                  liste[n] = (int*)realloc(liste[n],(m+1)*sizeof(int));
                  liste[n][m] = i-ni;
                  check[k][i-ni] = 1;
                  m++;     
               }
               if ( (i-ni+1) > 0 && ((i+1)%ni) != 0 && WW[k][i-ni+1] < seuil4 && Severe_Updraft[i-ni+1] == 0 && check[k][i-ni+1] == 0 ) 
               {
                  liste[n] = (int*)realloc(liste[n],(m+1)*sizeof(int));
                  liste[n][m] = i-ni+1;
                  check[k][i-ni+1] = 1;
                  m++;     
               }

               for ( o = 1; o < m; o++)
               {
                  if ( (liste[n][o]+1) < (ni*nj) && ((liste[n][o]+1)%ni) != 0 && WW[k][liste[n][o]+1] < seuil4 && Severe_Updraft[liste[n][o]+1] == 0 && check[k][liste[n][o]+1] == 0 )
                  {
                     liste[n] = (int*)realloc(liste[n],(m+1)*sizeof(int));
                     liste[n][m] = liste[n][o]+1;
                     check[k][liste[n][m]] = 1;
                     m++;
                  }
                  if ( (liste[n][o]-1) > 0 && (liste[n][o]%ni) != 0 && WW[k][liste[n][o]-1] < seuil4 && Severe_Updraft[liste[n][o]-1] == 0 && check[k][liste[n][o]-1] == 0 )
                  {
                     liste[n] = (int*)realloc(liste[n],(m+1)*sizeof(int));
                     liste[n][m] = liste[n][o]-1;
                     check[k][liste[n][m]] = 1;
                     m++;
                  }
                  if ( (liste[n][o]+ni-1) < (ni*nj) && (liste[n][o]%ni) != 0 && WW[k][liste[n][o]+ni-1] < seuil4 && Severe_Updraft[liste[n][o]+ni-1] == 0 && check[k][liste[n][o]+ni-1] == 0 )
                  {
                     liste[n] = (int*)realloc(liste[n],(m+1)*sizeof(int));
                     liste[n][m] = liste[n][o]+ni-1;
                     check[k][liste[n][m]] = 1;
                     m++;
                  }
                  if ( (liste[n][o]+ni) < (ni*nj) && WW[k][liste[n][o]+ni] < seuil4 && Severe_Updraft[liste[n][o]+ni] == 0 && check[k][liste[n][o]+ni] == 0 )
                  {
                     liste[n] = (int*)realloc(liste[n],(m+1)*sizeof(int));
                     liste[n][m] = liste[n][o]+ni;
                     check[k][liste[n][m]] = 1;
                     m++;
                  }
                  if ( (liste[n][o]+ni+1) < (ni*nj) && ((liste[n][o]+1)%ni) != 0 && WW[k][liste[n][o]+ni+1] < seuil4 && Severe_Updraft[liste[n][o]+ni+1] == 0 && check[k][liste[n][o]+ni+1] == 0 )
                  {
                     liste[n] = (int*)realloc(liste[n],(m+1)*sizeof(int));
                     liste[n][m] = liste[n][o]+ni+1;
                     check[k][liste[n][m]] = 1;
                     m++;
                  }
                  if ( (liste[n][o]-ni-1) > 0 && (liste[n][o]%ni) != 0 && WW[k][liste[n][o]-ni-1] < seuil4 && Severe_Updraft[liste[n][o]-ni-1] == 0 && check[k][liste[n][o]-ni-1] == 0 )
                  {  
                     liste[n] = (int*)realloc(liste[n],(m+1)*sizeof(int));
                     liste[n][m] = liste[n][o]-ni-1;
                     check[k][liste[n][m]] = 1;
                     m++;
                  }
                  if ( (liste[n][o]-ni) > 0 && WW[k][liste[n][o]-ni] < seuil4 && Severe_Updraft[liste[n][o]-ni] == 0 && check[k][liste[n][o]-ni] == 0 )
                  {
                     liste[n] = (int*)realloc(liste[n],(m+1)*sizeof(int));
                     liste[n][m] = liste[n][o]-ni;
                     check[k][liste[n][m]] = 1;
                     m++;
                  }
                  if ( (liste[n][o]-ni+1) > 0 && ((liste[n][o]+1)%ni) != 0 && WW[k][liste[n][o]-ni+1] < seuil4 && Severe_Updraft[liste[n][o]-ni+1] == 0 && check[k][liste[n][o]-ni+1] == 0 ) 
                  {
                     liste[n] = (int*)realloc(liste[n],(m+1)*sizeof(int));
                     liste[n][m] = liste[n][o]-ni+1;
                     check[k][liste[n][m]] = 1;
                     m++;
                  }
               }   
               size[n] = m;    
               while ( (GZ[k][i] - GZ[0][i]) < 1200. && m != 0 && k < (nk-1) )
               {
                  n++;
                  k++;
                  liste = (int**)realloc(liste,(n+1)*sizeof(int*));
                  liste[n] = (int*)(malloc)(1*sizeof(int));
                         
                  p = 0; 
                  for ( o = 0; o < m; o++)
                  {
                     if ( WW[k][liste[n-1][o]] < seuil4 && check[k][liste[n-1][o]] == 0 )
                     {
                        liste[n] = (int*)realloc(liste[n],(p+1)*sizeof(int));
                        liste[n][p] = liste[n-1][o];
                        check[k][liste[n][p]] = 1;
                        p++;
                     }
                     if ( (liste[n-1][o]+1) < (ni*nj) && ((liste[n-1][o]+1)%ni) != 0 && WW[k][liste[n-1][o]+1] < seuil4 && check[k][liste[n-1][o]+1] == 0 )
                     {
                        liste[n] = (int*)realloc(liste[n],(p+1)*sizeof(int));
                        liste[n][p] = liste[n-1][o]+1;
                        check[k][liste[n][p]] = 1;
                        p++;
                     }
                     if ( (liste[n-1][o]-1) > 0 && (liste[n-1][o]%ni) != 0 && WW[k][liste[n-1][o]-1] < seuil4 && check[k][liste[n-1][o]-1] == 0 )
                     {
                        liste[n] = (int*)realloc(liste[n],(p+1)*sizeof(int));
                        liste[n][p] = liste[n-1][o]-1;
                        check[k][liste[n][p]] = 1;
                        p++;
                     }
                     if ( (liste[n-1][o]+ni-1) < (ni*nj) && (liste[n-1][o]%ni) != 0 && WW[k][liste[n-1][o]+ni-1] < seuil4 && check[k][liste[n-1][o]+ni-1] == 0 )
                     {
                        liste[n] = (int*)realloc(liste[n],(p+1)*sizeof(int));
                        liste[n][p] = liste[n-1][o]+ni-1;
                        check[k][liste[n][p]] = 1;
                        p++;
                     }
                     if ( (liste[n-1][o]+ni) < (ni*nj) && WW[k][liste[n-1][o]+ni] < seuil4 && check[k][liste[n-1][o]+ni] == 0 )
                     {
                        liste[n] = (int*)realloc(liste[n],(p+1)*sizeof(int));
                        liste[n][p] = liste[n-1][o]+ni;
                        check[k][liste[n][p]] = 1;
                        p++;
                     }
                     if ( (liste[n-1][o]+ni+1) < (ni*nj) && ((liste[n-1][o]+1)%ni) != 0 && WW[k][liste[n-1][o]+ni+1] < seuil4 && check[k][liste[n-1][o]+ni+1] == 0 )
                     {
                        liste[n] = (int*)realloc(liste[n],(p+1)*sizeof(int));
                        liste[n][p] = liste[n-1][o]+ni+1;
                        check[k][liste[n][p]] = 1;
                        p++;
                     }
                     if ( (liste[n-1][o]-ni-1) > 0 && (liste[n-1][o]%ni) != 0 && WW[k][liste[n-1][o]-ni-1] < seuil4 && check[k][liste[n-1][o]-ni-1] == 0 )
                     {
                        liste[n] = (int*)realloc(liste[n],(p+1)*sizeof(int));
                        liste[n][p] = liste[n-1][o]-ni-1;
                        check[k][liste[n][p]] = 1;
                        p++;
                     }
                     if ( (liste[n-1][o]-ni) > 0 && WW[k][liste[n-1][o]-ni] < seuil4 && check[k][liste[n-1][o]-ni] == 0 )
                     {
                        liste[n] = (int*)realloc(liste[n],(p+1)*sizeof(int));
                        liste[n][p] = liste[n-1][o]-ni;
                        check[k][liste[n][p]] = 1;
                        p++;
                     }
                     if ( (liste[n-1][o]-ni+1) > 0 && ((liste[n-1][o]+1)%ni) != 0 && WW[k][liste[n-1][o]-ni+1] < seuil4 && check[k][liste[n-1][o]-ni+1] == 0 )
                     {
                        liste[n] = (int*)realloc(liste[n],(p+1)*sizeof(int));
                        liste[n][p] = liste[n-1][o]-ni+1;
                        check[k][liste[n][p]] = 1;
                        p++;
                     }
                  }

                  for ( o = 0; o < p; o++)
                  {
                     if ( (liste[n][o]+1) < (ni*nj) && ((liste[n][o]+1)%ni) != 0 && WW[k][liste[n][o]+1] < seuil4 && check[k][liste[n][o]+1] == 0 )
                     {
                        liste[n] = (int*)realloc(liste[n],(p+1)*sizeof(int));
                        liste[n][p] = liste[n][o]+1;
                        check[k][liste[n][p]] = 1;
                        p++;
                     }
                     if ( (liste[n][o]-1) > 0 && (liste[n][o]%ni) != 0 && WW[k][liste[n][o]-1] < seuil4 && check[k][liste[n][o]-1] == 0 )
                     {
                        liste[n] = (int*)realloc(liste[n],(p+1)*sizeof(int));
                        liste[n][p] = liste[n][o]-1;
                        check[k][liste[n][p]] = 1;
                        p++;
                     }
                     if ( (liste[n][o]+ni-1) < (ni*nj) && (liste[n][o]%ni) != 0 && WW[k][liste[n][o]+ni-1] < seuil4 && check[k][liste[n][o]+ni-1] == 0 )
                     {
                        liste[n] = (int*)realloc(liste[n],(p+1)*sizeof(int));
                        liste[n][p] = liste[n][o]+ni-1;
                        check[k][liste[n][p]] = 1;
                        p++;
                     }
                     if ( (liste[n][o]+ni) < (ni*nj) && WW[k][liste[n][o]+ni] < seuil4 && check[k][liste[n][o]+ni] == 0 )
                     {
                        liste[n] = (int*)realloc(liste[n],(p+1)*sizeof(int));
                        liste[n][p] = liste[n][o]+ni;
                        check[k][liste[n][p]] = 1;
                        p++;
                     }
                     if ( (liste[n][o]+ni+1) < (ni*nj) && ((liste[n][o]+1)%ni) != 0 && WW[k][liste[n][o]+ni+1] < seuil4 && check[k][liste[n][o]+ni+1] == 0 )
                     {
                        liste[n] = (int*)realloc(liste[n],(p+1)*sizeof(int));
                        liste[n][p] = liste[n][o]+ni+1;
                        check[k][liste[n][p]] = 1;
                        p++;
                     }
                     if ( (liste[n][o]-ni-1) > 0 && (liste[n][o]%ni) != 0 && WW[k][liste[n][o]-ni-1] < seuil4 && check[k][liste[n][o]-ni-1] == 0 )
                     {
                        liste[n] = (int*)realloc(liste[n],(p+1)*sizeof(int));
                        liste[n][p] = liste[n][o]-ni-1;
                        check[k][liste[n][p]] = 1;
                        p++;
                     }
                     if ( (liste[n][o]-ni) > 0 && WW[k][liste[n][o]-ni] < seuil4 && check[k][liste[n][o]-ni] == 0 )
                     {
                        liste[n] = (int*)realloc(liste[n],(p+1)*sizeof(int));
                        liste[n][p] = liste[n][o]-ni;
                        check[k][liste[n][p]] = 1;
                        p++;
                     }
                     if ( (liste[n][o]-ni+1) > 0 && ((liste[n][o]+1)%ni) != 0 && WW[k][liste[n][o]-ni+1] < seuil4 && check[k][liste[n][o]-ni+1] == 0 )
                     {
                        liste[n] = (int*)realloc(liste[n],(p+1)*sizeof(int));
                        liste[n][p] = liste[n][o]-ni+1;
                        check[k][liste[n][p]] = 1;
                        p++;
                     }
                  }

                  size = (int*)realloc(size,(n+1)*sizeof(int));
                  size[n] = p;
                  m = p;
               }

               top = k;

               size1 = (int*)malloc(1*sizeof(int));
               liste1 = (int**)malloc(1*sizeof(int*));
               liste1[0] = (int*)malloc(size[0]*sizeof(int));
         
               size1[0] = size[0];
               for ( o = 0; o < size1[0]; o++)
               {
                  liste1[0][o] = liste[0][o];
               }
               
               n = 0;
               k = start;
               m = size1[0];
               while ( m != 0 && k > 0 )
               {
                  n++;
                  k--;
                  liste1 = (int**)realloc(liste1,(n+1)*sizeof(int*));
                  liste1[n] = (int*)(malloc)(1*sizeof(int));

                  p = 0;
                  for ( o = 0; o < m; o++)
                  {
                     if ( WW[k][liste1[n-1][o]] < seuil4 && check[k][liste1[n-1][o]] == 0 )
                     {
                        liste1[n] = (int*)realloc(liste1[n],(p+1)*sizeof(int));
                        liste1[n][p] = liste1[n-1][o];
                        check[k][liste1[n][p]] = 1;
                        p++;
                     }
                     if ( (liste1[n-1][o]+1) < (ni*nj) && ((liste1[n-1][o]+1)%ni) != 0 && WW[k][liste1[n-1][o]+1] < seuil4 && check[k][liste1[n-1][o]+1] == 0 )
                     {
                        liste1[n] = (int*)realloc(liste1[n],(p+1)*sizeof(int));
                        liste1[n][p] = liste1[n-1][o]+1;
                        check[k][liste1[n][p]] = 1;
                        p++;
                     }        
                     if ( (liste1[n-1][o]-1) > 0 && (liste1[n-1][o]%ni) != 0 && WW[k][liste1[n-1][o]-1] < seuil4 && check[k][liste1[n-1][o]-1] == 0 )
                     {
                        liste1[n] = (int*)realloc(liste1[n],(p+1)*sizeof(int));
                        liste1[n][p] = liste1[n-1][o]-1;
                        check[k][liste1[n][p]] = 1;
                        p++;
                     }
                     if ( (liste1[n-1][o]+ni-1) < (ni*nj) && (liste1[n-1][o]%ni) != 0 && WW[k][liste1[n-1][o]+ni-1] < seuil4 && check[k][liste1[n-1][o]+ni-1] == 0 )
                     {
                        liste1[n] = (int*)realloc(liste1[n],(p+1)*sizeof(int));
                        liste1[n][p] = liste1[n-1][o]+ni-1;
                        check[k][liste1[n][p]] = 1;
                        p++;
                     }
                     if ( (liste1[n-1][o]+ni) < (ni*nj) && WW[k][liste1[n-1][o]+ni] < seuil4 && check[k][liste1[n-1][o]+ni] == 0 )
                     {
                        liste1[n] = (int*)realloc(liste1[n],(p+1)*sizeof(int));
                        liste1[n][p] = liste1[n-1][o]+ni;
                        check[k][liste1[n][p]] = 1;
                        p++;
                     }
                     if ( (liste1[n-1][o]+ni+1) < (ni*nj) && ((liste1[n-1][o]+1)%ni) != 0 && WW[k][liste1[n-1][o]+ni+1] < seuil4 && check[k][liste1[n-1][o]+ni+1] == 0 )
                     {
                        liste1[n] = (int*)realloc(liste1[n],(p+1)*sizeof(int));
                        liste1[n][p] = liste1[n-1][o]+ni+1;
                        check[k][liste1[n][p]] = 1;
                        p++;
                     }
                     if ( (liste1[n-1][o]-ni-1) > 0 && (liste1[n-1][o]%ni) != 0 && WW[k][liste1[n-1][o]-ni-1] < seuil4 && check[k][liste1[n-1][o]-ni-1] == 0 )
                     {        
                        liste1[n] = (int*)realloc(liste1[n],(p+1)*sizeof(int));
                        liste1[n][p] = liste1[n-1][o]-ni-1;
                        check[k][liste1[n][p]] = 1;
                        p++;
                     }        
                     if ( (liste1[n-1][o]-ni) > 0 && WW[k][liste1[n-1][o]-ni] < seuil4 && check[k][liste1[n-1][o]-ni] == 0 )
                     {
                        liste1[n] = (int*)realloc(liste1[n],(p+1)*sizeof(int));
                        liste1[n][p] = liste1[n-1][o]-ni;
                        check[k][liste1[n][p]] = 1;
                        p++;
                     }
                     if ( (liste1[n-1][o]-ni+1) > 0 && ((liste1[n-1][o]+1)%ni) != 0 && WW[k][liste1[n-1][o]-ni+1] < seuil4 && check[k][liste1[n-1][o]-ni+1] == 0 )
                     {
                        liste1[n] = (int*)realloc(liste1[n],(p+1)*sizeof(int));
                        liste1[n][p] = liste1[n-1][o]-ni+1;
                        check[k][liste1[n][p]] = 1;
                        p++;
                     }
                  }
                  for ( o = 0; o < p; o++)
                  {
                     if ( (liste1[n][o]+1) < (ni*nj) && ((liste1[n][o]+1)%ni) != 0 && WW[k][liste1[n][o]+1] < seuil4 && check[k][liste1[n][o]+1] == 0 )
                     {
                        liste1[n] = (int*)realloc(liste1[n],(p+1)*sizeof(int));
                        liste1[n][p] = liste1[n][o]+1;
                        check[k][liste1[n][p]] = 1;
                        p++;
                     }
                     if ( (liste1[n][o]-1) > 0 && (liste1[n][o]%ni) != 0 && WW[k][liste1[n][o]-1] < seuil4 && check[k][liste1[n][o]-1] == 0 )
                     {
                        liste1[n] = (int*)realloc(liste1[n],(p+1)*sizeof(int));
                        liste1[n][p] = liste1[n][o]-1;
                        check[k][liste1[n][p]] = 1;
                        p++;
                     }
                     if ( (liste1[n][o]+ni-1) < (ni*nj) && (liste1[n][o]%ni) != 0 && WW[k][liste1[n][o]+ni-1] < seuil4 && check[k][liste1[n][o]+ni-1] == 0 )
                     {  
                        liste1[n] = (int*)realloc(liste1[n],(p+1)*sizeof(int));
                        liste1[n][p] = liste1[n][o]+ni-1;
                        check[k][liste1[n][p]] = 1;
                        p++;
                     }
                     if ( (liste1[n][o]+ni) < (ni*nj) && WW[k][liste1[n][o]+ni] < seuil4 && check[k][liste1[n][o]+ni] == 0 )
                     {
                        liste1[n] = (int*)realloc(liste1[n],(p+1)*sizeof(int));
                        liste1[n][p] = liste1[n][o]+ni;
                        check[k][liste1[n][p]] = 1;
                        p++;
                     }
                     if ( (liste1[n][o]+ni+1) < (ni*nj) && ((liste1[n][o]+1)%ni) != 0 && WW[k][liste1[n][o]+ni+1] < seuil4 && check[k][liste1[n][o]+ni+1] == 0 )
                     {
                        liste1[n] = (int*)realloc(liste1[n],(p+1)*sizeof(int));
                        liste1[n][p] = liste1[n][o]+ni+1;
                        check[k][liste1[n][p]] = 1;
                        p++;
                     }
                     if ( (liste1[n][o]-ni-1) > 0 && (liste1[n][o]%ni) != 0 && WW[k][liste1[n][o]-ni-1] < seuil4 && check[k][liste1[n][o]-ni-1] == 0 )
                     {         
                        liste1[n] = (int*)realloc(liste1[n],(p+1)*sizeof(int));
                        liste1[n][p] = liste1[n][o]-ni-1;
                        check[k][liste1[n][p]] = 1;
                        p++;
                     }
                     if ( (liste1[n][o]-ni) > 0 && WW[k][liste1[n][o]-ni] < seuil4 && check[k][liste1[n][o]-ni] == 0 )
                     {
                        liste1[n] = (int*)realloc(liste1[n],(p+1)*sizeof(int));
                        liste1[n][p] = liste1[n][o]-ni;
                        check[k][liste1[n][p]] = 1;
                        p++;
                     }
                     if ( (liste1[n][o]-ni+1) > 0 && ((liste1[n][o]+1)%ni) != 0 && WW[k][liste1[n][o]-ni+1] < seuil4 && check[k][liste1[n][o]-ni+1] == 0 ) 
                     {
                        liste1[n] = (int*)realloc(liste1[n],(p+1)*sizeof(int));
                        liste1[n][p] = liste1[n][o]-ni+1;
                        check[k][liste1[n][p]] = 1;
                        p++;
                     }
                  }
                  size1 = (int*)realloc(size1,(n+1)*sizeof(int));
                  size1[n] = p;
                  m = p;
               }
               
               base = k;

               size2 = (int*)malloc((top-base-1)*sizeof(int));
               liste2 = (int**)malloc((top-base-1)*sizeof(int*));
               
               for ( o = (start-base-1); o > 0; o--)
               {
                  liste2[start-base-1-o] = (int*)malloc(size1[o]*sizeof(int));
                  size2[start-base-1-o] = size1[o];
                  for ( p = 0; p < size1[o]; p++ )
                  {
                     liste2[start-base-1-o][p] = liste1[o][p];
                  }
               }
               for ( o = 0; o < (top-start); o++)
               {
                  liste2[start-base-1+o] = (int*)malloc(size[o]*sizeof(int));
                  size2[start-base-1+o] = size[o];
                  for ( p = 0; p < size[o]; p++ )
                  {
                     liste2[start-base-1+o][p] = liste[o][p];
                  }
               }
        
               for ( o = 0; o < (top-start+1); o++)
               {
                  free(liste[o]);
               }
               for ( o = 0; o < (start-base+1); o++)
               {
                  free(liste1[o]);
               }
	       
               free(size);
               free(size1);
               free(liste);
               free(liste1);
               size = NULL;
               size1 = NULL;
               liste = NULL;
               liste1 = NULL;                 
  
               if ( (GZ[top-1][i] - GZ[base+1][i]) > 200. )
               {
                  for ( o = 0; o < (top-base-1); o++ )
                  {
                     for ( p = 0; p < size2[o]; p++ )
                     { 
                        Severe_Updraft[liste2[o][p]] = 1;
                     }
                  }
               }  
               for ( o = 0; o < (top-base-1); o++)
               {
                  free(liste2[o]);
               }
               free(size2);
               free(liste2);
               size2 = NULL;
               liste2 = NULL;   
               
               break;
            }
         }
         else
            break;
      }
   }
      
   //---------------------
   // Liberer la memoire
   //---------------------
   for ( k = 0; k < nk; k++ )
   {
      free(check[k]);
   }
   free(check);
   check = NULL;
   
   //--------------------------
   // Retourner les resultats
   //--------------------------
   *Updraft = Severe_Updraft;
   
   return SUCCESS;
}

/*!
===============================================================================
\brief Permet de reperer les zones ou il y a du surplomb. 


\param[in]  ZET             : donnees ZET          (dBz)
\param[in]  ZEC             : MAXR                 (dBz)
\param[in]  indiceMaxR      : hauteur des MAXR     (m)
\param[in]  GZ              : hauteur              (m)
\param[in]  ni              : dimension en x
\param[in]  nj              : dimension en y
\param[in]  nk              : dimension en z
\param[in]  parametre       : structure qui contient les parametres du programme.
\param[out] ozoneSurplomb   : grille contenant nombre de surplombs
\param[out] oBWER           : vecteur contenant les points ou une voute d'echo faible
                              a ete trouver


\date 10 mars 2010
\version 0.1
\author Anna-Belle Filion
-------------------------------------------------------------------------------
 */
int repererSurplomb(float **ZET, float **GZ, float *MaxR, int ni, int nj, 
		    int nk, float res, float **ozoneSurplomb, float **oBWER)
{
   int i, k, l, m, n ,o, p, k_max, i_max = 0;
   int base, top, start;
   int z = 0;
   
   int **check   = (int**)NULL;
   int **liste   = (int**)NULL;
   int **liste1 = (int**)NULL;
   int **liste2 = (int**)NULL;
   int *size     = (int*)NULL;
   int *size1   = (int*)NULL;
   int *size2   = (int*)NULL;
   
   float seuil = 45., seuil2 = 30.;
   //float seuil = 45., seuil2 = 30.;
   
   float *BWER = (float*)NULL;
   float *zoneSurplomb = (float*)NULL;         
   
   int Surplomb;
   
   //-----------------------
   // Allocation de memoire
   //-----------------------
   check = (int**)malloc(nk*sizeof(int*));
   for ( k = 0; k < nk; k++ )
   {
      check[k] = (int*)malloc(ni*nj*sizeof(int));
   }
   BWER  = (float*)malloc(ni*nj*sizeof(float));
   zoneSurplomb = (float*)malloc(ni*nj*sizeof(float));
   
   //------------------------------
   // Initialiser les variables
   //------------------------------
   for ( k = 0; k < nk; k++ )
   {
      for ( i = 0; i < (ni*nj); i++)
      {
         check[k][i] = 0;
         if ( k == 0 )
         {
            BWER[i] = 0.;
            zoneSurplomb[i] = 0;
         }
      }
   }
   
   //-----------------------------------------------------
   // Trouver les surplombs et les voutes d'echo faible
   //-----------------------------------------------------
   for ( i = 0; i < (ni*nj); i++)
   {
      if ( MaxR[i] >= seuil )
      {

         k = 0;
         while ( ZET[k][i] < seuil2 )
         {
            k++;
         }
         start = k;
         base = k;
	 top = k;
	 
         m = 0; 
         n = 0;
         l = 0;
         
         size   = (int*)malloc(1*sizeof(int));
         liste   = (int**)malloc(1*sizeof(int*));
         
         liste[n] = (int*)(malloc)(1*sizeof(int));
         liste[n][m] = i;
         check[k][i] = 1;
         m++;
         
         if ( (i+1) < (ni*nj) && ((i+1)%ni) != 0 && ZET[k][i+1] > seuil2 && check[k][i+1] == 0 )
         {
            liste[n] = (int*)realloc(liste[n],(m+1)*sizeof(int));
            liste[n][m] = i+1;
            check[k][i+1] = 1;
            m++;
	 }
         if ( (i-1) > 0 && (i%ni) != 0 && ZET[k][i-1] > seuil2 && check[k][i-1] == 0 )
         {
            liste[n] = (int*)realloc(liste[n],(m+1)*sizeof(int));
            liste[n][m] = i-1;
            check[k][i-1] = 1;
            m++;   
         }
         if ( (i+ni-1) < (ni*nj) && (i%ni) != 0 && ZET[k][i+ni-1] > seuil2 && check[k][i+ni-1] == 0 )
         {
            liste[n] = (int*)realloc(liste[n],(m+1)*sizeof(int));
            liste[n][m] = i+ni-1;
            check[k][i+ni-1] = 1;
            m++;
         }
         if ( (i+ni) < (ni*nj) && ZET[k][i+ni] > seuil2 && check[k][i+ni] == 0 )
         {
            liste[n] = (int*)realloc(liste[n],(m+1)*sizeof(int));
            liste[n][m] = i+ni;
            check[k][i+ni] = 1;
            m++;   
         }
         if ( (i+ni+1) < (ni*nj) && ((i+1)%ni) != 0 && ZET[k][i+ni+1] > seuil2 && check[k][i+ni+1] == 0 )
         {
            liste[n] = (int*)realloc(liste[n],(m+1)*sizeof(int));
            liste[n][m] = i+ni+1;
            check[k][i+ni+1] = 1;
            m++; 
         }
         if ( (i-ni-1) > 0 && (i%ni) != 0 && ZET[k][i-ni-1] > seuil2 && check[k][i-ni-1] == 0 )
         {
            liste[n] = (int*)realloc(liste[n],(m+1)*sizeof(int));
            liste[n][m] = i-ni-1;
            check[k][i-ni-1] = 1;
            m++;   
         }
         if ( (i-ni) > 0 && ZET[k][i-ni] > seuil2 && check[k][i-ni] == 0 )
         {
            liste[n] = (int*)realloc(liste[n],(m+1)*sizeof(int));
            liste[n][m] = i-ni;
            check[k][i-ni] = 1;
            m++;  
         }
         if ( (i-ni+1) > 0 && ((i+1)%ni) != 0 && ZET[k][i-ni+1] > seuil2 && check[k][i-ni+1] == 0 ) 
         {
            liste[n] = (int*)realloc(liste[n],(m+1)*sizeof(int));
            liste[n][m] = i-ni+1;
            check[k][i-ni+1] = 1;
            m++;   
         }
         
         for ( o = 1; o < m; o++)
         {
            if ( (liste[n][o]+1) < (ni*nj) && ((liste[n][o]+1)%ni) != 0 && ZET[k][liste[n][o]+1] > seuil2 && check[k][liste[n][o]+1] == 0 )
            {
               liste[n] = (int*)realloc(liste[n],(m+1)*sizeof(int));
               liste[n][m] = liste[n][o]+1;
               check[k][liste[n][m]] = 1;
               m++;
            }
            if ( (liste[n][o]-1) > 0 && (liste[n][o]%ni) != 0 && ZET[k][liste[n][o]-1] > seuil2 && check[k][liste[n][o]-1] == 0 )
            {
               liste[n] = (int*)realloc(liste[n],(m+1)*sizeof(int));
               liste[n][m] = liste[n][o]-1;
               check[k][liste[n][m]] = 1;
               m++;
            }
            if ( (liste[n][o]+ni-1) < (ni*nj) && (liste[n][o]%ni) != 0 && ZET[k][liste[n][o]+ni-1] > seuil2 && check[k][liste[n][o]+ni-1] == 0 )
            {
               liste[n] = (int*)realloc(liste[n],(m+1)*sizeof(int));
               liste[n][m] = liste[n][o]+ni-1;
               check[k][liste[n][m]] = 1;
               m++;
            }
            if ( (liste[n][o]+ni) < (ni*nj) && ZET[k][liste[n][o]+ni] > seuil2 && check[k][liste[n][o]+ni] == 0 )
            {
               liste[n] = (int*)realloc(liste[n],(m+1)*sizeof(int));
               liste[n][m] = liste[n][o]+ni;
               check[k][liste[n][m]] = 1;
               m++;
            }
            if ( (liste[n][o]+ni+1) < (ni*nj) && ((liste[n][o]+1)%ni) != 0 && ZET[k][liste[n][o]+ni+1] > seuil2 && check[k][liste[n][o]+ni+1] == 0 )
            {
               liste[n] = (int*)realloc(liste[n],(m+1)*sizeof(int));
               liste[n][m] = liste[n][o]+ni+1;
               check[k][liste[n][m]] = 1;
               m++;
            }
            if ( (liste[n][o]-ni-1) > 0 && (liste[n][o]%ni) != 0 && ZET[k][liste[n][o]-ni-1] > seuil2 && check[k][liste[n][o]-ni-1] == 0 )
            {  
               liste[n] = (int*)realloc(liste[n],(m+1)*sizeof(int));
               liste[n][m] = liste[n][o]-ni-1;
               check[k][liste[n][m]] = 1;
               m++;
            }
            if ( (liste[n][o]-ni) > 0 && ZET[k][liste[n][o]-ni] > seuil2 && check[k][liste[n][o]-ni] == 0 )
            {
               liste[n] = (int*)realloc(liste[n],(m+1)*sizeof(int));
               liste[n][m] = liste[n][o]-ni;
               check[k][liste[n][m]] = 1;
               m++;
            }
            if ( (liste[n][o]-ni+1) > 0 && ((liste[n][o]+1)%ni) != 0 && ZET[k][liste[n][o]-ni+1] > seuil2 && check[k][liste[n][o]-ni+1] == 0 ) 
            {
               liste[n] = (int*)realloc(liste[n],(m+1)*sizeof(int));
               liste[n][m] = liste[n][o]-ni+1;
               check[k][liste[n][m]] = 1;
               m++;
            }
         }
	 
         size[n] = m;
         while ( m != 0 && k < (nk-1) )
         {
            n++;
            k++;
            liste = (int**)realloc(liste,(n+1)*sizeof(int*));
            liste[n] = (int*)(malloc)(1*sizeof(int));
            p = 0;
            
            for ( o = 0; o < m; o++)
            {
               if ( ZET[k][liste[n-1][o]] > seuil2 && check[k][liste[n-1][o]] == 0 )
               {
                  liste[n] = (int*)realloc(liste[n],(p+1)*sizeof(int));
                  liste[n][p] = liste[n-1][o];
                  check[k][liste[n][p]] = 1;
                  p++;
               }
               if ( (liste[n-1][o]+1) < (ni*nj) && ((liste[n-1][o]+1)%ni) != 0 && ZET[k][liste[n-1][o]+1] > seuil2 && check[k][liste[n-1][o]+1] == 0 )
               {
                  liste[n] = (int*)realloc(liste[n],(p+1)*sizeof(int));
                  liste[n][p] = liste[n-1][o]+1;
                  check[k][liste[n][p]] = 1;
                  p++;
               }
               if ( (liste[n-1][o]-1) > 0 && ((liste[n-1][o])%ni) != 0 && ZET[k][liste[n-1][o]-1] > seuil2 && check[k][liste[n-1][o]-1] == 0 )
               {
                  liste[n] = (int*)realloc(liste[n],(p+1)*sizeof(int));
                  liste[n][p] = liste[n-1][o]-1;
                  check[k][liste[n][p]] = 1;
                  p++;
               }
               if ( (liste[n-1][o]+ni-1) < (ni*nj) && (liste[n-1][o]%ni) != 0 && ZET[k][liste[n-1][o]+ni-1] > seuil2 && check[k][liste[n-1][o]+ni-1] == 0 )
               {
                  liste[n] = (int*)realloc(liste[n],(p+1)*sizeof(int));
                  liste[n][p] = liste[n-1][o]+ni-1;
                  check[k][liste[n][p]] = 1;
                  p++;
               }
               if ( (liste[n-1][o]+ni) < (ni*nj) && ZET[k][liste[n-1][o]+ni] > seuil2 && check[k][liste[n-1][o]+ni] == 0 )
               {
                  liste[n] = (int*)realloc(liste[n],(p+1)*sizeof(int));
                  liste[n][p] = liste[n-1][o]+ni;
                  check[k][liste[n][p]] = 1;
                  p++;
               }
               if ( (liste[n-1][o]+ni+1) < (ni*nj) && ((liste[n-1][o]+1)%ni) != 0 && ZET[k][liste[n-1][o]+ni+1] > seuil2 && check[k][liste[n-1][o]+ni+1] == 0 )
               {
                  liste[n] = (int*)realloc(liste[n],(p+1)*sizeof(int));
                  liste[n][p] = liste[n-1][o]+ni+1;
                  check[k][liste[n][p]] = 1;
                  p++;
               }
               if ( (liste[n-1][o]-ni-1) > 0 && (liste[n-1][o]%ni) != 0 && ZET[k][liste[n-1][o]-ni-1] > seuil2 && check[k][liste[n-1][o]-ni-1] == 0 )
               {
                  liste[n] = (int*)realloc(liste[n],(p+1)*sizeof(int));
                  liste[n][p] = liste[n-1][o]-ni-1;
                  check[k][liste[n][p]] = 1;
                  p++;
               }
               if ( (liste[n-1][o]-ni) > 0 && ZET[k][liste[n-1][o]-ni] > seuil2 && check[k][liste[n-1][o]-ni] == 0 )
               {
                  liste[n] = (int*)realloc(liste[n],(p+1)*sizeof(int));
                  liste[n][p] = liste[n-1][o]-ni;
                  check[k][liste[n][p]] = 1;
                  p++;
               }
               if ( (liste[n-1][o]-ni+1) > 0 && ((liste[n-1][o]+1)%ni) != 0 && ZET[k][liste[n-1][o]-ni+1] > seuil2 && check[k][liste[n-1][o]-ni+1] == 0 )
               {
                  liste[n] = (int*)realloc(liste[n],(p+1)*sizeof(int));
                  liste[n][p] = liste[n-1][o]-ni+1;
                  check[k][liste[n][p]] = 1;
                  p++;
               }
            }
            
            for ( o = 0; o < p; o++)
            {
               if ( (liste[n][o]+1) < (ni*nj) && ((liste[n][o]+1)%ni) != 0 && ZET[k][liste[n][o]+1] > seuil2 && check[k][liste[n][o]+1] == 0 )
               {
                  liste[n] = (int*)realloc(liste[n],(p+1)*sizeof(int));
                  liste[n][p] = liste[n][o]+1;
                  check[k][liste[n][p]] = 1;
                  p++;
               }
               if ( (liste[n][o]-1) > 0 && (liste[n][o]%ni) != 0 && ZET[k][liste[n][o]-1] > seuil2 && check[k][liste[n][o]-1] == 0 )
               {
                  liste[n] = (int*)realloc(liste[n],(p+1)*sizeof(int));
                  liste[n][p] = liste[n][o]-1;
                  check[k][liste[n][p]] = 1;
                  p++;
               }
               if ( (liste[n][o]+ni-1) < (ni*nj) && (liste[n][o]%ni) != 0 && ZET[k][liste[n][o]+ni-1] > seuil2 && check[k][liste[n][o]+ni-1] == 0 )
               {  
                  liste[n] = (int*)realloc(liste[n],(p+1)*sizeof(int));
                  liste[n][p] = liste[n][o]+ni-1;
                  check[k][liste[n][p]] = 1;
                  p++;
               }
               if ( (liste[n][o]+ni) < (ni*nj) && ZET[k][liste[n][o]+ni] > seuil2 && check[k][liste[n][o]+ni] == 0 )
               {
                  liste[n] = (int*)realloc(liste[n],(p+1)*sizeof(int));
                  liste[n][p] = liste[n][o]+ni;
                  check[k][liste[n][p]] = 1;
                  p++;
               }
               if ( (liste[n][o]+ni+1) < (ni*nj) && ((liste[n][o]+1)%ni) != 0 && ZET[k][liste[n][o]+ni+1] > seuil2 && check[k][liste[n][o]+ni+1] == 0 )
               {
                  liste[n] = (int*)realloc(liste[n],(p+1)*sizeof(int));
                  liste[n][p] = liste[n][o]+ni+1;
                  check[k][liste[n][p]] = 1;
                  p++;
               }
               if ( (liste[n][o]-ni-1) > 0 && (liste[n][o]%ni) != 0 && ZET[k][liste[n][o]-ni-1] > seuil2 && check[k][liste[n][o]-ni-1] == 0 )
               {         
                  liste[n] = (int*)realloc(liste[n],(p+1)*sizeof(int));
                  liste[n][p] = liste[n][o]-ni-1;
                  check[k][liste[n][p]] = 1;
                  p++;
               }
               if ( (liste[n][o]-ni) > 0 && ZET[k][liste[n][o]-ni] > seuil2 && check[k][liste[n][o]-ni] == 0 )
               {
                  liste[n] = (int*)realloc(liste[n],(p+1)*sizeof(int));
                  liste[n][p] = liste[n][o]-ni;
                  check[k][liste[n][p]] = 1;
                  p++;
               }
               if ( (liste[n][o]-ni+1) > 0 && ((liste[n][o]+1)%ni) != 0 && ZET[k][liste[n][o]-ni+1] > seuil2 && check[k][liste[n][o]-ni+1] == 0 ) 
               {
                  liste[n] = (int*)realloc(liste[n],(p+1)*sizeof(int));
                  liste[n][p] = liste[n][o]-ni+1;
                  check[k][liste[n][p]] = 1;
                  p++;
               }
            }
            
            size = (int*)realloc(size,(n+1)*sizeof(int));
            size[n] = p;
            m = p;
         }
         
         top = k;

         size1 = (int*)malloc(1*sizeof(int));
         liste1 = (int**)malloc(1*sizeof(int*));
         liste1[0] = (int*)malloc(size[0]*sizeof(int));
         
         size1[0] = size[0];
         for ( o = 0; o < size1[0]; o++)
         {
            liste1[0][o] = liste[0][o];
         }
        
         n = 0;
         k = start;
         m = size1[0];
         while ( m != 0 && k > 0 )
         {
            n++;
            k--;
            liste1 = (int**)realloc(liste1,(n+1)*sizeof(int*));
            liste1[n] = (int*)(malloc)(1*sizeof(int));

            p = 0;
            for ( o = 0; o < m; o++)
            {
               if ( ZET[k][liste1[n-1][o]] > seuil2 && check[k][liste1[n-1][o]] == 0 )
               {
                  liste1[n] = (int*)realloc(liste1[n],(p+1)*sizeof(int));
                  liste1[n][p] = liste1[n-1][o];
                  check[k][liste1[n][p]] = 1;
                  p++;
               }
               if ( (liste1[n-1][o]+1) < (ni*nj) && ((liste1[n-1][o]+1)%ni) != 0 && ZET[k][liste1[n-1][o]+1] > seuil2 && check[k][liste1[n-1][o]+1] == 0 )
               {
                  liste1[n] = (int*)realloc(liste1[n],(p+1)*sizeof(int));
                  liste1[n][p] = liste1[n-1][o]+1;
                  check[k][liste1[n][p]] = 1;
                  p++;
               }
               if ( (liste1[n-1][o]-1) > 0 && (liste1[n-1][o]%ni) != 0 && ZET[k][liste1[n-1][o]-1] > seuil2 && check[k][liste1[n-1][o]-1] == 0 )
               {
                  liste1[n] = (int*)realloc(liste1[n],(p+1)*sizeof(int));
                  liste1[n][p] = liste1[n-1][o]-1;
                  check[k][liste1[n][p]] = 1;
                  p++;
               }
               if ( (liste1[n-1][o]+ni-1) < (ni*nj) && (liste1[n-1][o]%ni) != 0 && ZET[k][liste1[n-1][o]+ni-1] > seuil2 && check[k][liste1[n-1][o]+ni-1] == 0 )
               {
                  liste1[n] = (int*)realloc(liste1[n],(p+1)*sizeof(int));
                  liste1[n][p] = liste1[n-1][o]+ni-1;
                  check[k][liste1[n][p]] = 1;
                  p++;
               }
               if ( (liste1[n-1][o]+ni) < (ni*nj) && ZET[k][liste1[n-1][o]+ni] > seuil2 && check[k][liste1[n-1][o]+ni] == 0 )
               {
                  liste1[n] = (int*)realloc(liste1[n],(p+1)*sizeof(int));
                  liste1[n][p] = liste1[n-1][o]+ni;
                  check[k][liste1[n][p]] = 1;
                  p++;
               }
               if ( (liste1[n-1][o]+ni+1) < (ni*nj) && ((liste1[n-1][o]+1)%ni) != 0 && ZET[k][liste1[n-1][o]+ni+1] > seuil2 && check[k][liste1[n-1][o]+ni+1] == 0 )
               {
                  liste1[n] = (int*)realloc(liste1[n],(p+1)*sizeof(int));
                  liste1[n][p] = liste1[n-1][o]+ni+1;
                  check[k][liste1[n][p]] = 1;
                  p++;
               }
               if ( (liste1[n-1][o]-ni-1) > 0 && (liste1[n-1][o]%ni) != 0 && ZET[k][liste1[n-1][o]-ni-1] > seuil2 && check[k][liste1[n-1][o]-ni-1] == 0 )
               {
                  liste1[n] = (int*)realloc(liste1[n],(p+1)*sizeof(int));
                  liste1[n][p] = liste1[n-1][o]-ni-1;
                  check[k][liste1[n][p]] = 1;
                  p++;
               }
               if ( (liste1[n-1][o]-ni) > 0 && ZET[k][liste1[n-1][o]-ni] > seuil2 && check[k][liste1[n-1][o]-ni] == 0 )
               {
                  liste1[n] = (int*)realloc(liste1[n],(p+1)*sizeof(int));
                  liste1[n][p] = liste1[n-1][o]-ni;
                  check[k][liste1[n][p]] = 1;
                  p++;
               }
               if ( (liste1[n-1][o]-ni+1) > 0 && ((liste1[n-1][o]+1)%ni) != 0 && ZET[k][liste1[n-1][o]-ni+1] > seuil2 && check[k][liste1[n-1][o]-ni+1] == 0 )
               {
                  liste1[n] = (int*)realloc(liste1[n],(p+1)*sizeof(int));
                  liste1[n][p] = liste1[n-1][o]-ni+1;
                  check[k][liste1[n][p]] = 1;
                  p++;
               }
            }
            for ( o = 0; o < p; o++)
            {
               if ( (liste1[n][o]+1) < (ni*nj) && ((liste1[n][o]+1)%ni) != 0 && ZET[k][liste1[n][o]+1] > seuil2 && check[k][liste1[n][o]+1] == 0 )
               {
                  liste1[n] = (int*)realloc(liste1[n],(p+1)*sizeof(int));
                  liste1[n][p] = liste1[n][o]+1;
                  check[k][liste1[n][p]] = 1;
                  p++;
               }
               if ( (liste1[n][o]-1) > 0 && (liste1[n][o]%ni) != 0 && ZET[k][liste1[n][o]-1] > seuil2 && check[k][liste1[n][o]-1] == 0 )
               {
                  liste1[n] = (int*)realloc(liste1[n],(p+1)*sizeof(int));
                  liste1[n][p] = liste1[n][o]-1;
                  check[k][liste1[n][p]] = 1;
                  p++;
               }
               if ( (liste1[n][o]+ni-1) < (ni*nj) && ((liste1[n][o]-1)%ni) != 0 && ZET[k][liste1[n][o]+ni-1] > seuil2 && check[k][liste1[n][o]+ni-1] == 0 )
               {  
                  liste1[n] = (int*)realloc(liste1[n],(p+1)*sizeof(int));
                  liste1[n][p] = liste1[n][o]+ni-1;
                  check[k][liste1[n][p]] = 1;
                  p++;
               }
               if ( (liste1[n][o]+ni) < (ni*nj) && ZET[k][liste1[n][o]+ni] > seuil2 && check[k][liste1[n][o]+ni] == 0 )
               {
                  liste1[n] = (int*)realloc(liste1[n],(p+1)*sizeof(int));
                  liste1[n][p] = liste1[n][o]+ni;
                  check[k][liste1[n][p]] = 1;
                  p++;
               }
               if ( (liste1[n][o]+ni+1) < (ni*nj) && ((liste1[n][o]+1)%ni) != 0 && ZET[k][liste1[n][o]+ni+1] > seuil2 && check[k][liste1[n][o]+ni+1] == 0 )
               {
                  liste1[n] = (int*)realloc(liste1[n],(p+1)*sizeof(int));
                  liste1[n][p] = liste1[n][o]+ni+1;
                  check[k][liste1[n][p]] = 1;
                  p++;
               }
               if ( (liste1[n][o]-ni-1) > 0 && (liste1[n][o]%ni) != 0 && ZET[k][liste1[n][o]-ni-1] > seuil2 && check[k][liste1[n][o]-ni-1] == 0 )
               {         
                  liste1[n] = (int*)realloc(liste1[n],(p+1)*sizeof(int));
                  liste1[n][p] = liste1[n][o]-ni-1;
                  check[k][liste1[n][p]] = 1;
                  p++;
               }
               if ( (liste1[n][o]-ni) > 0 && ZET[k][liste1[n][o]-ni] > seuil2 && check[k][liste1[n][o]-ni] == 0 )
               {
                  liste1[n] = (int*)realloc(liste1[n],(p+1)*sizeof(int));
                  liste1[n][p] = liste1[n][o]-ni;
                  check[k][liste1[n][p]] = 1;
                  p++;
               }
               if ( (liste1[n][o]-ni+1) > 0 && ((liste1[n][o]+1)%ni) != 0 && ZET[k][liste1[n][o]-ni+1] > seuil2 && check[k][liste1[n][o]-ni+1] == 0 ) 
               {
                  liste1[n] = (int*)realloc(liste1[n],(p+1)*sizeof(int));
                  liste1[n][p] = liste1[n][o]-ni+1;
                  check[k][liste1[n][p]] = 1;
                  p++;
               }
            }
            size1 = (int*)realloc(size1,(n+1)*sizeof(int));
            size1[n] = p;
            m = p;
         }
         base = k;
         
         size2 = (int*)malloc((top-base-1)*sizeof(int));
         liste2 = (int**)malloc((top-base-1)*sizeof(int*));
         for ( o = (start-base-1); o > 0; o--)
         {
            liste2[start-base-1-o] = (int*)malloc(size1[o]*sizeof(int));
            size2[start-base-1-o] = size1[o];
            for ( p = 0; p < size1[o]; p++ )
            {
               liste2[start-base-1-o][p] = liste1[o][p];
            }
         }
         for ( o = 0; o < (top-start); o++)
         {
            liste2[start-base-1+o] = (int*)malloc(size[o]*sizeof(int));
            size2[start-base-1+o] = size[o];

            for ( p = 0; p < size[o]; p++ )
            {
               liste2[start-base-1+o][p] = liste[o][p];
            }
         }
         
         for ( o = 0; o < (top-start+1); o++)
         {
            free(liste[o]);
         }
         for ( o = 0; o < (start-base+1); o++)
         {
            free(liste1[o]);
         }
        
         free(size);
         free(size1);
         free(liste);
         free(liste1);
         size = NULL;
         size1 = NULL;
         liste = NULL;
         liste1 = NULL;

         Surplomb = 0;
         if ( (GZ[top-1][i] - GZ[base+1][i]) > 500. )
         {
            for ( o = 0; o < (top-base-1); o++ )
            {
               for ( p = 0; p < size2[o]; p++ )
               {
                  k = base;
		  Surplomb = 1;
                  while ( check[k][liste2[o][p]] != 1 && (GZ[k][liste2[o][p]]- GZ[base][liste2[o][p]]) < 300. )
                  {
                     if ( ZET[k][liste2[o][p]] > 10. )
                        Surplomb = 0;
                     k++;
                  }
                  
                  if ( Surplomb == 1 )
		    goto end;
               }
            }
	 }
         end:
         if (  Surplomb == 1 )
         {
            k_max = 0;
            i_max = 0;

            for (o = 0; o < n; o++ )
            {
               for ( p = 0; p < size2[o]; p++ )
               {
                  zoneSurplomb[liste2[o][p]] = 1;
                  if (  MaxR[liste2[o][p]] >= 40. && res <= 1.0 )
                  {
                     k = 0;
                     while ( ZET[k][liste2[o][p]] <= 0. )
                     {
                        if ( k > k_max )
                        {
                           k_max = k;
                           i_max = liste2[o][p];
                        }
                        k++;
                     }
                  }
	       }
            }
            z++;
             
	    if ( k_max > 3 && res <= 1.0 )
            {
               k = k_max;
               while ( ZET[k][i_max] <= 10. )
               {
                  if ( k > k_max )
                     k_max = k;
                  k++;
               }
               
               for ( o = 1; o < 5; o++ )
               {     
                  if ( (i_max - o) >= (floor(i_max/ni)*ni) && ZET[k_max][i_max - o] >= 30. && check[k_max][i_max - o] == 0 )
                  {
                     BWER[i_max] = 1; 
                     break;
                  }
                  if ( (i_max + o) < (ceil(i_max/ni)*ni) && ZET[k_max][i_max + o] >= 30. && check[k_max][i_max + o] == 0 )
                  {
                     BWER[i_max] = 1; 
                     break;
                  }
                  if ( (i_max + (o*ni)) < (ni*nj) && ZET[k_max][i_max + (o*ni)] >= 30. && check[k_max][i_max + (o*ni)] == 0 )
                  {
                     BWER[i_max] = 1; 
                     break;
                  }
                  if ( (i_max - (o*ni)) >= 0 && ZET[k_max][i_max - (o*ni)] >= 30. && check[k_max][i_max - (o*ni)] == 0 )
                  {
                     BWER[i_max] = 1; 
                     break;
                  }
                  if ( (i_max + (o*ni) + o) < (ni*nj) && (i_max + o) < (ceil(i_max/ni)*ni) && ZET[k_max][i_max + (o*ni) + o] >= 30. && check[k_max][i_max + (o*ni) + o] == 0 )
                  {
                     BWER[i_max] = 1; 
                     break;
                  }
                  if ( (i_max + (o*ni) - o) < (ni*nj) && (i_max - o) >= (floor(i_max/ni)*ni) && ZET[k_max][i_max + (o*ni) - o] >= 30. && check[k_max][i_max + (o*ni) - o] == 0 )
                  {
                     BWER[i_max] = 1; 
                     break;
                  }
                  if ( (i_max - (o*ni) + o) >= 0 && (i_max + o) < (ceil(i_max/ni)*ni) && ZET[k_max][i_max - (o*ni) + o] >= 30. && check[k_max][i_max - (o*ni) + o] == 0 )
                  {
                     BWER[i_max] = 1; 
                     break;
                  }
                  if ( (i_max - (o*ni) - o) >= 0 && (i_max - o) >= (floor(i_max/ni)*ni) && ZET[k_max][i_max - (o*ni) - o] >= 30. && check[k_max][i_max - (o*ni) - o] == 0 )
                  {
                     BWER[i_max] = 1; 
                     break;
                  }
                  if ( o != 1 )
                  {
                     p = 1;
                     while ( p != o )
                     {
                        if ( (i_max + (p*ni) + o) < (ni*nj) && (i_max +(p*ni) +o) < (ceil((i_max+(p*ni))/ni)*ni) && ZET[k_max][i_max + (p*ni) + o] >= 30. && check[k_max][i_max + (p*ni) + o] == 0 )
                        {
                           BWER[i_max] = 1; 
                           break;
                        }
                        if ( (i_max + (p*ni) - o) < (ni*nj) && (i_max +(p*ni) -o) >= (floor((i_max+(p*ni))/ni)*ni) && ZET[k_max][i_max + (p*ni) - o] >= 30. && check[k_max][i_max + (p*ni) - o] == 0 )
                        {
                           BWER[i_max] = 1; 
                           break;
                        }
                        if ( (i_max - (p*ni) + o) >= 0 && (i_max - (p*ni) + o) < (ceil((i_max-(p*ni))/ni)*ni) && ZET[k_max][i_max - (p*ni) + o] >= 30. && check[k_max][i_max - (p*ni) + o] == 0 )
                        {
                           BWER[i_max] = 1; 
                           break;
                        }
                        if ( (i_max - (p*ni) - o) >= 0 && (i_max -(p*ni) -o) >= (floor((i_max-(p*ni))/ni)*ni) && ZET[k_max][i_max - (p*ni) - o] >= 30. && check[k_max][i_max - (p*ni) - o] == 0 )
                        {
                           BWER[i_max] = 1; 
                           break;
                        }
                        if ( (i_max + (o*ni) - p) < (ni*nj) && (i_max + (o*ni) - p) >= (floor((i_max + (o*ni))/ni)*ni) && ZET[k_max][i_max + (o*ni) - p] >= 30. && check[k_max][i_max + (o*ni) - p] == 0 )
                        {
                           BWER[i_max] = 1; 
                           break;
                        }
                        if ( (i_max - (o*ni) - p) >= 0 && (i_max - (o*ni) - p) >= (floor((i_max - (o*ni))/ni)*ni) && ZET[k_max][i_max - (o*ni) - p] >= 30. && check[k_max][i_max - (o*ni) - p] == 0 )
                        {
                           BWER[i_max] = 1; 
                           break;
                        }
                        if ( (i_max + (o*ni) + p) < (ni*nj) && (i_max + (o*ni) + p) < (ceil((i_max+(o*ni))/ni)*ni) && ZET[k_max][i_max + (o*ni) + p] >= 30. && check[k_max][i_max + (o*ni) + p] == 0 )
                        {
                           BWER[i_max] = 1; 
                           break;
                        }
                        if ( (i_max - (o*ni) + p) >= 0 && (i_max - (o*ni) + p) < (ceil((i_max-(o*ni))/ni)*ni) && ZET[k_max][i_max - (o*ni) + p] >= 30. && check[k_max][i_max - (o*ni) + p] == 0 )
                        {
                           BWER[i_max] = 1; 
                           break;
                        }
                        p++;
                     }
                  }
               }
            }
         }
         
         for ( o = 0; o < (top-base-1); o++)
         {
            free(liste2[o]);
         }
         free(size2);
         free(liste2);
         size2 = NULL;
         liste2 = NULL;
      }
   }

   //-----------------------
   // Liberer la memoire
   //-----------------------
   for ( k = 0; k < nk; k++ )
   {
      free(check[k]);
   }
   free(check);
   check = NULL;
   
   //------------------
   // Envoyer resultats
   //------------------
   *oBWER = BWER;
   *ozoneSurplomb = zoneSurplomb;

   return SUCCESS;
}


/*!
==============================================================================
\brief Calculer le transport vertical d'humidite
\date 08 mars 2010
\version 0.1
\author Anna-Belle Filion

param[in]  WW             : Tableau de mouvement verticaux (Pa/s)
param[in]  HU             : Tableu d'humidite specifique (kg/Kg)
param[in]  ni             : Dimension en x
param[in]  nj             : Dimension en y
param[in]  nbNiveau       : Dimension en z
param[out] VMF            : Tableau de transport vertical d'humidite 

------------------------------------------------------------------------------
*/
int CalculerVerticalMoistureFlux(float **WW, float **HU, int ni, int nj, int nbNiveau,
                                 float **VMF, float **VMF_NEW)
{
   int i, k;
   
   float *VerticalMoistureFlux = (float*)NULL;
   float *VerticalMoistureFlux_NEW = (float*)NULL;
   
   //----------------------
   // Allouee la memoire
   //----------------------
   VerticalMoistureFlux = (float*)malloc(ni*nj*sizeof(float));
   VerticalMoistureFlux_NEW = (float*)malloc(ni*nj*sizeof(float));
   
   //--------------------------------------------------
   // Trouver le transport vertical total d'humidite
   //--------------------------------------------------
   for( i = 0; i < (ni*nj); i++ )
   {
      VerticalMoistureFlux[i] = 0;
      VerticalMoistureFlux_NEW[i] = 0;
      for ( k = 0; k < nbNiveau; k++)
      {
	 if ( WW[k][i] < 0. && HU[k][i] > 0 )  
               VerticalMoistureFlux[i] += fabs((double)(WW[k][i] * HU[k][i] *1000.));
      }
      for ( k = 0; k < nbNiveau; k++)
      {
	 if ( HU[k][i] < 0.005 )
	    break;
	 else if ( WW[k][i] < 0. && HU[k][i] > 0 )
               VerticalMoistureFlux_NEW[i] += fabs((double)(WW[k][i] * HU[k][i] *1000.));
      }
      if ( k > 1 )
      VerticalMoistureFlux_NEW[i] = VerticalMoistureFlux_NEW[i]/(k-1); 
   }

   //--------------------------
   // Retourner les resultats
   //--------------------------
   *VMF = VerticalMoistureFlux;
   *VMF_NEW = VerticalMoistureFlux_NEW;
   
   return SUCCESS;
}



/*!
=============================================================================
\brief Calculer le flux vertical de glace ( graupel + cristaux + hail)
\date 07 mars 2010
\version 0.1
\author Anna-Belle Filion

param[in]  TT             : Tableau de temperature (C)
param[in]  WW             : Tableau de mouvement verticaux (Pa/s)
param[in]  QJT1           : Tableau de rapport de melange de graupel (kg/Kg)
param[in]  QIT1           : Tableau de rapport de melange de cristaux de glace (kg/Kg)
param[in]  QHT1           : Tableau de rapport de melange de grele (kg/Kg)
param[in]  ni             : Dimension en x
param[in]  nj             : Dimension en y
param[in]  nbNiveau       : Dimension en z
param[out] VIF            : Vecteur qui contient le flu vertical de glace

-----------------------------------------------------------------------------
*/
int CalculerVerticalIceFlux(float **TT,float **WW, float **QJT1, float **QIT1, 
                    float **QHT1, int ni, int nj, int nbNiveau, float **VIF)
{
   int i, k;
   
   float *VerticalIceFlux;
   
   //----------------------
   // Allouee la memoire
   //----------------------
   VerticalIceFlux = (float*)malloc(ni*nj*sizeof(float));
   
   //------------------------------------------------
   // Calculer le transport vertical total de glace
   //------------------------------------------------
   for ( i = 0; i < (ni*nj); i++)
   {
      VerticalIceFlux[i] = 0;
      for ( k = 1; k < nbNiveau; k++)
      {
         if ( TT[k-1][i] > -15. && TT[k][i] < -15.)
         {
            if ( WW[k][i] < 0 && ( QIT1[k][i] != 0 || QJT1[k][i] != 0 || QHT1[k][i] != 0) )
               VerticalIceFlux[i] += fabs((double)((QJT1[k][i] + QIT1[k][i] + QHT1[k][i]) * fabs((double)(WW[k][i])) * 1000.));
         }
      }
   }
   
   //--------------------------
   // Retourner les resultats
   //--------------------------
   *VIF = VerticalIceFlux;
   
   return SUCCESS;
}


/*!
===============================================================================
\brief Calculer le cisaillement entre la surface et 3km
       Calculer le cisaillement entre la surface et 6km
\date 06 mars 2010
\version 0.1
\author Anna-Belle Filion

param[in]  GZ               : Tableau de hauteur geopotentielle (DAM)
param[in]  WD               : Tableau de la direction du vent (degree meteo)
param[in]  UV               : Tableau du module du vent (noeuds)
param[in]  ni               : Dimension en x
param[in]  nj               : Dimension en y
param[in]  nbNiveau         : Dimension en z
param[out] cisaillement_3km : Vecteur de cisaillement entre la surface et 3km
param[out] cisaillement_6km : Vecteur de cisaillement entre la surface et 6km

-------------------------------------------------------------------------------
*/
int CalculerCisaillement_3kmET6km( float **GZ, float **WD, float **UV, int ni, int nj,
                                   int nbNiveau, float **cisaillement_3km, 
                                   float **cisaillement_6km)
{
   int i, k;
   
   float A, Carre_Norme, Ci;
   
   float *Cisaillement3km = (float*)NULL;
   float *Cisaillement6km = (float*)NULL;
   
   //----------------------
   // Allouee la memoire
   //----------------------
   Cisaillement3km = (float*)malloc(ni*nj*sizeof(float));
   Cisaillement6km = (float*)malloc(ni*nj*sizeof(float));
   
   //------------------------------------------------------
   // Calculer le cisaillement total (direction+vitesse)
   // entre la surface et 3km AGL et la surface et 6km AGL
   //------------------------------------------------------
   for ( i = 0; i < (ni*nj); i++)
   {
      Cisaillement3km[i] = 0.;
      Cisaillement6km[i] = 0.;
      
      for( k = 1; k < nbNiveau; k++)
      {
         if ( (GZ[k][i]-GZ[0][i]) <= 600. )
         {
            A = WD[k][i] - WD[k-1][i];
            
            Carre_Norme = pow(UV[k-1][i], 2.) + pow(UV[k][i], 2.) - 2.* UV[k-1][i] * UV[k][i] * cos(A*PI/180.);
            
            if ( Carre_Norme > 0.0 )
               Ci = sqrt( Carre_Norme );
            else 
               Ci = 0.;
            
            Cisaillement6km[i] = Cisaillement6km[i] + Ci;
            if ( (GZ[k][i]-GZ[0][i]) <= 300. )
               Cisaillement3km[i] = Cisaillement3km[i] + Ci;
         }
      }
      //Ne pas diviser par la hauteur car on veut le cisaillement totaux dans la couche si > 35 pour 6 km ==> cisaillement elever et cisaillemt doit etre en kts
      Cisaillement6km[i] = Cisaillement6km[i];
      Cisaillement3km[i] = Cisaillement3km[i];
   }
   
   //--------------------------
   // Retourner les resultats
   //--------------------------
   *cisaillement_3km = Cisaillement3km;
   *cisaillement_6km = Cisaillement6km;
   
   return SUCCESS;
}



/*!
==========================================================================
\brief Rechercher la reflectivite maximum dans une colonne qui depasse 
       un certain seuil et se trouve au dessus du point de congelation. 

\date 23 fevrier 2010
\version 0.2
\author  Anna-Belle Filion

param[in] ZET          : Tableau de reflectivite           (dBZ)
param[in] GZ           : Tableau de hauteur geopotentielle (DAM)
param[in] indiceT0     : Vecteur d'indices du point de congelation
                         dans la colonne
param[in] ni           : Dimension en x
param[in] nj           : Dimension en y
param[in] nbNiveau     : Dimension en z
param[in] seuil        : Seuil de reflectivite minimale    (dBZ)
param[out] maxR        : Vecteur de reflectivite maximum qui depasse le seuil (dBZ)
param[out] indiceMaxR  : Vecteur d'indice de la hauteur de reflectivite 
                         maximum dans la colonne
param[out] hauteurMaxR : Vecteur de hauteur maximum ou la reflectivite 
                         depasse le seuil (m)

--------------------------------------------------------------------------
*/
int TrouverMaxRAboveT0( float **ZET, float **GZ, int *indiceT0, int ni, int nj,
                        int nbNiveau, float seuil, float **maxR, int **indiceMaxR,
                        float **hauteurMaxR)
{
   int i, k;
   
   int *IndiceMaxR = (int*)NULL;
   
   float *MaxR        = (float*)NULL;
   float *HauteurMaxR = (float*)NULL;
   
   //--------------------------
   // Allocation de memoire
   //--------------------------
   MaxR        = (float*)malloc(ni*nj*sizeof(float));
   HauteurMaxR = (float*)malloc(ni*nj*sizeof(float));
   IndiceMaxR  = (int*)malloc(ni*nj*sizeof(int));
   
   //-----------------------------------------------
   // Trouver la reflectivite maximale au dessus 
   // du niveau de congelation
   //-----------------------------------------------
   for ( i = 0; i < (ni*nj); i++)
   {
      MaxR[i] = 0.;
      HauteurMaxR[i] = 0.;
      IndiceMaxR[i] = 0;
      
      for( k = (indiceT0[i]+1); k < nbNiveau; k++)
      {
         if ( ZET[k][i] > MaxR[i] && ZET[k][i] >= seuil)
         {
            MaxR[i] = ZET[k][i];
            IndiceMaxR[i] = k;
            HauteurMaxR[i] = (GZ[k][i] - GZ[0][i]) * 10.;
         }
      }
   }
   
   //--------------------------
   // Retourner les resultats
   //--------------------------
   *maxR = MaxR;
   *hauteurMaxR = HauteurMaxR;
   *indiceMaxR  = IndiceMaxR;

   return SUCCESS;
}

/*!
=========================================================================
\brief Trouver le niveau de congelation pour chaque point
\date 08 mars 2010
\version 0.1
\author Anna-Belle Filion

param[in]  TT         : Tableau de temperatures (C)
param[in]  ni         : Dimension en x
param[in]  nj         : Dimension en y
param[in]  nbNiveau   : Dimension en z
param[out] indiceT0   : Vecteur d'indices qui correspond a l'indice
                        du niveau de congelation

-------------------------------------------------------------------------
*/
int TrouverNiveauCongelation(float **TT, int ni, int nj, int nbNiveau,
                             int **indiceT0)
{
   int i, k;
   
   int *IndiceT0 = (int*)NULL;
   
   //--------------------------
   // Allocation de memoire
   //--------------------------
   IndiceT0 = (int*)malloc(ni*nj*sizeof(int));
   
   //-------------------------------------
   // Trouver le niveau de congelation
   //-------------------------------------
   for ( i = 0; i < (ni*nj); i++)
   {
      IndiceT0[i] = 0;
      
      for( k = 1; k < nbNiveau; k++)
      {
         if ( TT[k][i] == 0. )
         {
            IndiceT0[i] = k;
            break;
         }
         else if ( TT[k-1][i] > 0. && TT[k][i] < 0.0 )
         {
            IndiceT0[i] = k;
            break;
         }
      }
   }
   
   //--------------------------
   // Retourner les resultats
   //--------------------------
   *indiceT0 = IndiceT0;
   
   return SUCCESS;
}


