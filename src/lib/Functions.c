#include <stdio.h>
#include <string.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

#include "Structures.h"
#include "Define.h"

#define     OUI                             1

 /*!
=================================================================================
\brief Lire le fichier de configuration

\date 11 juin 2013

\version 0.1
\author Anna-Belle Filion

\param[in]  config_file  : nom et chemin du fichier de configuration
\param[out] parametre    : Parametre a configurer

---------------------------------------------------------------------------------
 */
int read_config_file_( int DEBUG, char *fichier_cfg, CFG *cfg)
{
	if ( DEBUG == OUI )
		fprintf(stderr,"\nFCN ==> Read_Config_File\n\n");
	
   FILE *Config = (FILE*)NULL;
   
   char line[CARAC_MAX]    = "";
   char element[CARAC_MAX/5] = "";
 
   if ( (Config = fopen(fichier_cfg, "r")) != NULL ) 
   {
      while( ! feof(Config) )
      {
         fgets(line, CARAC_MAX, Config);
         if ( line == "" || line[0] == '#' )
				continue;
			else
			{
            sscanf(line, "%s %*s", element);
            if ( strcmp("IP2",element) == 0 )
              sscanf(line, "%*s %d", &cfg->ip2);
            else if ( strcmp("NPAS",element) == 0 )
              sscanf(line, "%*s %d", &cfg->npas);
            else if ( strcmp("HH_PREVU",element) == 0 )
              sscanf(line, "%*s %d", &cfg->hh_fcst);
            else if ( strcmp("NBR_NIVEAU_DISPO",element) == 0 )
              sscanf(line, "%*s %d", &cfg->nbr_level_avail);
            else if ( strcmp("NBR_NIVEAU_VOULU",element) == 0 )
              sscanf(line, "%*s %d", &cfg->nbr_level_2keep);
            else if ( strcmp("RESOLUTION",element) == 0 )
              sscanf(line, "%*s %f", &cfg->resolution);
            else if ( strcmp("SEUIL_MESO",element) == 0 )
              sscanf(line, "%*s %f", &cfg->thres_meso);
            else if ( strcmp("SEUIL_TORNADE",element) == 0 )
              sscanf(line, "%*s %f", &cfg->thres_tornade);
            else if ( strcmp("SEUIL_SEVERE_UPDRAFT_MAX",element) == 0 )
              sscanf(line, "%*s %f", &cfg->thres_severe_updraft[0]);
            else if ( strcmp("SEUIL_SEVERE_UPDRAFT_ENVELOPPE",element) == 0 )
              sscanf(line, "%*s %f", &cfg->thres_severe_updraft[1]);
            else if ( strcmp("SEUIL_VORTEX_NEGATIVE_MAX",element) == 0 )
              sscanf(line, "%*s %f", &cfg->thres_couplet_tourbillon[0][0]);
            else if ( strcmp("SEUIL_VORTEX_NEGATIVE_ENVELOPPE",element) == 0 )
              sscanf(line, "%*s %f", &cfg->thres_couplet_tourbillon[0][1]);
            else if ( strcmp("SEUIL_VORTEX_POSITIF_MAX",element) == 0 )
              sscanf(line, "%*s %f", &cfg->thres_couplet_tourbillon[1][0]);
            else if ( strcmp("SEUIL_VORTEX_POSITIF_ENVELOPPE",element) == 0 )
              sscanf(line, "%*s %f", &cfg->thres_couplet_tourbillon[1][1]);
				else if ( strcmp("YYYYMMDDHH",element) == 0 )
					sscanf(line, "%*s %d", &cfg->yyyymmddhh);
            else if ( strcmp("HEURE_INITIALISTION",element) == 0 )
              sscanf(line, "%*s %d", &cfg->hh_init);
				else if ( strcmp("REGION",element) == 0 )
					sscanf(line, "%*s %s", cfg->region);
				else if ( strcmp("EXTENSION_FICHIER_SORTIE",element) == 0 )
					sscanf(line, "%*s %s", cfg->ext_file_out);
            else if ( strcmp("EXTENSION_FICHIER_ENTREE",element) == 0 )
              sscanf(line, "%*s %s", cfg->ext_file_in);
				else if ( strcmp("PATH_FICHIER_MODELE",element) == 0 )
					sscanf(line, "%*s %s", cfg->path_file_in);
				else if ( strcmp("PATH_FICHIER_SORTIE",element) == 0 )
					sscanf(line, "%*s %s", cfg->path_file_out);
				else if ( strcmp("LISTENOMVAR3D",element) == 0 )
					sscanf(line, "%*s %s", cfg->varname3D);
				else if ( strcmp("LISTENOMVAR2D",element) == 0 )
					sscanf(line, "%*s %s", cfg->varname2D);
         }
      }
   } 
   else 
   {
      fprintf(stderr,"ERREUR lors de la lecture du fichier de configuration <%s>\n",fichier_cfg);
      exit(1);
   }
  
  if ( DEBUG == OUI )
		fprintf(stderr,"npas = %d\nip2 = %d\nhh_fcst = %d\nnbr_level_2keep = %d\nnbr_level_avail = %d\nnbr_level_2keep = %d\ncfg->resolution = %f\n\
thres_meso = %f\nthres_tornade = %f\nyyyymmddhh = %d\nhh_init = %d\nregion = %s\next_file_out = %s\next_file_in = %s\npath_file_in = %s\n\
path_file_out = %s\nvarname3D = %s\nvarname2D = %s\nthres_severe_updraft[0] = %f\nthres_severe_updraft[1] = %f\nthres_couplet_tourbillon[0][0] = %f\n\
thres_couplet_tourbillon[0][1] = %f\nthres_couplet_tourbillon[1][0] = %f\nthres_couplet_tourbillon[1][1] = %f\n\n", cfg->npas, \
		cfg->ip2, cfg->hh_fcst, cfg->nbr_level_2keep, cfg->nbr_level_avail, cfg->nbr_level_2keep, cfg->resolution, cfg->thres_meso, cfg->thres_tornade, \
		cfg->yyyymmddhh, cfg->hh_init, cfg->region, cfg->ext_file_out, cfg->ext_file_in, cfg->path_file_in, cfg->path_file_out, cfg->varname3D, cfg->varname2D, \
		cfg->thres_severe_updraft[0], cfg->thres_severe_updraft[1], cfg->thres_couplet_tourbillon[0][0], cfg->thres_couplet_tourbillon[0][1], \
		cfg->thres_couplet_tourbillon[1][0], cfg->thres_couplet_tourbillon[1][1]);

  if ( DEBUG == OUI )
		fprintf(stderr,"FCN ==> Read_Config_File completed iwth success\n\n");
  
  return SUCCESS;
}
