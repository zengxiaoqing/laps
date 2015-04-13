#ifndef CONSTANTE_MATH_H
#define CONSTANTE_MATH_H


#define		T_ZERO_C	               	273.15
#define         T_TRIPLE_POINT                  273.16
#define	        PI			        3.141592653589793238462643
#define         PI_DIV_2                        1.570796326794896619231322
#define         PI_DIV_4                        0.7853981633974483096156608
#define         DEGRE_A_RADIAN                  0.017453292519943295769237
#define         RADIAN_A_DEGRE                  57.295779513082320876798155
#define	        METRES_PAR_SEC_EN_NOEUDS        1.94254
#define         NOEUDS_A_METRES_PAR_SEC         1./METRES_PAR_SEC_EN_NOEUDS
#define         PIEDS_A_METRES                  .3048
#define         METRES_A_PIEDS                  1./PIEDS_A_METRES
#define         PRES_POUCE_MB                   33.86395
#define         NOEUDS_A_KMH                    1.853244
#define         KMH_A_NOEUDS                    1./NOEUDS_A_KMH
#define         KMH_A_METRES_PAR_SEC            0.2777777778
#define         METRES_PAR_SEC_A_KMH            1./KMH_A_METRES_PAR_SEC
#define         PRES_MB_POUCE                   .02952993
#define		CONSTANTE_P0	                1013.246 	/* mb */
#define	        CONSTANTE_RD	 	        2.8705e02      /*J Kg^-1 K^-1*/
#define	        CONSTANTE_CPD                   1005.   /*J Kg^-1 K^-1*/
#define         CONSTANTE_RV                    461.51          /* J/kg K */
#define		CONSTANTE_R_SUR_CP	        0.2859063745	/* sans unite */
#define		CONSTANTE_G		        9.80665		/* m s^-2  */
#define         CONSTANTE_LS                    2.8345e+6        /* J kg^-1 */
#define         CONSTANTE_LV                    2.5008e+6        /* J kg^-1 */
#define         CONSTANTE_LF                    0.3337e+6        /* J kg^-1 */
#define		EPSILON		                0.62197     /* Rd/Rv */
#define		GAMMA_SEC		        -CONSTANTE_G/CONSTANTE_CPD
#define		UNITE_LV	                1.0e+06		/* en J */
#define		UNITE_CALORIE_A_JOULE	        4.1868 		/* cal -> J */
#define		METRE_A_METRE_GEOPOTENTIEL      1./0.980665 	/* m -> gpm */
#define         VITESSE_LUMIERE                 2.997925e+8     /* m/s/s  */
#define         CONSTANTE_GRAVITE_G             6.670e-11       /* N*m*m/kg/kg */
#define         PASCAL_A_MBAR                   .01
#define         MBAR_A_PASCAL                   100.
#define         NM_A_METRE                      1852.
#define         METRE_A_NM                      1./NM_A_METRE
#define         CTE_AVOGADRO                    6.02252e+23     /* mole */
#define         CTE_LOG_NEPER                   2.718281828459045235360287
#define         LAPSE_RATE_DIV_TEMP_SEA_LEVEL   2.255e-5
#define         G0_DIV_RD_B0                    5.256
#define         INV_G0_DIV_RD_B0                0.19025875  /*  1/ 5.256 */
#define	        DEFAUT_VALEUR_NUMERIQUE	   	-99999.0

/*
--- Definitions des types de parametres attendus et des references pour les
    interpolations
*/
#define		PARAMETRE_TT			0
#define		PARAMETRE_TD			1
#define		PARAMETRE_VM			2
#define		PARAMETRE_UU			3
#define		PARAMETRE_VV			4

#define		REFERENCE_PP			0
#define		REFERENCE_PP_WW			1
#define		REFERENCE_HH			2

/*
--- Definition des nombres de niveaux dans les tephis ou pour le calcul
*/
#define		MAX_NIVEAUX_CALCUL		 7000
#define		LARGEUR_NIVEAU_CALCUL_EN_PIEDS	   10
#define		LARGEUR_MOYENNE_CALCUL_EN_MANAIR   10

#define		MAX_NIVEAUX_TEPHI		  200

#define		MAX_PRECISION_CALCUL		 0.01
#define		MAX_PRECISION_ITERATION		0.001



/*
--- Definitions des limites entre lesquelles on fait les calculs
*/
#define		LIMITE_PRESSION_MIN		 100.0	/* en mb */
#define		LIMITE_PRESSION_MAX		1100.0  /* en mb */
#define		LIMITE_ANGLE_MIN		   0.0	/* degres */
#define		LIMITE_ANGLE_MAX		 360.0  /* degres */
#define		THETA_REF_SCALE			210.0   /* degres */
#define		T_REF_SCALE			-60.0   /* degres celsius */
#define		T_LIMITE_ADIABATIQUE_HUMIDE	-50.0   /* degres celsius */

/*
--- Definitions de valeurs generales pour les erreurs ou de decision
*/
#define 	SUCCESS 			0
#define         NO_ERROR                        0
#define 	ERROR				1
#define         FALSE                           0
#define         TRUE                            1



#endif
