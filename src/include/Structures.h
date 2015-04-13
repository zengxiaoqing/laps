
#define CARAC_MAX       200

//----------------------------------------
// structure des parametres a configurer
//----------------------------------------
typedef struct CFG{
  int ip2;
  int npas;
  int yyyymmddhh;
  int hh_init;
  int hh_fcst;
  int nbr_level_2keep;
  int nbr_level_avail;
  char region[10];
  char ext_file_out[CARAC_MAX];
  char ext_file_in[CARAC_MAX];
  char path_file_in[CARAC_MAX];
  char path_file_out[CARAC_MAX];
  char varname3D[CARAC_MAX/2];
  char varname2D[CARAC_MAX/4];
  float resolution;
  float thres_meso;
  float thres_tornade;
  float thres_severe_updraft[2];
  float thres_couplet_tourbillon[2][2];
} CFG;


