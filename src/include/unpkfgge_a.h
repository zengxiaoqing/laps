/*---------------------------------------------------------------------------*/
typedef struct fgge_head
  { long Q,S1,F1,T,C1,E1,M,X,S2,F2,N,C2,E2,CD,CM,KS,K,U1;
    long NW,JJ,MM,YY,GG,R,G,J,B,Z,A,U2,SN,MN,RT,DO,U3;
  } FGGE_HEAD ;
/*---------------------------------------------------------------------------*/
void open_fgge_file(char *fgge_name, int *status);
void read_fgge_record(char *obuf, int *status);
void decode_fgge_header(char *ibuf, FGGE_HEAD *header,int *status);
void decode_fgge_data(char *ibuf, FGGE_HEAD *header, float *data, int *status);
void close_fgge_file(int *status);
/*---------------------------------------------------------------------------*/
