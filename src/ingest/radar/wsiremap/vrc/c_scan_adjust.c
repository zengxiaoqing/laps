#ifdef FORTRANUNDERSCORE
#define c_scan_adjust c_scan_adjust_
#endif
#ifdef FORTRANDOUBLEUNDERSCORE
#define c_scan_adjust c_scan_adjust__
#endif
#ifdef FORTRANCAPS
#define c_scan_adjust C_SCAN_ADJUST
#endif

#ifdef __STDC__
void c_scan_adjust(unsigned char *ByteArray, long *lines, long *elem, long *miss)
#else
void c_scan_adjust(ByteArray, lines, elem, miss)
unsigned char *ByteArray; 
long *lines; 
long *elem;
long *miss;
#endif
{
	int i,j, tot_elem;
        unsigned char *ptr;

        ptr = ByteArray;
        tot_elem = ((*lines)*(*elem));
        for (i = 0; i < tot_elem; i++) {
          switch ((int)*ptr)
          {
            case 0:
              break;
            case 1:
              break;
            case 2:
              break;
            case 3:
              break;
            case 4:
              break;
            case 5:
              break;
            case 6:
              break;
            case 7:
              break;
            case 8:
              break;
            case 9:
              break;
            case 10:
              break;
            case 11:
              break;
            case 12:
              break;
            case 13:
              break;
            case 14:
              break;
            case 15:
              break;
            case 16:
              *ptr = (unsigned char)0;
              break;
            case 17:
              *ptr = (unsigned char)1;
              break;
            case 18:
              *ptr = (unsigned char)2;
              break;
            case 19:
              *ptr = (unsigned char)3;
              break;
            case 20:
              *ptr = (unsigned char)4;
              break;
            case 21:
              *ptr = (unsigned char)5;
              break;
            case 22:
              *ptr = (unsigned char)6;
              break;
            case 23:
              *ptr = (unsigned char)7;
              break;
            case 24:
              *ptr = (unsigned char)8;
              break;
            case 25:
              *ptr = (unsigned char)9;
              break;
            case 26:
              *ptr = (unsigned char)10;
              break;
            case 27:
              *ptr = (unsigned char)11;
              break;
            case 28:
              *ptr = (unsigned char)12;
              break;
            case 29:
              *ptr = (unsigned char)13;
              break;
            case 30:
              *ptr = (unsigned char)14;
              break;
            case 31:
              *ptr = (unsigned char)15;
              break;
            case 32:
              *ptr = (unsigned char)0;
              break;
            case 33:
              *ptr = (unsigned char)1;
              break;
            case 34:
              *ptr = (unsigned char)2;
              break;
            case 35:
              *ptr = (unsigned char)3;
              break;
            case 36:
              *ptr = (unsigned char)4;
              break;
            case 37:
              *ptr = (unsigned char)5;
              break;
            case 38:
              *ptr = (unsigned char)6;
              break;
            case 39:
              *ptr = (unsigned char)7;
              break;
            case 40:
              *ptr = (unsigned char)8;
              break;
            case 41:
              *ptr = (unsigned char)9;
              break;
            case 42:
              *ptr = (unsigned char)10;
              break;
            case 43:
              *ptr = (unsigned char)11;
              break;
            case 44:
              *ptr = (unsigned char)12;
              break;
            case 45:
              *ptr = (unsigned char)13;
              break;
            case 46:
              *ptr = (unsigned char)14;
              break;
            case 47:
              *ptr = (unsigned char)15;
              break;
            case 48:
              *ptr = (unsigned char)0;
              break;
            case 49:
              *ptr = (unsigned char)1;
              break;
            case 50:
              *ptr = (unsigned char)2;
              break;
            case 51:
              *ptr = (unsigned char)3;
              break;
            case 52:
              *ptr = (unsigned char)4;
              break;
            case 53:
              *ptr = (unsigned char)5;
              break;
            case 54:
              *ptr = (unsigned char)6;
              break;
            case 55:
              *ptr = (unsigned char)7;
              break;
            case 56:
              *ptr = (unsigned char)8;
              break;
            case 57:
              *ptr = (unsigned char)9;
              break;
            case 58:
              *ptr = (unsigned char)10;
              break;
            case 59:
              *ptr = (unsigned char)11;
              break;
            case 60:
              *ptr = (unsigned char)12;
              break;
            case 61:
              *ptr = (unsigned char)13;
              break;
            case 62:
              *ptr = (unsigned char)14;
              break;
            case 63:
              *ptr = (unsigned char)15;
              break;
            case 64:
              *ptr = (unsigned char)0;
              break;
            default:
              *ptr = (unsigned char)miss;
              break;
          }
          ptr++;
        }




}
