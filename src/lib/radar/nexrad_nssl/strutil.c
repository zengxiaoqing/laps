/**************************************************************************
 *				STRUTIL.C				  *
 *									  *
 *  Version  Description					  Date	  *
 *  -------  -------------------------------------------------	--------  *
 *  1.00.00  Preliminary Program Development			07.19.91  *
 *									  *
 **************************************************************************/
#include <config.h>
#include <ctype.h>
#include <string.h>
#include <stdlib.h>
#include <malloc.h>

/*

#include <stdio.h>

main()
{
FILE *fopen();
FILE *fpin;

char buffer[2048];
char buff[2048];
char *tmp;
char *selpairfield();

int count;
int loop;

fpin=fopen("mesodata.dat","rb");

while (fgets(buffer,2048,fpin) != '\0') {
	count=pairfields(buffer,"()");

	for (loop=0; loop < count+1; loop++) {
		strcpy(buff,selpairfield(buffer,"()",loop));
		}

	}

fclose(fpin);
}
*/

/**************************************************************************
				trimright
 **************************************************************************/

char *trimright(string)
char *string;
{
char *tmp;

if (strlen(string) == 0)
	return(string);

tmp=string+(strlen(string)-1);
while (isspace(*tmp)) {
	*tmp='\0';
	tmp--;
	}

return(string);
}

/**************************************************************************
				trimcharr
 **************************************************************************/


char *trimcharr(string,letter)
char *string;
int letter;
{
char *tmp;

if (strlen(string) == 0)
	return(string);

tmp=string+(strlen(string)-1);

while (strlen(string) > 0 && *tmp == letter) {
	*tmp='\0';
	tmp--;
	}

return(string);
}

/**************************************************************************
				trimleft
 **************************************************************************/

char *trimleft(string)
char *string;
{
char *tmp;


if (strlen(string) == 0)
	return(string);

tmp=string;
while (isspace(*tmp)) {
	tmp++;
	}

strcpy(string,tmp);
return(string);
}

/**************************************************************************
				trimcharl
 **************************************************************************/

char *trimcharl(string,letter)
char *string;
int letter;
{
char *tmp;



if (strlen(string) == 0)
	return(string);

tmp=string;
while (strlen(tmp) > 0 && *tmp == letter) {
	tmp++;
	}

strcpy(string,tmp);
return(string);
}

/**************************************************************************
				trimboth
 **************************************************************************/

char *trimboth(string)
char *string;
{
char *trimright();
char *trimleft();

if (strlen(string) == 0)
	return(string);

return(trimleft(trimright(string)));
}

/**************************************************************************
				trimchar
 **************************************************************************/

char *trimchar(string,letter)
char *string;
int letter;
{
char *trimcharr();
char *trimcharl();

if (strlen(string) == 0)
	return(string);

return(trimcharl(trimcharr(string,letter),letter));
}


/**************************************************************************
				strtoupper
 **************************************************************************/

char *strtoupper(string)
char *string;
{
char *tmp;

if (strlen(string) == 0)
	return(string);

tmp=string;
while(*tmp != '\0') {
	if (islower(*tmp))
		*tmp=(*tmp-'a')+'A';

	tmp++;
	}

return(string);
}

/**************************************************************************
				strtolower
 **************************************************************************/

char *strtolower(string)
char *string;
{
char *tmp;

if (strlen(string) == 0)
	return(string);


tmp=string;
while(*tmp != '\0') {
	if (isupper(*tmp))
		*tmp=(*tmp-'A')+'a';

	tmp++;
	}

return(string);
}

/**************************************************************************
				padleft
 **************************************************************************/

char *padleft(string,length,value)
char *string;
int length;
char value;
{
char *temp;



int count;
int loop;
int size;

if ((size=strlen(string)) > length)
	return('\0');

temp=malloc(length+1);

count=length-size;
for (loop=0; loop < count; loop++)
	*(temp+loop)=value;

strcpy(temp+loop,string);

memcpy(string,temp,length);
*(string+length)='\0';

free(temp);

return(string);
}


/**************************************************************************
				padright
 **************************************************************************/

char *padright(string,length,value)
char *string;
int length;
char value;
{
char *temp;



int count;
int loop;
int size;

if ((size=strlen(string)) > length)
	return('\0');

temp=malloc(length+1);

strcpy(temp,string);

count=length-size;
for (loop=0; loop < count; loop++)
	*(temp+size+loop)=value;


memcpy(string,temp,length);
*(string+length)='\0';

free(temp);

return(string);
}

/**************************************************************************
				space
 **************************************************************************/

char *space(cnt)
int cnt;
{
char buffer[513];

int loop;

if (cnt > 512)
	return('\0');


for (loop=0; loop < cnt; loop++)
	buffer[loop]=' ';

buffer[loop]='\0';

return(buffer);
}

/**************************************************************************
				strcompress
 **************************************************************************/

char *strcompress(string)
char *string;
{
char *tmpa;
char *tmpb;
char *trimboth();


int spc;
int pnt;
int skip;

tmpa=trimboth(string);
skip=1;
tmpb=tmpa;

while(*tmpa != '\0') {
	if (ispunct(*tmpa))
		skip=3;

	else if (!isspace(*tmpa))
		skip=2;

	else if (isspace(*tmpa) && skip==0)
		strcpy(tmpa,tmpa+1);

	else
		skip--;

	if (skip > 0)
		tmpa++;
	}

return(tmpb);
}

/**************************************************************************
				strstrip
 **************************************************************************/

char *strstrip(string)
char *string;
{
char *tmpa;
char *tmpb;
char *trimboth();


int spc;
int pnt;
int skip;

tmpa=trimboth(string);
skip=0;
tmpb=tmpa;

while(*tmpa != '\0') {
	if (isspace(*tmpa)) {
		if (skip == 0) {
			*tmpa=' ';
			tmpa++;
			skip++;
			continue;
			}
		else {
			strcpy(tmpa, tmpa+1);
			continue;
			}
		}
	skip=0;
	tmpa++;
	}

return(tmpb);
}

/**************************************************************************
				strcollapse
 **************************************************************************/

char *strcollapse(string)
char *string;
{
char *tmpa;
char *tmpb;
char *trimboth();


int spc;
int pnt;
int skip;

tmpa=trimboth(string);
skip=1;
tmpb=tmpa;

while(*tmpa != '\0') {
	if (!isspace(*tmpa))
		skip=2;

	else if (isspace(*tmpa) && skip==0)
		strcpy(tmpa,tmpa+1);

	else {
		*tmpa=' ';
		skip--;
		}

	if (skip > 0)
		tmpa++;
	}

return(tmpb);
}

/**************************************************************************
				strrept
 **************************************************************************/

char *strrept(buff,str,cnt)
char *buff;
char *str;
int cnt;
{



int loop;

for (loop=0; loop < cnt; loop++)
	strcpy(buff+(loop*strlen(str)),str);

return(buff);
}


/**************************************************************************
				strjust
 **************************************************************************/


char *strjust(buff,length,str,mode)
char *buff;
int  length;
char *str;
int  mode;
{
int remainder;
int offset;
int str_length;
int middle;
int loop;
int words;
int spc;
int point;



char *space();
char *selword();
char *tmp;

tmp=(char *)malloc(strlen(str+1));
strcpy(tmp,str);

trimboth(tmp);
str_length=strlen(tmp);
remainder=length-str_length;

switch (mode) {
	case 1:	strcpy(buff,tmp);		/* left justified	  */
		strcat(buff,space(remainder));
		free(tmp);
		return(buff);

	case 2: strcpy(buff,space(remainder));	/* right justified	  */
		strcat(buff,tmp);
		free(tmp);
		return(buff);

	case 3: offset=remainder/2;		/* centered		  */
		strcpy(buff,space(offset));
		strcat(buff,tmp);
		strcat(buff,space(remainder-offset));
		free(tmp);
		return(buff);

	case 4:	offset=remainder/((words=strwords(tmp))-1);  /* word fill */
		strcpy(buff,selword(tmp,1));
		strcat(buff,space(offset));
		for (loop=1; loop < words-1; loop++) {
			strcat(buff,selword(tmp,loop+1));
			strcat(buff,space(offset));
			}
		strcat(buff,selword(tmp,words));
		free(tmp);
		return(buff);


	case 5: if (remainder < (str_length*3)) {	/* do word fill ? */
			offset=remainder/((words=strwords(tmp))-1);
			strcpy(buff,selword(tmp,1));
			strcat(buff,space(offset));
			for (loop=1; loop < words-1; loop++) {
				strcat(buff,selword(tmp,loop+1));
				strcat(buff,space(offset));
				}
			strcat(buff,selword(tmp,words));
			free(tmp);
			return(buff);
			}

		else {				/* character fill	  */
			str_length=strlen(strcollapse(tmp));
			remainder=length-str_length;
			offset=remainder/(str_length+(((words=strwords(tmp))-1)*2));
			spc=(remainder%(str_length+(words*2)))/(words-1);

			buff[0]=tmp[0];
			strcpy(buff+1,space(offset));
			for(loop=1; loop < str_length-1; loop++) {
				point=strlen(buff);
				buff[point]=str[loop];
				strcpy(buff+point+1,space(offset));
				if (isspace(str[loop]))
					strcat(buff,space(spc));
				}
			strcat(buff,tmp+loop);
			free(tmp);
			return(buff);
			}
	}
free(tmp);
return('\0');
}

/**************************************************************************
				strwords
 **************************************************************************/

int strwords(str)
char *str;
{
char *buff;

char *strcollapse();


int count;
int loop;

if((buff=(char *)malloc(strlen(str)+1)) == '\0')
	return(-1);

count=1;

strcpy(buff,str);
strcollapse(buff);

for (loop=0; loop < strlen(buff); loop++) {
	if (buff[loop] == ' ')
		count++;
        }

free(buff);
return(count);
}


/**************************************************************************
				selword
 **************************************************************************/

char *selword(str,word)
char *str;
int   word;
{
char wd_buff[256];

char *buff;
char *wd;

char *strcollapse();


int count;
int loop;
int words;
int loopend;

if (word==0)
	return(str);

if (word > (words=strwords(str)))
	return('\0');


if((buff=(char *)malloc(strlen(str)+1)) == '\0')
	return('\0');


strcpy(buff,str);
strcollapse(buff);

wd=buff;
loop=0;
loopend=strlen(buff);

for (count=0; count < words; count++) {

	for (; loop < loopend; loop++) {
		if (buff[loop] == ' ')
			break;
		}

	if (word-1 == count) {
		buff[loop] = '\0';
		strcpy(wd_buff,wd);
                free(buff);
		return(wd_buff);
		}
	wd=buff+loop+1;
	loop++;
	}

free(buff);
return('\0');
}

/**************************************************************************
				strins
 **************************************************************************/

char *strins(dest,source,offset)
char *dest;
int source;
int offset;
{
char *temp;



temp=(char *)malloc(strlen(dest)+2);

if (offset > strlen(dest))
	return(dest);

strcpy(temp,dest);

*(temp+offset)=source;
strcpy(temp+offset+1,dest+offset);
strcpy(dest,temp);
return(dest);
}

/**************************************************************************
				strinsert
 **************************************************************************/

char *strinsert(dest,source,offset)
char *dest;
char *source;
int offset;
{
char *temp;



temp=(char *)malloc(strlen(dest)+strlen(source)+1);

if (offset > strlen(dest))
	return(dest);

strcpy(temp,dest);
strcpy(temp+offset,source);
strcat(temp,dest+offset);
strcpy(dest,temp);
return(dest);
}


/**************************************************************************
				strdel
 **************************************************************************/

char *strdel(source,sub,occur)
char *source;
char *sub;
int  occur;
{
char *tmp_a;
char *tmp_b;



int count;
int loop;

if (occur == 0) {					 /* all occurance */
        tmp_a=source;
	while((tmp_b=strchr(tmp_a,*sub)) != '\0') {
		if (strncmp(tmp_b,sub,strlen(sub)) == 0)
                	strcpy(tmp_b,tmp_b+strlen(sub));
                }
	return(source);
        }

else if (occur == -1) {					/* last occurance */
	tmp_a=source;
        while ((tmp_b=strchr(tmp_a,*sub)) != '\0') {
        	if (strncmp(tmp_b,sub,strlen(sub)) == 0)
                	tmp_a=tmp_b+1;
                }
        if (tmp_a == source)
        	return(source);
        tmp_a--;
        strcpy(tmp_a,tmp_a+strlen(sub));
        return(source);
        }

else {							/* nth occurance */
	count=0;
	tmp_a=source;
        while ((tmp_b=strchr(tmp_a,*sub)) != '\0') {
        	count++;
                if (count == occur) {
                	strcpy(tmp_b,tmp_b+strlen(sub));
                        return(source);
                        }
                else {
                	tmp_a=tmp_b+1;
                        }
                }
        return(source);
        }

}

/**************************************************************************
				strdelete
 **************************************************************************/

char *strdelete(buff,offset,count)
char *buff;
int offset;
int count;
{
int size;
char *temp;

offset--;
size=strlen(buff);

if (offset >= size)
	return(buff);

if (offset+count > size)
        count=size-offset;

temp=buff;
temp+=offset;
*temp='\0';

strcpy(temp,temp+count);
return(buff);
}

/**************************************************************************
				buffdelete
 **************************************************************************/

void *buffdelete(buff,size,offset,count)
void *buff;
int   size;
int   offset;
int   count;
{
char *temp;

if (offset > size)
	return(buff);

if (offset+count > size)
	return(buff);

temp=(char *)buff;
temp=temp+offset;

memmove((void *)temp,(void *)(temp+count),size-(offset+count));

temp=temp+count;
initbuff(temp,size-(offset+count),'\0');
return(buff);
}

/**************************************************************************
				buffinsert
 **************************************************************************/

void *buffinsert(dest,source,size,offset,count)
void *dest;
void *source;
int size;
int offset;
int count;
{
char *temp;
char *tmp;

if ((tmp=(char *)malloc(size-offset)) == '\0')
	return('\0');

temp=(char *)dest;
temp=temp+offset;
memmove((void *)tmp,(void *)temp,size-offset);
memmove((void *)temp,source,count);
memmove((void *)(temp+count),tmp,size-offset);

free(tmp);

return(dest);
}

/**************************************************************************
				buff_format
 **************************************************************************/

int buff_format(array,string,width)
char *array[];
char *string;
int   width;
{
char *buff;
char temp[256];
char *tmpa;
char *tmpb;
char  tmpc;

char  *strcollapse();
char **strtable();

int   lines;
int   loop;
int   length;

strcollapse(string);
if ((length=strlen(string)) == 0)
	return(0);

loop=0;
lines=0;
tmpa=string;

while (loop < length) {
	if ((loop+width) > length)
		tmpb=string+length;
	else
		tmpb=tmpa+width;

        while(!isspace((int)*tmpb) && *tmpb != '\0')
        	tmpb--;

        *tmpb='\0';
	array[lines]=tmpa;

	loop+=strlen(tmpa)+1;
        lines++;
	tmpa=tmpb+1;
        }

return(lines);
}

/**************************************************************************
				textfields
 **************************************************************************/

int textfields(buffer,delimiter)
char *buffer;
int delimiter;
{
int length;
int fields;
int loop;

char *temp;

if (strlen(buffer) == 0)
	return(0);

if (delimiter == ' ') {	       /* if the delimiter is the space character */
	temp=(char *)malloc(strlen(buffer)+1);	/* collapse the string	  */
        strcpy(temp,buffer);
        strcollapse(temp);
        }
else
	temp=buffer;

length=strlen(temp);
fields=1;

for (loop=0; loop < length; loop++) {
	if ((int)(temp[loop]) == delimiter  && (loop == 0 || temp[loop-1] != '\\'))
        	fields++;
	}

return(fields);
}

/**************************************************************************
				selfield
 **************************************************************************/

char *selfield(str,word,delimiter)
char *str;
int   word;
int   delimiter;
{
char wd_buff[256];

char *buff;
char *wd;

char *strcollapse();


int count;
int loop;
int words;
int loopend;

if (word==0)
	return(str);

if (word > (words=textfields(str,delimiter)))
	return('\0');


if((buff=(char *)malloc(strlen(str)+1)) == '\0')
	return('\0');


strcpy(buff,str);

if (delimiter == ' ')
	strcollapse(buff);

wd=buff;
loop=0;
loopend=strlen(buff);

for (count=0; count < words; count++) {

	for (; loop < loopend; loop++) {
		if ((int)(buff[loop]) == delimiter && (loop == 0 || buff[loop-1] != '\\'))
			break;
		}

	if (word-1 == count) {
		buff[loop] = '\0';
		strcpy(wd_buff,wd);
                free(buff);
		return(wd_buff);
		}
	wd=buff+loop+1;
	loop++;
	}

free(buff);
return('\0');
}


/**************************************************************************
				textlines
 **************************************************************************/

int textlines(buffer)
char *buffer;
{
int length;
int lines;
int loop;


lines=0;
length=strlen(buffer);

if (length <= 0)
	return(0);

for (loop=0; loop < length; loop++) {
	if (loop == 0 && (buffer[loop] == '\r' || buffer[loop] == '\n'))
        	lines++;

        if (loop != 0 && buffer[loop] == '\r' && buffer[loop-1] != '\n')
        	lines++;

        if (loop != 0 && buffer[loop] == '\n' && buffer[loop-1] != '\r')
        	lines++;
        }


if (lines == 0 || (buffer[length-1] != '\r' && buffer[length-1] != '\n'))
	lines++;


return(lines);
}

/**************************************************************************
				selline
 **************************************************************************/

char *selline(buffer,line)
char *buffer;
int line;
{
char buff[1024];
int length;
int lines;
int loop;
int offset;

if ((length=strlen(buffer)) < 0)
	return('\0');

if (length == 0) {
	buff[0] = '\0';
        return(buff);
	}

lines=0;
offset=0;

for (loop=0; loop < length; loop++) {
	if (loop == 0 && (buffer[loop] == '\r' || buffer[loop] == '\n'))
        	lines++;

        if (loop != 0 && buffer[loop] == '\r' && buffer[loop-1] != '\n')
        	lines++;

        if (loop != 0 && buffer[loop] == '\n' && buffer[loop-1] != '\r')
        	lines++;

        if (line == lines && buffer[loop] != '\r' && buffer[loop] != '\n') {
        	buff[offset]=buffer[loop];
                offset++;
                }

        if (lines > line) {
        	buff[offset]='\0';
        	return(buff);
                }
        }

if (offset == 0)
	return('\0');
else {
	buff[offset]='\0';
        return(buff);
        }
}

/**************************************************************************
				fitedit
 **************************************************************************/

char *fitedit(string,edit)
char *string;
char *edit;
{
char *tempstr;
char *tempedit;
char *editarray;
char *tmp_a;
char *tmp_b;
char *tmp_c;

char *_fitedit();

if ((tempstr=malloc(strlen(string))) == '\0')
	return('\0');

if ((tempedit=malloc(strlen(edit))) == '\0') {
	free(tempstr);
        return('\0');
        }

if ((editarray=malloc(strlen(edit))) == '\0') {
	free(tempstr);
        free(tempedit);
        return('\0');
        }

strcpy(tempstr,string);
strcpy(tempedit,edit);

if (*(tmp_a=trimchar(tempedit,'*')) == '\0')
	return(string);

strcpy(tempedit,edit);
tmp_a=tempedit;

if ((tmp_b=strchr(tmp_a,'|')) != '\0') {
	while((tmp_b=strchr(tmp_a,'|')) != '\0') {
		*tmp_b='\0';
                strcpy(editarray,tmp_a);

                if ((tmp_c=_fitedit(tempstr,editarray)) != '\0') {
                	strcpy(string,tmp_c);
                        free(tempstr);
                        free(tempedit);
                        free(editarray);
                        return(string);
                        }

                tmp_a=tmp_b+1;
                }
        
        return('\0');
        }

if ((tmp_c=_fitedit(tempstr,tmp_a)) != '\0') {
	strcpy(string,tmp_c);
        free(tempstr);
        free(tempedit);
        free(editarray);
        return(string);
        }

return('\0');
}

/**************************************************************************
				_fitedit
 **************************************************************************/

char *_fitedit(string,edit)
char *string;
char *edit;
{
int loop_str;
int loop_edt;
int loop_buf;

int space_flag;
int dec_flag;
int sign_flag;

char buff[256];
char default_chr;

loop_str=0;
loop_buf=0;

space_flag=0;
dec_flag=0;
sign_flag=0;

for (loop_edt=0; *(edit+loop_edt) != '\0'; loop_edt++) {
	switch(*(edit+loop_edt)) {
        	case '1':
                case '2':
                case '3':
                case '4':
                case '5':
                case '6':
                case '7':
                case '8':
                case '9':  if (*(string+loop_str) >= '0' && *(string+loop_str) <= *(edit+loop_edt)) {
                		buff[loop_buf]=(*(string+loop_str));
                                loop_buf++;
                                loop_str++;
                                }

                	   else {
                           	return('\0');
                                }
                	   break;

		case 'A':  buff[loop_buf]=toupper(*(string+loop_str));
			   loop_buf++;
			   loop_str++;
			   break;
			   
		case 'a':  buff[loop_buf]=tolower(*(string+loop_str));
			   loop_buf++;
			   loop_str++;
			   break;
		
		case '+':  if (*(string+loop_str) == ' ' && space_flag == 0) {
                		buff[loop_buf]=(*(string+loop_str));
                                loop_buf++;
                                loop_str++;
                                break;
                                }

                           if ((*(string+loop_str) == '+' || *(string+loop_str) == '-') && sign_flag == 0) {
                           	buff[loop_str]=(*(string+loop_str));
                                loop_buf++;
                                loop_str++;
                                space_flag=1;
                                sign_flag=1;
                                break;
                                }

                           if (*(string+loop_str) == '.' && dec_flag == 0) {
                           	buff[loop_buf]=(*(string+loop_str));
                           	loop_buf++;
                                loop_str++;
                                space_flag=1;
                                sign_flag=1;
                                dec_flag=1;
                                break;
                                }

                           if (isdigit(*(string+loop_str))) {
                           	buff[loop_buf]=(*(string+loop_str));
                                loop_buf++;
                                loop_str++;
                                space_flag=1;
                                sign_flag=1;
                                break;
                                }
			   
                           return('\0');
                           

                case '?':  if (*(string+loop_str) == ' ' && space_flag == 0) {
                		buff[loop_buf]=(*(string+loop_str));
                                loop_buf++;
                                loop_str++;
                                break;
                                }

		case '#':  if (isdigit(*(string+loop_str)) || (*(string+loop_str) == '.' && dec_flag == 0))  {
                                space_flag=1;

                                if (*(string+loop_str) == '.')
                                	dec_flag=1;

                                buff[loop_buf]=(*(string+loop_str));
				loop_buf++;
				loop_str++;                 
                                }

                  	   else {
                           	return('\0');
                                }
                           break;

                case '^':  if (loop_buf == 0 || buff[loop_buf-1] == ' ')
                		buff[loop_buf]=toupper(*(string+loop_str));
			   else
                           	buff[loop_buf]=tolower(*(string+loop_str));
                           loop_buf++;
                           loop_str++;
                           break;

		case '*':  buff[loop_buf]=(*(string+loop_str));
			   loop_buf++;
                           loop_str++;
                           break;

                case '$':  if (isalpha((int)(*(string+loop_str))) != 0) {
				buff[loop_buf]=(*(string+loop_str));
                	   	loop_buf++;
                                loop_str++;
                                }
                           else
                           	return('\0');
                           break;

                case '[':  default_chr=(*(edit+loop_edt+1));
			   while ((*(string+loop_str) != *(edit+loop_edt)) && (*(edit+loop_edt) != ']' && *(edit+loop_edt) != '\0'))
                		loop_edt++;

		           if (*(edit+loop_edt) == '\0' || *(edit+loop_edt) == ']') {
                                buff[loop_buf]=default_chr;
				loop_buf++;
                                }

			   else {
			   	buff[loop_buf]=(*(string+loop_str));
				loop_buf++;
				loop_str++;

				while (*(edit+loop_edt) != ']' && *(edit+loop_edt+1) != '\0')
                                	loop_edt++;
                                }
                           break;

		case '\\': loop_edt++;

                default:   if (*(string+loop_str) == *(edit+loop_edt)) {
                		buff[loop_buf]=(*(string+loop_str));
                                loop_buf++;
                                loop_str++;
                                }

		           else {
                           	buff[loop_buf]=(*(edit+loop_edt));
                                loop_buf++;
                                }
                }
/*
        if (*(string+loop_str) == '\0' && *(edit+loop_edt+1) != '\0') {
        	return('\0');
                }
*/
        }


if (*(edit+loop_edt) == '\0' && *(string+loop_str) == '\0') {
	buff[loop_buf]='\0';
        return(buff);
        }

return('\0');
}

/**************************************************************************
				checkedit
 **************************************************************************/

int checkedit(string,edit)
char *string;
char *edit;
{
char *tempstr;
char *tempedit;
char *editarray;
char *tmp_a;
char *tmp_b;
char *tmp_c;

char *_fitedit();

int test;

if ((tempstr=malloc(strlen(string))) == '\0')
	return(-1);

if ((tempedit=malloc(strlen(edit))) == '\0') {
	free(tempstr);
        return(-1);
        }

if ((editarray=malloc(strlen(edit))) == '\0') {
	free(tempstr);
        free(tempedit);
        return(-1);
        }

strcpy(tempstr,string);
strcpy(tempedit,edit);

if (*(tmp_a=trimchar(tempedit,'*')) == '\0')
	return(0);

strcpy(tempedit,edit);
tmp_a=tempedit;

if ((tmp_b=strchr(tmp_a,'|')) != '\0') {
	while((tmp_b=strchr(tmp_a,'|')) != '\0') {
		*tmp_b='\0';
                tmp_b++;
                strcpy(editarray,tmp_a);

                if ((test=_checkedit(tempstr,editarray)) == 0) {
                        free(tempstr);
                        free(tempedit);
                        free(editarray);
                        return(0);
                        }

                tmp_a=tmp_b;
                }
	}

if ((test=_checkedit(tempstr,tmp_a)) != '\0') {
	strcpy(string,tmp_c);
        free(tempstr);
        free(tempedit);
        free(editarray);
        return(0);
        }

return(1);
}

/**************************************************************************
				_checkedit
 **************************************************************************/

int _checkedit(string,edit)
char *string;
char *edit;
{
int loop_str;
int loop_edt;
int loop_buf;

int space_flag;
int dec_flag;
int sign_flag;

char buff[256];
char default_chr;

loop_str=0;
loop_buf=0;

space_flag=0;
dec_flag=0;
sign_flag=0;

for (loop_edt=0; *(edit+loop_edt) != '\0'; loop_edt++) {
	switch(*(edit+loop_edt)) {
        	case '1':
                case '2':
                case '3':
                case '4':
                case '5':
                case '6':
                case '7':
                case '8':
                case '9':  if (*(string+loop_str) >= '0' && *(string+loop_str) <= *(edit+loop_edt)) {
                		buff[loop_buf]=(*(string+loop_str));
                                loop_buf++;
                                loop_str++;
                                }

                	   else {
                           	return(1);
                                }
                	   break;

		case 'A':  buff[loop_buf]=toupper(*(string+loop_str));
			   loop_buf++;
			   loop_str++;
			   break;
			   
		case 'a':  buff[loop_buf]=tolower(*(string+loop_str));
			   loop_buf++;
			   loop_str++;
			   break;
		
		case '+':  if (*(string+loop_str) == ' ' && space_flag == 0) {
                		buff[loop_buf]=(*(string+loop_str));
                                loop_buf++;
                                loop_str++;
                                break;
                                }

                           if ((*(string+loop_str) == '+' || *(string+loop_str) == '-') && sign_flag == 0) {
                           	buff[loop_str]=(*(string+loop_str));
                                loop_buf++;
                                loop_str++;
                                space_flag=1;
                                sign_flag=1;
                                break;
                                }

                           if (*(string+loop_str) == '.' && dec_flag == 0) {
                           	buff[loop_buf]=(*(string+loop_str));
                           	loop_buf++;
                                loop_str++;
                                space_flag=1;
                                sign_flag=1;
                                dec_flag=1;
                                break;
                                }

                           if (isdigit(*(string+loop_str))) {
                           	buff[loop_buf]=(*(string+loop_str));
                                loop_buf++;
                                loop_str++;
                                space_flag=1;
                                sign_flag=1;
                                break;
                                }
			   
                           return(1);


                case '?':  if (*(string+loop_str) == ' ' && space_flag == 0) {
                		buff[loop_buf]=(*(string+loop_str));
                                loop_buf++;
                                loop_str++;
                                break;
                                }

		case '#':  if (isdigit(*(string+loop_str)) || (*(string+loop_str) == '.' && dec_flag == 0))  {
                                space_flag=1;

                                if (*(string+loop_str) == '.')
                                	dec_flag=1;

                                buff[loop_buf]=(*(string+loop_str));
				loop_buf++;
				loop_str++;                 
                                }

                  	   else {
                           	return(1);
                                }
                           break;

                case '^':  if (loop_buf == 0 || buff[loop_buf-1] == ' ')
                		buff[loop_buf]=toupper(*(string+loop_str));
			   else
                           	buff[loop_buf]=tolower(*(string+loop_str));
                           loop_buf++;
                           loop_str++;
                           break;

		case '*':  buff[loop_buf]=(*(string+loop_str));
			   loop_buf++;
                           loop_str++;
                           break;

                case '$':  if (isalpha((int)(*(string+loop_str))) != 0) {
				buff[loop_buf]=(*(string+loop_str));
                	   	loop_buf++;
                                loop_str++;
                                }
                           else
                           	return(1);
                           break;

                case '[':  default_chr=(*(edit+loop_edt+1));
			   while ((*(string+loop_str) != *(edit+loop_edt)) && (*(edit+loop_edt) != ']' && *(edit+loop_edt) != '\0'))
                		loop_edt++;

		           if (*(edit+loop_edt) == '\0' || *(edit+loop_edt) == ']') {
                                return(1);
                                }

			   else {
			   	buff[loop_buf]=(*(string+loop_str));
				loop_buf++;
				loop_str++;

				while (*(edit+loop_edt) != ']' && *(edit+loop_edt+1) != '\0')
                                	loop_edt++;
                                }
                           break;

		case '\\': loop_edt++;

                default:   if (*(string+loop_str) == *(edit+loop_edt)) {
                		buff[loop_buf]=(*(string+loop_str));
                                loop_buf++;
                                loop_str++;
                                }

		           else {
                           	return(1);
                                }
		}

        if (*(string+loop_str) == '\0' && *(edit+loop_edt+1) != '\0')
        	return(1);

        }


if (*(edit+loop_edt) == '\0' && *(string+loop_str) == '\0') {
	buff[loop_buf]='\0';
        return(0);
        }

return(1);
}


/**************************************************************************
				strdecode
 **************************************************************************/

char *strdecode(string)
char *string;
{
int loop;
int length;
int tm_loop;

char *temp;
char *tmp;

length=strlen(string);			/* get the length of the string	  */
								  
if ((temp=(char *)malloc(length+1)) == '\0')	/* allocate working space */
	return('\0');				/* error allocating space */

strcpy(temp,string);			/* copy to working space	  */

tmp=string;				/* setup point for decoded string */

for (loop=0; loop < length; loop++) {	/* loop through the string	  */
	if (temp[loop] == '\\')	{	/* special character		  */
        	loop++;
                switch(temp[loop]) {
                	case 'b': *tmp='\b';	/* match backspace	  */
                        	  break;

			case 'E':
                        case 'e': *tmp='\033';	/* match escape		  */
                        	  break;

                        case 'f': *tmp='\f';	/* match formfeed	  */
                        	  break;

                        case 'n': *tmp='\n';	/* match newline	  */
                        	  break;

                        case 'r': *tmp='\r';	/* match carriage return  */
                        	  break;

                        case 't': *tmp='\t';	/* match tab		  */
                        	  break;

                        case 'x': *tmp='\0';  	/* decode hex charcter    */
				  for (tm_loop=1; tm_loop < 3; tm_loop++) {
                                        (*tmp)<<=4;
                        	  	switch(temp[loop+tm_loop]) {
                                        	case '0':
                                                case '1':
                                                case '2':
                                                case '3':
                                                case '4':
                                                case '5':
                                                case '6':
                                                case '7':
                                                case '8':
                                                case '9': (*tmp)|=(temp[loop+tm_loop]-'0');
                                                	  break;
                                                case 'a':
                                                case 'b':
                                                case 'c':
                                                case 'd':
                                                case 'e':
                                                case 'f': (*tmp)|=((temp[loop+tm_loop]-'a')+10);
                                                	  break;

                                                case 'A':
                                                case 'B':
                                                case 'C':
                                                case 'D':
                                                case 'E':
                                                case 'F': (*tmp)|=((temp[loop+tm_loop]-'A')+10);
                                                	  break;

                                                default:  strcpy(string,temp);
                                                	  free(temp);
                                                          return('\0');
                                                }
                                        }
                                  loop+=2;
                                  break;

                        case '0':		/* octal character	  */
                        case '1':
                        case '2':
                        case '3':
                        case '4':
                        case '5':
                        case '6':
                        case '7': *tmp='\0';
                        	  for (tm_loop=0; tm_loop < 3; tm_loop++) {
                                  	(*tmp)<<=3;
                                  	if (temp[loop+tm_loop] < '0' || temp[loop+tm_loop] > '7') {
                                        	strcpy(string,temp);
                                                free(temp);
                                                return('\0');
                                                }
                                        (*tmp)|=(temp[loop+tm_loop]-'0');
                                        }
                                  loop+=2;
                                  break;

                        default:  *tmp=temp[loop];    /* match literal	  */
				  break;

                        }
                }

	else
        	*tmp=temp[loop];
        tmp++;
        }

*tmp='\0';
free(temp);
return(string);
}

/**************************************************************************
				pairfields
 **************************************************************************/

int pairfields(string,delimiter)
char *string;
char *delimiter;
{
int left;
int right;
char buffer[8];

if (strlen(delimiter) != 2)
	return(-1);

strcpy(buffer,delimiter);
buffer[1]='\0';

left=stroccur(string,buffer);
strcpy(buffer,delimiter+1);
right=stroccur(string,buffer);

if (right < left)
	return(right);
return(left);
}

/**************************************************************************
				selpairfield
 **************************************************************************/

char *selpairfield(string,delimiter,occurance)
char *string;
char *delimiter;
int occurance;
{
char result[1024];
char tmp_string[2048];
char *left;
char *right;
char *seloccur();
int count;
int left_count;
int right_count;

char tmp_left[8];

strcpy(tmp_string,string);

if (occurance == 0)
	return(tmp_string);

if ((right=seloccur(tmp_string,delimiter+1,occurance)) == '\0') {
	return('\0');
	}

*(right+1)='\0';
strcpy(tmp_left,delimiter);
*(tmp_left+1)='\0';

left=right-1;
left_count=0;
right_count=1;

while (left > tmp_string) {
	if (*left == *delimiter)
		left_count++;
	if (*left == *(delimiter+1))
		right_count++;

	if (left_count == right_count) {
		strcpy(result,left+1);
		*(result+strlen(result)-1)='\0';
		return(result);
		}
	left--;
	}

return('\0');
}

/**************************************************************************
				stroccur
 **************************************************************************/

int stroccur(string,substr)
char *string;
char *substr;
{
int count;
int result;


char *tmp;

tmp=string;
count=0;

while ((tmp=strchr(tmp,*substr)) != '\0') {
	if (strncmp(tmp,substr,strlen(substr)) == 0) {
		count++;
		tmp+=strlen(substr);
		}
	else {
		tmp++;
		}
	}

return(count);
}

/**************************************************************************
				seloccur
 **************************************************************************/

char *seloccur(string,substr,occurance)
char *string;
char *substr;
int occurance;
{
char *tmp;

int count;

count=stroccur(string,substr);
if (occurance < 1 || occurance > count)
	return('\0');

tmp=string;
count=0;

while ((tmp=strchr(tmp,*substr)) != '\0') {
	if (strncmp(tmp,substr,strlen(substr)) == 0) {
		count++;
		if (count == occurance)
			return(tmp);
		tmp+=strlen(substr);
		}
	else
		tmp++;
	}

return('\0');
}

/**************************************************************************
				strreplace
 **************************************************************************/

char *strreplace(buf,sub,rep,occurance)
char *buf;
char *sub;
char *rep;
int occurance;
{
char buffer[4096];
int loop;
int offset;
int count;
int event;

if (strlen(buf) == 0)
	return('\0');

if (strlen(sub) == 0)
	return('\0');

if ((count=strcount(buf,sub)) == 0) 
	return(buf);

event=0;
offset=0;
for (loop=0; loop < (strlen(buf)-strlen(sub)); loop++) {
	if (*(buf+loop) == *sub) {
		if (strncmp(buf+loop,sub,strlen(sub)) == 0) {
			event++;
			if (occurance == 0 || event == occurance || ((count+occurance)+1) == event) {
				strcpy(buffer+offset,rep);
				offset+=strlen(rep);
				loop+=strlen(sub);
				if (occurance == 0){
					loop--;
					continue;
					}
				strcpy(buffer+offset,buf+loop);
				loop=strlen(buf);
				offset=strlen(buffer);
				break;
				}
			else {
				*(buffer+offset)=(*(buf+loop));
				offset++;
				continue;
				}
			}
		else {
			*(buffer+offset)=(*(buf+loop));
			offset++;
			continue;
			}
		}
	else {
		*(buffer+offset)=(*(buf+loop));
		offset++;
		continue;
		}
	}
*(buffer+offset)='\0';
strcat(buffer,buf+loop);
strcpy(buf,buffer);
return(buf);
}


/***************************************************************************
				strchrcount
 ***************************************************************************

int strchrcnt(buff,letter)
char *buff;
int letter;
{
char *pnt;
int count;

count=0;

for (pnt=buff;pnt != '\0' && *pnt != '\0';pnt++) {
	if ((int)(*pnt) == letter)
		count++;
	}

return(count);
}

/***************************************************************************
				strcount
 ***************************************************************************/

int strcount(buff,sub)
char *buff;
char *sub;
{
return(stroccur(buff,sub));
}

/**************************************************************************
				strchrcpy
 **************************************************************************/

char *strchrcpy(dest,source,limit)
char *dest;
char *source;
int limit;
{
char *tmp_d;
char *tmp_s;

tmp_d=dest;
tmp_s=source;

if (tmp_d == '\0')
	return(0);

if (tmp_s == '\0')
	return('\0');

while(*tmp_s != limit && *tmp_s != '\0') {
	*tmp_d=(*tmp_s);
	tmp_d++;
	tmp_s++;
	}

return(dest);
}

/**************************************************************************
				strchrcat
 **************************************************************************/

char *strchrcat(dest,source,limit)
char *dest;
char *source;
int limit;
{
char *tmp_d;
char *tmp_s;
char *strchrcpy();

tmp_d=dest;
tmp_s=source;

if (tmp_d == '\0')
	return('\0');

tmp_d=dest+strlen(dest);

return(strchrcpy(tmp_d,tmp_s,limit));
}

/**************************************************************************
				strreverse
 **************************************************************************/

char *strreverse(buffer)
char *buffer;
{
int loop;
int length;

char *temp_a;
char *temp_b;
char temp_c;

length=strlen(buffer);
temp_a=buffer;
temp_b=buffer+strlen(buffer)-1;

for (loop=0; loop < length/2; loop++) {
	temp_a=buffer+loop;
	temp_b=buffer+length-(loop+1);
	temp_c=(*temp_a);
	*temp_a=(*temp_b);
	*temp_b=temp_c;
	}

return(buffer);
}
