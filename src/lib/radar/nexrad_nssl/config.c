/********************************************************************************************************************************
 *							openconfig.c								*
 *																*
 * Version  Description												  Date	  Int	*
 * -------  --------------------------------------------------------------------------------------------------	--------  ---	*
 * 1.00.00  Preliminary Program Development									04.03.90  DTR	*
 *																*
 *																*
 * Processes the named file for configuration varibles.										*
 *																*
 ********************************************************************************************************************************/
#include <config.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
char *strtoupper(char *string);
/*
main()
{
char **table;
char **openconfig();
char **initconfig();
char *getvalue();
char *temp;

table=initconfig();
table=openconfig("disk$user:[lsdb]lsdbconf.dat",table);

temp=getvalue("LSDB_TAPE_DEVICE",table);
}
*/
/********************************************************************************************************************************
 *							openconfig()								*
 ********************************************************************************************************************************/

char **openconfig(filename,config_table)
char *filename;							/* File name containing the configuration information		*/
char **config_table;						/* Table containing the configuration information		*/
{
FILE *conf_file;						/* File pointer to configuration file				*/
FILE *fopen();

char  config_ent[256];						/* buffer for holding the configuration entry 			*/
char  config_val[256];						/* buffer for holding the configuration value			*/
char  config_var[256];						/* buffer for holding the configuration name			*/
char  file_buffer[256];						/* buffer for holding information read from the file		*/
char *temp;


char *trimboth();						/* delete the white space before and after the ascii string	*/



char **setconfig();						/* function to maintain the configuration table			*/

if ((conf_file = fopen(filename,"r")) == '\0')			/* check to see if the file can be read				*/
	return(config_table);						/* the configuration file cannot be read			*/

while (fgets(file_buffer,256,conf_file) != '\0') {		/* read the file until EOF or error, no distiction made		*/
	if ((temp=strchr(file_buffer,'!')) != '\0')		/* test to see if this line has comments			*/
		*temp='\0';					/* delete the comments						*/

	if (strlen(file_buffer) < 3)
		continue;

	if ((temp=strchr(file_buffer,'=')) == '\0')		/* test to see if the line is in the correct format		*/
		continue;					/* the line is not in the correct format			*/

	*((temp+=1)-1)='\0';					/* split the line into the varible name and the variable value	*/
	strcpy(config_var,strtoupper(trimboth(file_buffer)));	/* convert the raw variable name into a cooked variable		*/
	strcpy(config_val,trimboth(temp));			/* convert the variable value 					*/

	strcpy(config_ent,config_var);				/* copy the variable name to the entry holding buffer		*/
	strcat(config_ent,"=");					/* concatenate a delimiter to the holding buffer		*/
	strcat(config_ent,config_val);				/* concatenate the variable value to the entry holding buffer	*/

	config_table=setconfig(config_ent,config_var,config_table);	/* add the varible to the configuration table		*/
	}

fclose(conf_file);						/* close the configuration file					*/
return(config_table);						/* return the pointer to the configuration table		*/
}

/********************************************************************************************************************************
 *							setconfig								*
 ********************************************************************************************************************************/

char **setconfig(entry,variable,table)
char *entry;							/* entry to place in the configuration table			*/
char *variable;							/* variable name, determined previously no need to do it again	*/
char **table;							/* table to place variable in					*/
{
int entries;							/* the number of entries within the current config table	*/
int count;
int replace;							/* flag to indicate that the varible has been replaced		*/

char **new_table;						/* pointer to the new table when it is created			*/
char **new_temp;						/* safe house for the new table					*/
char **table_temp;						/* temp pointer for the current table				*/

char *getconfigname();						/* routine to process the entry title out of the table entry	*/


table_temp=table;						/* set the temp pointer for the table				*/
entries=0;							/* set the entries counter to zero				*/

while(*table_temp != '\0') {					/* loop through the table for entries	    			*/
	if (strcmp(variable,getconfigname(*table_temp)) != 0)	/* does the varible exist					*/
		entries++;					/* no, add the varible				      		*/
	table_temp++;						/* get the next entry						*/
	}

table_temp=table;						/* reset the temp pointer to the table				*/
new_table=(char **)malloc(sizeof(char *)*(entries+2));		/* set up the size of the new table				*/
new_temp = new_table;						/* store the table pointer in a safe place			*/

replace=0;

for(count=0; count < entries; count++) {			/* loop through the table to copy the contents to the new table */
	if (strcmp(variable,getconfigname(*table_temp)) == 0) {	/* does the variable exist					*/
		*new_table = (char *)malloc(strlen(entry)+1);	/* yes, allocate the space for the new variable	value		*/
		strcpy(*new_table,entry);			/* copy the contents to the new varible				*/
		free(*table_temp);				/* free up the space used by the old varible			*/
		replace=1;					/* set the replaced flag					*/
		}
	else {
		(*new_table)=(*table_temp);			/* copy the old contents to the new table			*/
		}

	table_temp++;						/* increment the old table					*/
        new_table++;						/* increment the new table					*/
	}

if (replace == 0) {						/* the varible was not added to the table add it now		*/
	*new_table = (char *)malloc(strlen(entry)+1);		/* yes, allocate the space for the new variable	value		*/
	strcpy(*new_table,entry);				/* copy the contents to the new table				*/
	new_table++;						/* increment the new table					*/
        }

*new_table = '\0';						/* terminate the table						*/
free(table);							/* free up the space occupied by the old table			*/
table = new_temp;						/* reset the table pointer to the new table			*/

return(new_temp);						/* return the table pointer					*/
}

/********************************************************************************************************************************
 *							getconfigname								*
 ********************************************************************************************************************************/

char *getconfigname(string)					/* return the name of a variable				*/
char *string;							/* Entry from the configuration table				*/
{
char buffer[128];

strcpy(buffer,string);
*strchr(buffer,'=') = '\0';
return(buffer);
}

/********************************************************************************************************************************
 *							initconfig								*
 ********************************************************************************************************************************/

char **initconfig()
{
char **config_table;

config_table = (char **)malloc(sizeof(char *));			/* set up the room for the null pointer				*/
*config_table = '\0';						/* init the configuration table 				*/

return(config_table);
}

/********************************************************************************************************************************
 *							getvalue								*
 ********************************************************************************************************************************/

char *getvalue(value,table)
char *value;
char **table;
{
static char buffer[256];
char *trimboth();
static char *temp;
char *getvalue();

static char *holder;

char **loop;

loop=table;							/* set the temporary pointer to that of the table		*/
while(*loop != '\0') {						/* loop while the table pointer is not null			*/
	strcpy(buffer,*loop);					/* copy the table varible to the temp buffer			*/
	*((temp=strchr(buffer,'=')+1)-1) = '\0';		/* seperate the title from the value				*/

	if (strcmp(strtoupper(buffer),trimboth(strtoupper(value))) == 0) {
		trimboth(temp);
		if (*temp == '@') {
			holder=malloc(strlen(temp));
			strcpy(holder,temp+1);			
      			temp=getvalue(holder,table);
			free(holder);
			}
		return(temp);
		}
	loop++;							/* increment the loop value					*/
	}
temp ='\0';
return(temp);
}

/********************************************************************************************************************************
 *							closeconfig								*
 ********************************************************************************************************************************/

closeconfig(table)						/* free up the space occupied by the configuration table	*/
char **table;							/* table be deleted and whose space is to be released		*/
{
char **loop;							/* temporary pointer for the config table			*/
loop = table;							/* set the pointer to the table					*/

while(*loop != '\0') {						/* loop through the table 					*/
	free(*loop);						/* free up the space used by this table parameter		*/
	loop++;							/* get the next loop parameter					*/
	}

free(table);							/* free up the space used by the table 				*/
*table='\0';
}
