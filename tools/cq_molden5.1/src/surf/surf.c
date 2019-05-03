/* Copyright amitabh varshney@cs.unc.edu. All Rights Reserved. March 16, 1994 */
/*--------------------------------------------------------------------------------
surf.c

This is the main file that defines all the globals and starts everything.
----------------------------------------------------------------------------------*/
#define	EXTERN	
#include "surf.h"

void
main(ac,av)
int	ac;
char*	av[];
{

   char                 filename[128], outfile[128];

   Checks_On = FALSE;
   Write_Option = 0;
   Num_atoms = 0;
   Probe_radius = -1;
   Max_Tess_Len = -1;
   Max_Gp_Polys = 0;

   deal_wif_options(ac,av, filename);  /* set the filename */

   if ((Write_Option == 2)&&(Max_Gp_Polys == 0))
   { printf("Maximum expected size of the output surface (in thousands of tris): ");
     scanf("%d", &Max_Gp_Polys);
     Max_Gp_Polys *= 1000;
   }
   sprintf(outfile, "%s.tri", filename);

   /* read in the dataset */
   input_dataset(filename);

   if (Write_Option)
     begin_output_dataset(outfile);

   printf("Constructing solvent-accessible surface ..\n");

   /* start timing the molecular surface computations here */
   START

   /* initialize and compute the molecular surface */
   init_and_compute(); 

   /* stop timing the molecular surface computations here */
   STOP

   if (Write_Option == 1)
     printf("Surface construction + writing time %4.2f seconds\n", et);
   else
     printf("Surface construction time %4.2f seconds\n", et);

   if (Write_Option == 2) output_dataset();

   if (Write_Option) end_output_dataset();
}


/*---------------------------------------------------------------------------------
usage Prints out the command line options 
---------------------------------------------------------------------------------*/
usage(prog_name)
char	*prog_name;
{
   fprintf(stderr,"usage: %s [options] <data_file>\n", prog_name);
   fprintf(stderr,"-R flpt    initial probe radius (default 1.4)\n");
   fprintf(stderr,"-E flpt    maximum triangle edge length (default 1.4)\n");
   fprintf(stderr,"-W 0       don't output any tris, for timing (default)\n");
   fprintf(stderr,"-W 1       write the triangles as they are generated\n");
   fprintf(stderr,"-W 2       store up the tris till the end & then write\n");
   fprintf(stderr,"-T int     size of the tri buffer (in thousands)\n");
   fprintf(stderr,"-C         turn on the checks for robustness\n");
   fprintf(stderr,"Output triangles are written to <data_file.tri>\n");
   exit(1);
}

/*---------------------------------------------------------------------------------
deal_wif_options intializes the various command line parameters mentioned above
---------------------------------------------------------------------------------*/
deal_wif_options(ac,av, filename)
int	  	ac;
char*	  	av[];
char* 	  	filename;
{
   char 	*s;
   int	 	argumentOk=0;
   char		*prog_name;

   prog_name = av[0];
   while (--ac > 0 )
   {
      if ((*++av)[0]=='-'){
          for (s = av[0]+1; *s; s++)  switch (*s) {
           case 'R':
	       (void) sscanf(*++av,"%f",&Probe_radius);
	       ac-=1;
	       break;
           case 'E':
	       (void) sscanf(*++av,"%lf",&Max_Tess_Len);
	       ac-=1;
	       break;
           case 'C':
	       Checks_On = TRUE;
	       break;
           case 'W':
	       (void) sscanf(*++av,"%d",&Write_Option);
	       ac-=1;
	       break;
           case 'T':
	       (void) sscanf(*++av,"%d",&Max_Gp_Polys);
               Max_Gp_Polys *= 1000;
	       ac-=1;
	       break;
	   case '?':
	       usage(prog_name);
	       break;
	   default:
	       break;
         };     /*of switch */ 
       }   /* of if */ 
       else {
	   sprintf((char *) filename, "%s", *av); 
	   argumentOk+=1;
       }
    }
    if (argumentOk!=1) usage(prog_name);
}
