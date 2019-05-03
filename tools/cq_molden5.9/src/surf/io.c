/* Copyright amitabh varshney@cs.unc.edu. All Rights Reserved. March 16, 1994 */
/*--------------------------------------------------------------------------------
io.c

This file contains the procedures for input and output.
Any interface to the .pdb files, as well as to the output graphics libraries
would go here.
----------------------------------------------------------------------------------*/
#include "surf.h"

/*--------------------------------------------------------------------------------
input_dataset() reads in the atoms and their radii from the input file. At present
the format of this input file is as follows:

<atom_id> <atom_radius> <atom_center_x> <atom_center_y> <atom_center_z>

The program right now ignores all the <atom_id> values in the input file 
but in the output it tags each triangle with the sequential location of its 
atom position in the input file.
----------------------------------------------------------------------------------*/
input_dataset(filename)
char		*filename;
{ FILE  	*fp;
  char  	line[128];
  int		i = 0;
  int		dummy;

  { if (!(fp = fopen(filename, "r")))
    { fprintf(stderr,"Unable to open input file %s\n", filename);
      exit(-1) ;
    }
    /* get the number of atoms */
    while (fgets(line, sizeof(line), fp)) i++;
    fseek(fp, 0L, 0);

    Num_atoms = i;

    ALLOCN(atoms, Gp_Atom, Num_atoms);

    printf("Reading %d atoms..", Num_atoms);

    /* input all atoms */
    for(i = 0; i < Num_atoms; i++)
    { fgets(line, sizeof(line), fp);
      sscanf(line, "%d %f %f %f %f", &dummy, &atoms[i].radius, 
      &atoms[i].center[X], &atoms[i].center[Y], &atoms[i].center[Z]);
      atoms[i].type = 0;
    }

    printf("done\n");

    if (Probe_radius == -1) Probe_radius = 1.4;
    printf("Probe radius = %2.3f\n", Probe_radius);
     
    fclose(fp);
  }
}

/*--------------------------------------------------------------------------------
 *output_dataset() write out the triangles associated with the molecular surface. 

the program outputs triangles in the following format:
<atom_id>
<coord0_x> <coord0_y> <coord0_z> <normal0_x> <normal0_y> <normal0_z>
<coord1_x> <coord1_y> <coord1_z> <normal1_x> <normal1_y> <normal1_z>
<coord2_x> <coord2_y> <coord2_z> <normal2_x> <normal2_y> <normal2_z>
----------------------------------------------------------------------------------*/

/*---------------------------------------------------------------------------------
begin_output_dataset initializes/opens the output file for writing
---------------------------------------------------------------------------------*/
begin_output_dataset(outfilename)
char		outfilename[128];
{
  if (!(Opf = fopen(outfilename, "w")))
  { fprintf(stderr,"Unable to open output file %s\n", outfilename);
    exit(-1) ;
  }
/* ---------------------archive format-------------- 
  fprintf(Opf,"\nstructure Surface posted {");
  fprintf(Opf,"\ncolor polygon {");
  fprintf(Opf,"\ndiff 100 100 250 0 100 100 250 0");
  fprintf(Opf,"\nspec 255 255 255 1 255 255 255 1\n};");
 -----------------------------------------------------*/
}

/*---------------------------------------------------------------------------------
output_dataset writes out all the triangles in one go
---------------------------------------------------------------------------------*/
output_dataset()
{ 
  int   i;

  printf("Writing out %d triangles", Num_polys);
  for(i = 0; i < Num_polys; i++)
  { if (i%1000==0)
    { printf(".");
      fflush(stdout);
    }
    write_archive_tri(&verts[i][0], &verts[i][1], &verts[i][2], (int)atom_type[i]);
  }
}


/*---------------------------------------------------------------------------------
write_archive_tri writes out the triangles 
---------------------------------------------------------------------------------*/
write_archive_tri(pt0, pt1, pt2, atype)
VertexType	*pt0, *pt1, *pt2;
int		atype;
{ 

#if 0 /*--------------archive format---------*/
  NORMALIZE3(pt0->Normal); NORMALIZE3(pt1->Normal); NORMALIZE3(pt2->Normal);
  fprintf(Opf,"\npolygon 3 {");
  fprintf(Opf,"\n%3.4lf %3.4lf %3.4lf %3.4lf %3.4lf %3.4lf;",
                pt0->Coord[X],  pt0->Coord[Y],  pt0->Coord[Z], 
                pt0->Normal[X], pt0->Normal[Y], pt0->Normal[Z]);
  fprintf(Opf,"\n%3.4lf %3.4lf %3.4lf %3.4lf %3.4lf %3.4lf;",
                pt1->Coord[X],  pt1->Coord[Y],  pt1->Coord[Z], 
                pt1->Normal[X], pt1->Normal[Y], pt1->Normal[Z]);
  fprintf(Opf,"\n%3.4lf %3.4lf %3.4lf %3.4lf %3.4lf %3.4lf;",
                pt2->Coord[X],  pt2->Coord[Y],  pt2->Coord[Z], 
                pt2->Normal[X], pt2->Normal[Y], pt2->Normal[Z]);
  fprintf(Opf,"\n};");

#else
#if 1 /*---------------default format-----------*/
  fprintf(Opf,"%d\n",atype);
  fprintf(Opf,"%3.4lf %3.4lf %3.4lf %3.4lf %3.4lf %3.4lf\n",
                pt0->Coord[X],  pt0->Coord[Y],  pt0->Coord[Z], 
                pt0->Normal[X], pt0->Normal[Y], pt0->Normal[Z]);
  fprintf(Opf,"%3.4lf %3.4lf %3.4lf %3.4lf %3.4lf %3.4lf\n",
                pt1->Coord[X],  pt1->Coord[Y],  pt1->Coord[Z], 
                pt1->Normal[X], pt1->Normal[Y], pt1->Normal[Z]);
  fprintf(Opf,"%3.4lf %3.4lf %3.4lf %3.4lf %3.4lf %3.4lf\n",
                pt2->Coord[X],  pt2->Coord[Y],  pt2->Coord[Z], 
                pt2->Normal[X], pt2->Normal[Y], pt2->Normal[Z]);

#else /*---------------poly_format-------------- */
  fprintf(Opf,"3\n");
  fprintf(Opf,"%3.4lf %3.4lf %3.4lf \n",
                pt0->Coord[X],  pt0->Coord[Y],  pt0->Coord[Z]);
  fprintf(Opf,"%3.4lf %3.4lf %3.4lf \n",
                pt1->Coord[X],  pt1->Coord[Y],  pt1->Coord[Z]);
  fprintf(Opf,"%3.4lf %3.4lf %3.4lf \n",
                pt2->Coord[X],  pt2->Coord[Y],  pt2->Coord[Z]);
#endif
#endif
}

/*---------------------------------------------------------------------------------
end_output_dataset closes the output file
---------------------------------------------------------------------------------*/
end_output_dataset()
{ 
/*
  fprintf(Opf,"\n};"); 
*/
  fflush(Opf);
  fclose(Opf);
  printf("done\n");
}

