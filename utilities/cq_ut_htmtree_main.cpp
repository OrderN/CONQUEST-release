// ------------------------------------------------------------------------------
// $Id$
// ------------------------------------------------------------------------------
// File cq_ut_htmtree_main.cpp
// ------------------------------------------------------------------------------

// ****h* Conquest/utilities/cq_ut_htmtree_main.cpp *
//
// NAME
//  cq_ut_htmtree
// PURPOSE
//  Creates an html tree from Conquest_out or Conquest time files
//   It reads the timing lines and reorders them in execution order
//   to facilitate their interpretation
// INPUTS
//  Either standard input (with no command line parameters or if the first one is -) 
//   or an input file given in the command line. The output file is execpage.html,
//   unless another one is specified as the second command line parameter
//   NOTE: The images "needed" by the html are under res/cq_ut_htmtree/
// USES
//  cq_ut_general_cTimerRegister
//  cq_ut_general_cTimerList
// AUTHOR
//  Antonio S. Torralba
// CREATION DATE
//  2008/07/28
// MODIFICATION HISTORY
//
// ****

#include <stdio.h>
#include <stdlib.h>
#include "cq_ut_general_cTimerList.hpp"

// -----------------------------------------------------------------------------
// Function main
// -----------------------------------------------------------------------------

// ****f* cq_ut_htmtree_main.cpp/main *
//
// NAME
//  main
// USAGE
//
// PURPOSE
//  Read Conquest output files with timing lines and create an html file
//   easier to interprete
// INPUTS
//
// AUTHOR
//  Antonio S. Torralba
// CREATION DATE
//  2008/07/28
// MODIFICATION HISTORY
//  
// ****

int main(int argc, char *argv[])
 {
    char inputfile[MAX_LIN];

    if(argc < 2)
     {
         // Print the information anyway
         fprintf(stderr,"Use: %s <Conquest_output_file> [<execpage.html>]\n", argv[0]);
         fprintf(stderr,"     The <Conquest_output_file> can be \"-\" to read from stdin\n");
         fprintf(stderr,"     Called without parameters. Reading from stdin\n\n");
//         exit(-1);
         strcpy(inputfile,"-");     // Without parameters, read from standard input
     }
    else
     {
         strcpy(inputfile,argv[1]);
     }

    // Create a list from the records in the file and print to html
    cTimerList TmrList(inputfile);
    if(argc >= 3)    TmrList.printHtml(argv[2]); 
    else             TmrList.printHtml("execpage.html"); 
 }
 
