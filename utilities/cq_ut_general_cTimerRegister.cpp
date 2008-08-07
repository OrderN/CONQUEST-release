// ------------------------------------------------------------------------------
// $Id$
// ------------------------------------------------------------------------------
// File cq_ut_general_cTimerRegister.cpp
// ------------------------------------------------------------------------------

// ****h* Conquest/utilities/cq_ut_general_cTimerRegister.cpp *
//
// NAME
//  cq_ut_general_cTimerRegister.cpp
// PURPOSE
//  Implementation of the cTimerRegister, a class that stores a timing line from Conquest
// USES
//  cq_ut_general_cTimerRegister.hpp
// AUTHOR
//  Antonio S. Torralba
// CREATION DATE
//  2008/07/28
// MODIFICATION HISTORY
//
// ****

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <iostream>
#include "cq_ut_general_cTimerRegister.hpp"

// -----------------------------------------------------------------------------
// Method cTimerRegister : char string
// -----------------------------------------------------------------------------

// ****m* cq_ut_general_cTimerRegister.cpp/cTimerRegister *
//
// NAME
//  cTimerRegister
// PURPOSE
//  Construction from a "Timing" line
// INPUTS
//
// AUTHOR
//  Antonio S. Torralba
// CREATION DATE
//  2008/07/28
// MODIFICATION HISTORY
//
// ****

cTimerRegister::cTimerRegister(const char *line)
 {
    int i;
    int tokcnt;
    char **tokenlist;

    tag[0]='\0';
    // Allocate space for the tokens
    tokenlist=(char **)malloc((size_t)MAX_TOK*sizeof(char *));
    for(i=0; i < MAX_TOK; ++i)
     {
        tokenlist[i]=(char *)malloc((size_t)MAX_LIN*sizeof(char));
     }
    // Make sure that the line is not too long
    if(strlen(line) > MAX_LIN)
     {
       char tmpline[MAX_LIN+1];

        // Use truncated local copy
        for(i=0; i < MAX_LIN; ++i)  tmpline[i]=line[i];
        tmpline[MAX_LIN]='\0';
        tokcnt=tokenize(tmpline,tokenlist);
     }
    else
     {
        tokcnt=tokenize(line,tokenlist);
     }

    // Assign values
    if(strcmp("Timing:",tokenlist[0]))   // This doesn't look like a Conquest timer register
     {
        // Use bad data
        processor=0;
        level=0;
        strcpy(tag,"BAD REGISTER");
        time=0.0;
     }
    else
     {
        // FOR THE MOMENT, CORRECT REGISTER FORMAT (WITH LEVEL) IS ASSUMED

        // The second token is 'Level' and the third, the value
        level=atoi(tokenlist[2]);

        // The fourth is a dash and the fifth and sixth are for the processor
        // But we need to remove the ':' from the number - THIS MAY CHANGE
        tokenlist[5][strlen(tokenlist[5])]='\0';
        processor=atoi(tokenlist[5]);

        // The next three can be ignored, and the rest are concatenated
        //   until an equal sign is found
        i=9;
        while(strcmp("=",tokenlist[i]))
         {
            strcat(tag,tokenlist[i]);
            strcat(tag," ");
            ++i;
            if(i==tokcnt)    // Time missing
             {
                strcpy(tag,"BAD REGISTER");
                time=0.0;
                break;
             }
         }
        if(i != tokcnt)   time=atof(tokenlist[i+1]);
     }

    // Free the token list
    for(i=0; i < MAX_TOK; ++i)
     {
        free(tokenlist[i]);
     }
    free(tokenlist);
 }

// -----------------------------------------------------------------------------
// Method cTimerRegister : data
// -----------------------------------------------------------------------------

// ****m* cq_ut_general_cTimerRegister.cpp/cTimerRegister *
//
// NAME
//  cTimerRegister
// PURPOSE
//  Construction from actual pieces of data: processor, level, tag, time
// INPUTS
//
// AUTHOR
//  Antonio S. Torralba
// CREATION DATE
//  2008/07/28
// MODIFICATION HISTORY
//
// ****

cTimerRegister::cTimerRegister(int p, int l, const char *t, double tm)
 {
    int i;
    int min,tmpl;

    processor=p;
    level=l;
    if((tmpl=strlen(t)) < MAX_LIN)   min=tmpl;
    else                           min=MAX_LIN;
    for(i=0; i < min; ++i)         tag[i]=t[i];
    tag[min]='\0';
    time=tm;
 }

// -----------------------------------------------------------------------------
// Method operator== : int
// -----------------------------------------------------------------------------

// ****m* cq_ut_general_cTimerRegister.cpp/operator== *
//
// NAME
//  operator==
// PURPOSE
//  Used implicitly by find algorithms that search by level
// INPUTS
//
// AUTHOR
//  Antonio S. Torralba
// CREATION DATE
//  2008/07/28
// MODIFICATION HISTORY
//
// ****

bool cTimerRegister::operator==(int v)
 {
     return (level==v);
 }

// -----------------------------------------------------------------------------
// Method operator== : long
// -----------------------------------------------------------------------------

// ****m* cq_ut_general_cTimerRegister.cpp/operator== *
//
// NAME
//  operator==
// PURPOSE
//  Used implicitly by find algorithms that search by processor
// INPUTS
//
// AUTHOR
//  Antonio S. Torralba
// CREATION DATE
//  2008/07/28
// MODIFICATION HISTORY
//
// ****

bool cTimerRegister::operator==(long v)
 {
     return (processor==v);
 }

// -----------------------------------------------------------------------------
// Method operator!= : int
// -----------------------------------------------------------------------------

// ****m* cq_ut_general_cTimerRegister.cpp/operator!= *
//
// NAME
//  operator!=
// PURPOSE
//  Used implicitly by find algorithms that search by level
// INPUTS
//
// AUTHOR
//  Antonio S. Torralba
// CREATION DATE
//  2008/07/28
// MODIFICATION HISTORY
//
// ****

bool cTimerRegister::operator!=(int v)
 {
     return (level!=v);
 }

// -----------------------------------------------------------------------------
// Method operator!= : long
// -----------------------------------------------------------------------------

// ****m* cq_ut_general_cTimerRegister.cpp/operator!= *
//
// NAME
//  operator!=
// PURPOSE
//  Used implicitly by find algorithms that search by processor
// INPUTS
//
// AUTHOR
//  Antonio S. Torralba
// CREATION DATE
//  2008/07/28
// MODIFICATION HISTORY
//
// ****

bool cTimerRegister::operator!=(long v)
 {
     return (processor!=v);
 }

// -----------------------------------------------------------------------------
// Method print
// -----------------------------------------------------------------------------

// ****m* cq_ut_general_cTimerRegister.cpp/print *
//
// NAME
//  print
// PURPOSE
//  Print a line in a semi-informative format (use only for debugging)
// INPUTS
//
// AUTHOR
//  Antonio S. Torralba
// CREATION DATE
//  2008/07/28
// MODIFICATION HISTORY
//
// ****

void cTimerRegister::print()
 {
    printf("P: %5d  L: %3d  Time: %12.6lf s   T: %s\n",processor,level,time,tag);
 }

// -----------------------------------------------------------------------------
// Method tokenize
// -----------------------------------------------------------------------------

// ****m* cq_ut_general_cTimerRegister.cpp/tokenize *
//
// NAME
//  tokenize
// PURPOSE
//  Break a Conquest "Timing" line in its words, with space as the separator
// INPUTS
//
// AUTHOR
//  Antonio S. Torralba
// CREATION DATE
//  2008/07/28
// MODIFICATION HISTORY
//
// ****

int cTimerRegister::tokenize(const char *lin, char **tok)
 {
   size_t length;
   int i,tokcnt,charcnt,ant;

   length=strlen(lin);
   tokcnt=0;
   charcnt=0;
   ant=0;
   for(i=0; i < length; ++i)
    {
       if(lin[i] != ' ')
        {
           tok[tokcnt][charcnt]=lin[i];
           ++charcnt;
           ant=1;
        }
       else
        { 
           if(ant==1)
            {
               tok[tokcnt][charcnt]='\0';
               ++tokcnt;
               ant=0;
            }
           charcnt=0;
        }
    }
   if(lin[length-1] != ' ')   // Close last token, if open
    {
      tok[tokcnt][charcnt]='\0';
      ++tokcnt;
    }

   return tokcnt;
 }

