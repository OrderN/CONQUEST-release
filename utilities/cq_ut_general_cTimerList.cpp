// ------------------------------------------------------------------------------
// $Id$
// ------------------------------------------------------------------------------
// File cq_ut_general_cTimerList.cpp
// ------------------------------------------------------------------------------

// ****h* Conquest/utilities/cq_ut_general_cTimerList.cpp *
//
// NAME
//  cq_ut_general_cTimerList.cpp
// PURPOSE
//  Implementation of cTimerList, a class that stores Conquest timing lines 
//    (cTimerRegister objects) in a list and sorts it in execution order. 
//    It can print an html with the execution tree
// USES
//  Standard Template Library: list, algorithms
//  cq_ut_general_cTimerRegister.hpp
//  cq_ut_general_cTimerList.hpp
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
#include <algorithm>
#include "cq_ut_general_cTimerRegister.hpp"
#include "cq_ut_general_cTimerList.hpp"

// -----------------------------------------------------------------------------
// Method cTimerList : char string
// -----------------------------------------------------------------------------

// ****m* cq_ut_general_cTimerList.cpp/cTimerList *
//
// NAME
//  cTimerList
// PURPOSE
//  Construction from a file. The list is ordered by processor and level at creation time
// INPUTS
//
// AUTHOR
//  Antonio S. Torralba
// CREATION DATE
//  2008/07/28
// MODIFICATION HISTORY
//
// ****

cTimerList :: cTimerList(const char *file)
 {
    FILE *fp0;
    char linea[MAX_LIN];
    char buffer[MAX_LIN];
    int i, no_nodes, done_nodes, next;
    list<cTimerRegister>::iterator it_ini;
    list<cTimerRegister>::iterator it_end;

    if(strcmp(file,"-"))
     {
      if((fp0=fopen(file,"r")) == NULL)
       {
         fprintf(stderr,"File '%s' couldn't be opened\n",file);
         exit(-1);
       }
     }
    else fp0=stdin;

    // Create a list from the records in the file
    printf("Reading Timing lines with level\n");
    while(fgets(linea,MAX_LIN+1,fp0) != NULL)
     {
        if(!strncmp(linea,"Timing: Level",13))
         {
           cTimerRegister reg(linea);   // Create register

           // Insert new element in the list of registers
           l.push_back(reg);
         }
     }
    fclose(fp0);

    no_nodes=orderByProcessor();               // Order the list in blocks by processor
    printf("NODES = %4d\n",no_nodes);
    done_nodes=0;
    i=0;
    while(done_nodes < no_nodes)
     {
       ++i;
       it_ini=find(l.begin(),l.end(),(long)(i));    // Conversion necessary to check by processor, not level       
       // Find the next processor number (there might be gaps)
       if((done_nodes+1)==no_nodes)   it_end=l.end();     // For the last processor, end of interval = end of list
       else
        {
          next=i+1;                                       // In principle, the next value
          it_end=find(l.begin(),l.end(),(long)(next));    // Conversion necessary to check by processor, not level       
          while(it_end == l.end())
           {
             ++next;
             it_end=find(l.begin(),l.end(),(long)(next)); // Conversion necessary to check by processor, not level       
           }
        }
       if(it_ini != l.end())                              // If the processor is present
        {
          ++done_nodes;
          printf("DOING NODE %d\n",i);
          printf("Fixing possible truncations\n");
          fixList(it_ini,it_end);                         // Fix possible gaps in the list
          printf("Sorting info by execution order\n");
          orderExecution(it_ini,it_end,1);                // Reorder the list in execution order
          it_ini=find(l.begin(),l.end(),(long)(i));       // We need to find the start of the node again
                                                          // Some registers could have been inserted before the first
          printf("Calculating subtotals and untimed differences\n");
          calculateSubtotals(it_ini,it_end);
          if(no_nodes > 1)                                // For more than one processor
           {
             sprintf(buffer,"PROCESSOR %4d",i);
             cTimerRegister tmp_reg(i,0,buffer,0.0);
             it_ini=find(l.begin(),l.end(),(long)(i));    // Find the first register of the node again       
             l.insert(it_ini,tmp_reg);                    // This will make printing easier
           }
        }
     }

// Show that the list is correctly formed
//list<cTimerRegister>::iterator iter;
//for(iter=l.begin(); iter != l.end(); ++iter)   iter->print();
//printf("%d and %d\n", l.size(), l.max_size());
 }
 
// -----------------------------------------------------------------------------
// Method comparePerProcessor
// -----------------------------------------------------------------------------

// ****m* cq_ut_general_cTimerList.cpp/comparePerProcessor *
//
// NAME
//  comparePerProcessor
// PURPOSE
//  Ordering method used implicitly by sorting algorithms
// INPUTS
//
// AUTHOR
//  Antonio S. Torralba
// CREATION DATE
//  2008/07/28
// MODIFICATION HISTORY
//
// ****

bool cTimerList :: comparePerProcessor(cTimerRegister r1, cTimerRegister r2)
 { 
    if(r1.processor < r2.processor)   return true;
    else                              return false;
 }

// -----------------------------------------------------------------------------
// Method orderByProcessor
// -----------------------------------------------------------------------------

// ****m* cq_ut_general_cTimerList.cpp/orderByProcessor *
//
// NAME
//  orderByProcessor
// PURPOSE
//  Orders the list in increasing processor order
// INPUTS
//
// AUTHOR
//  Antonio S. Torralba
// CREATION DATE
//  2008/07/28
// MODIFICATION HISTORY
//
// ****

int cTimerList :: orderByProcessor()
 {
    int no_process, prev;
    list<cTimerRegister>::iterator iter;

    l.sort(cTimerList::comparePerProcessor);
    no_process=0;
    prev=-1;
    iter=l.begin();
    while(iter != l.end())
     {
        if(iter->processor != prev)  ++no_process;
        prev=iter->processor;
        ++iter;
     }
    return no_process; 
 }

// -----------------------------------------------------------------------------
// Method orderExecution
// -----------------------------------------------------------------------------

// ****m* cq_ut_general_cTimerList.cpp/orderExecution *
//
// NAME
//  orderExecution
// PURPOSE
//  Sorts a fragment of the list (for one processor & level) in the order 
//    of execution of the timers, assuming that the lines were written
//    when each timer was stopped 
// INPUTS
//
// AUTHOR
//  Antonio S. Torralba
// CREATION DATE
//  2008/07/28
// MODIFICATION HISTORY
//
// ****

void cTimerList :: orderExecution(list<cTimerRegister>::iterator i1,list<cTimerRegister>::iterator i2,int level)
 {
    list<cTimerRegister>::iterator end_it;
 
    while(i1 != i2)
     {
       // Find the subinterval to order in the next level  
       end_it=find(i1,i2,level);
       if(end_it==i1)          // This block is just this register; move to next
        {
           ++i1;
        }
       else
        {      
            l.insert(i1,*end_it);                 // Move last element to the beggining
            orderExecution(i1,end_it,level+1);    // Reorder the rest of this block
            ++end_it;                             // Make a note of the first register after block
            i1=end_it;
            --end_it;
            l.erase(end_it);                      // Finally, remove the old register at the end
        }
     }
 }

// -----------------------------------------------------------------------------
// Method assignSubtotals
// -----------------------------------------------------------------------------

// ****m* cq_ut_general_cTimerList.cpp/assignSubtotals *
//
// NAME
//  assignSubtotals
// PURPOSE
//  Calculates the untimed remainders for complete levels and the total times spent
//   in completed tasks in truncated levels. Inserts lines in the list at the right places
// INPUTS
//
// AUTHOR
//  Antonio S. Torralba
// CREATION DATE
//  2008/07/28
// MODIFICATION HISTORY
//
// ****

// This expects execution order already
double cTimerList :: assignSubtotals(list<cTimerRegister>::iterator i1, list<cTimerRegister>::iterator i2, int level)
 {
    list<cTimerRegister>::iterator next_it;
    list<cTimerRegister>::iterator tmp_it;
    double tmptime,total, block_next,total_next;
 
    total=0.0;
/*
    if(i1->level != level)     // First level of the interval is not the requested
     {
       i1=find(i1,i2,level);
     }
*/
    while(i1 != i2)
     {
       tmptime=i1->time;
       tmp_it=i1;                 // Just in case we need it for a truncated record
       ++i1;                      // Find the next register of this level
       next_it=find(i1,i2,level);
       if(next_it!=i1)            // The previous block was not just one register
        {
          // Get the subtotal for the next level in this block
          block_next = assignSubtotals(i1,next_it,level+1);
          total_next += block_next;
          if(tmptime < 0.0)                                 // This means truncated. We want a total, not a difference
           {
             cTimerRegister new_reg(i1->processor,level+1,"TOTAL",block_next);    
             l.insert(next_it,new_reg);
             tmp_it->time=block_next;                       // Store the total time of the contents of a truncated record
           }
          else                                              // We want the non-timed difference
           {
             cTimerRegister new_reg(i1->processor,level+1,"ALL THE REST",tmptime-block_next);    
             l.insert(next_it,new_reg);
           }
          i1=next_it;
        }
       total += tmp_it->time;
     }
   return total;
 }

// -----------------------------------------------------------------------------
// Method calculateSubtotals
// -----------------------------------------------------------------------------

// ****m* cq_ut_general_cTimerList.cpp/calculateSubtotals *
//
// NAME
//  calculateSubtotals
// PURPOSE
//  This is the outer level of recursion, which calls assignSubtotals
//  At the moment, both methods are private. Only this one should/could be made public
// INPUTS
//
// AUTHOR
//  Antonio S. Torralba
// CREATION DATE
//  2008/07/28
// MODIFICATION HISTORY
//
// ****

// This function is the outer level, so it's the one to be called
void cTimerList :: calculateSubtotals(list<cTimerRegister>::iterator i1,list<cTimerRegister>::iterator i2)
 {
    double total;

    total=assignSubtotals(i1,i2,1);

    // This total is just appended at the interval
    cTimerRegister tail_reg(i1->processor,1,"TOTAL",total);    
    // Find the last register of the block again (it has changed due to insertions)
    i2=find(l.begin(),l.end(),(long)(i1->processor));   // This finds the first register for the processor
    while(i2 != l.end())                                // Now, we take advantage of the list being ordered
     {
        ++i2;
        if(i2->processor != i1->processor)    break;    // This is the beggining of the next interval
     }
    l.insert(i2,tail_reg);
 }

// -----------------------------------------------------------------------------
// Method fixList
// -----------------------------------------------------------------------------

// ****m* cq_ut_general_cTimerList.cpp/fixList *
//
// NAME
//  fixList
// PURPOSE
//  Annotates the list to indicate that some levels are missing (truncated)
// INPUTS
//
// AUTHOR
//  Antonio S. Torralba
// CREATION DATE
//  2008/07/28
// MODIFICATION HISTORY
//
// ****

void cTimerList :: fixList(list<cTimerRegister>::iterator i1,list<cTimerRegister>::iterator i2)
 {
    list<cTimerRegister>::iterator tmp_it,buf_it;    // For temporary use
    int last_level,ant_level;
    int i;

    // There are to possible things to fix
    //  1. Truncation at the end of the list, because Conquest did not fisish
    //     In a good list, the last level should be 1
    //  2. Truncation of intermediate levels because of the user-chosen level of output

    // Fix of 1
    tmp_it=i2;
    --tmp_it;
    cTimerRegister tmp_reg(tmp_it->processor,0,"TRUNCATED",-1.0);
    last_level=tmp_it->level;
    if(last_level != 1)
     {
        for(i=last_level-1; i >= 1; --i)
         {
            tmp_reg.level=i;
            l.insert(i2,tmp_reg);
         } 
     }
    
    // Fix 2: Make sure that all the intermediate levels are there
    //   Read the list backwards. If it works, the level will either decrease (always correct)
    //   or increase one by one. Truncations are signaled by jumps up by more than one level
    tmp_it=i2;
    --tmp_it;
    ant_level=tmp_it->level;
    while(tmp_it != i1)
     {
        --tmp_it;
        last_level=tmp_it->level;
        if(last_level > ant_level && last_level != ant_level+1)
         {
            ++tmp_it;
            for(i=last_level-1; i > ant_level; --i)
             {
                tmp_reg.level=i;
                l.insert(tmp_it,tmp_reg);
             }
            --tmp_it;
         }
        ant_level=tmp_it->level;
     }
 }

// -----------------------------------------------------------------------------
// Method printHtml
// -----------------------------------------------------------------------------

// ****m* cq_ut_general_cTimerList.cpp/printHtml *
//
// NAME
//  printHtml
// PURPOSE
//  Does something useful: Prints an html that easier to understand than 
//    the output from Conquest
// INPUTS
//
// AUTHOR
//  Antonio S. Torralba
// CREATION DATE
//  2008/07/28
// MODIFICATION HISTORY
//
// ****

// This should be done (more mantainable) using classes and configuration files
void cTimerList :: printHtml(const char *file)
 {
    list<cTimerRegister>::iterator it;
    FILE *fp0;
    int level1,level2;
    int labelcount;
    int openspans;
    char label[MAX_LIN];
    char class1[MAX_LIN],class2[MAX_LIN],image[MAX_LIN],branch[MAX_LIN];
    int i;
    int no_process;

    printf("Creating html file: %s\n", file);
    if((fp0=fopen(file,"w"))==NULL)
     {
         fprintf(stderr,"Couldn't create '%s'\n",file);
         exit(-1);
     }

// Headers (Title, CSS, Javascript, etc)
//   This could be (and was originally) automatically generated from an auxiliary html file
fprintf(fp0,"<html>\n");
fprintf(fp0,"<head>\n");
fprintf(fp0,"<META http-equiv=\"Content-Type\" content=\"text/html; charset=UTF-8\">\n");
fprintf(fp0,"<title>Conquest Execution Tree</title>\n");
fprintf(fp0,"<link href=\"exectree.css\" type=\"text/css\" rel=\"stylesheet\">\n");
fprintf(fp0,"<style type=\"text/css\">\n");
fprintf(fp0,"body{\n");
fprintf(fp0,"        font: 10pt Verdana,sans-serif;\n");
fprintf(fp0,"        color: navy;\n");
fprintf(fp0,"}\n");
fprintf(fp0,".title{\n");
fprintf(fp0,"        font: 25pt Verdana,sans-serif;\n");
fprintf(fp0,"        font-weight: bold\n");
fprintf(fp0,"        color: #000040\n");
fprintf(fp0,"}\n");
fprintf(fp0,".trigger{\n");
fprintf(fp0,"        cursor: pointer;\n");
fprintf(fp0,"        cursor: hand;\n");
fprintf(fp0,"        display: block;\n");
fprintf(fp0,"}\n");
fprintf(fp0,".branch{\n");
fprintf(fp0,"        display: none;\n");
fprintf(fp0,"        margin-left: 16px;\n");
fprintf(fp0,"}\n");
fprintf(fp0,".branch1{\n");
fprintf(fp0,"        display: block;\n");
fprintf(fp0,"        margin-left: 16px;\n");
fprintf(fp0,"}\n");
fprintf(fp0,"table.timer{\n");
fprintf(fp0,"text-align: center;\n");
fprintf(fp0,"font: 10pt Verdana,sans-serif;\n");
fprintf(fp0,"color: #404040;\n");
fprintf(fp0,"border: 0 px;\n");
fprintf(fp0,"border-spacing: 0px;\n");
fprintf(fp0,"}\n");
fprintf(fp0,".process{\n");
fprintf(fp0,"background-color: #f0f0f0;\n");
fprintf(fp0,"text-align: left;\n");
fprintf(fp0,"font-family: Verdana;\n");
fprintf(fp0,"font-weight: bold;\n");
fprintf(fp0,"font-size: 13pt;\n");
fprintf(fp0,"color: #702020;\n");
fprintf(fp0,"width: 170px;\n");
fprintf(fp0,"}\n");
fprintf(fp0,".empty{\n");
fprintf(fp0,"background-color: #f0f0f0;\n");
fprintf(fp0,"text-align: left;\n");
fprintf(fp0,"font-family: Verdana;\n");
fprintf(fp0,"font-weight: bold;\n");
fprintf(fp0,"font-size: 13pt;\n");
fprintf(fp0,"color: #702020;\n");
fprintf(fp0,"width: 5px;\n");
fprintf(fp0,"}\n");
fprintf(fp0,".tag{\n");
fprintf(fp0,"border-bottom: 1px solid #d79900;\n");
fprintf(fp0,"text-align: left;\n");
fprintf(fp0,"font-family: Verdana;\n");
fprintf(fp0,"font-weight: bold;\n");
fprintf(fp0,"font-size: 10pt;\n");
fprintf(fp0,"color: #404040;\n");
fprintf(fp0,"width: 300px;\n");
fprintf(fp0,"}\n");
fprintf(fp0,".total{\n");
fprintf(fp0,"border-bottom: 1px solid #d79900;\n");
fprintf(fp0,"background-color: #ffff55;\n");
fprintf(fp0,"text-align: left;\n");
fprintf(fp0,"font-family: Verdana;\n");
fprintf(fp0,"font-weight: bold;\n");
fprintf(fp0,"font-size: 10pt;\n");
fprintf(fp0,"color: #404040;\n");
fprintf(fp0,"width: 300px;\n");
fprintf(fp0,"}\n");
fprintf(fp0,".rest{\n");
fprintf(fp0,"border-bottom: 1px solid #d79900;\n");
fprintf(fp0,"background-color: #55ff55;\n");
fprintf(fp0,"text-align: left;\n");
fprintf(fp0,"font-family: Verdana;\n");
fprintf(fp0,"font-weight: bold;\n");
fprintf(fp0,"font-size: 10pt;\n");
fprintf(fp0,"color: #404040;\n");
fprintf(fp0,"width: 300px;\n");
fprintf(fp0,"}\n");
fprintf(fp0,".trunc{\n");
fprintf(fp0,"border-bottom: 1px solid #d79900;\n");
fprintf(fp0,"background-color: #e8f891;\n");
fprintf(fp0,"text-align: left;\n");
fprintf(fp0,"font-family: Verdana;\n");
fprintf(fp0,"font-weight: bold;\n");
fprintf(fp0,"font-size: 10pt;\n");
fprintf(fp0,"color: #404040;\n");
fprintf(fp0,"width: 300px;\n");
fprintf(fp0,"}\n");
fprintf(fp0,".time{\n");
fprintf(fp0,"border-bottom: 1px solid #d799ff;\n");
fprintf(fp0,"background-color: #f0f0f0;\n");
fprintf(fp0,"text-align: right;\n");
fprintf(fp0,"font-family: Verdana;\n");
fprintf(fp0,"font-weight: bold;\n");
fprintf(fp0,"font-size: 10pt;\n");
fprintf(fp0,"color: #404040;\n");
fprintf(fp0,"width: 100px;\n");
fprintf(fp0,"}\n");
fprintf(fp0,".calculated{\n");
fprintf(fp0,"border-bottom: 1px solid #d799ff;\n");
fprintf(fp0,"background-color: #a0e0ff;\n");
fprintf(fp0,"text-align: right;\n");
fprintf(fp0,"font-family: Verdana;\n");
fprintf(fp0,"font-weight: bold;\n");
fprintf(fp0,"font-size: 10pt;\n");
fprintf(fp0,"color: #404040;\n");
fprintf(fp0,"width: 100px;\n");
fprintf(fp0,"}\n");
fprintf(fp0,"a{\n");
fprintf(fp0,"        text-decoration: none;\n");
fprintf(fp0,"}\n");
fprintf(fp0,"a:hover{\n");
fprintf(fp0,"        text-decoration: underline;\n");
fprintf(fp0,"}\n");
fprintf(fp0,"</style>\n");
fprintf(fp0,"<script type=\"text/javascript\">\n");
fprintf(fp0,"var openImg = new Image();\n");
fprintf(fp0,"openImg.src = \"minus.gif\";\n");
fprintf(fp0,"var closedImg = new Image();\n");
fprintf(fp0,"closedImg.src = \"plus.gif\";\n");
fprintf(fp0,"\n");
fprintf(fp0,"function openFirst() {\n");
fprintf(fp0,"       var objFirst = document.getElementById('label0000000001').style;\n");
fprintf(fp0,"       objFirst.display=\"block\";\n");
fprintf(fp0,"}\n");
fprintf(fp0,"\n");
fprintf(fp0,"function showBranch(branch){\n");
fprintf(fp0,"        var objBranch = document.getElementById(branch).style;\n");
fprintf(fp0,"        if(objBranch.display==\"block\")\n");
fprintf(fp0,"                objBranch.display=\"none\";\n");
fprintf(fp0,"        else\n");
fprintf(fp0,"                objBranch.display=\"block\";\n");
fprintf(fp0,"        swapFolder('I' + branch);\n");
fprintf(fp0,"}\n");
fprintf(fp0,"\n");
fprintf(fp0,"function swapFolder(img){\n");
fprintf(fp0,"        objImg = document.getElementById(img);\n");
fprintf(fp0,"        if(objImg.src.indexOf('plus.gif')>-1)\n");
fprintf(fp0,"                objImg.src = openImg.src;\n");
fprintf(fp0,"        else\n");
fprintf(fp0,"                objImg.src = closedImg.src;\n");
fprintf(fp0,"}\n");
fprintf(fp0,"</script>\n");
fprintf(fp0,"</head>\n");
// End of headers

   // This works in an ordered list (only to know if there is more than one processor)
   no_process=((--l.end())->processor-(l.begin())->processor)+1;
   if(no_process != 1)                                   // More than one processor
      fprintf(fp0,"<body onLoad=\"openFirst();\">\n");   // Show the 1st node open
   else
      fprintf(fp0,"<body>\n");
   fprintf(fp0,"<div class=\"title\">Conquest Execution Tree</div>\n");
   fprintf(fp0,"<br>\n");
   fprintf(fp0,"   \n");

   // Now, the records, in tree form
   it=l.begin();
   labelcount=1;
   openspans=0;
   while(it != --l.end()) 
    {
      level1=it->level;
      ++it;
      level2=it->level;
      --it;
      if(!strcmp("TOTAL",it->tag))
       {
          strcpy(class1,"total");
          strcpy(class2,"calculated");
          strcpy(image,"sigma.gif");
       }
      else if(!strcmp("ALL THE REST",it->tag))
       {
          strcpy(class1,"rest");
          strcpy(class2,"calculated");
          strcpy(image,"delta.gif");
       }
      else if(!strncmp("PROCESSOR",it->tag,9))
       {
          strcpy(class1,"process");
          strcpy(class2,"empty");
       }
      else if(!strcmp("TRUNCATED",it->tag))
       {
          strcpy(class1,"trunc");
          strcpy(class2,"calculated");
       }
      else
       {
          strcpy(class1,"tag");
          strcpy(class2,"time");
       }
      
      if(level1 < level2)
       {
          // Open a <span> block
          sprintf(label,"label%010d",labelcount);
          strcpy(branch,"branch");
          if(strcmp("TOTAL",it->tag) && strcmp("ALL THE REST",it->tag))
           {
             strcpy(image,"plus.gif");
             if(no_process > 1 && labelcount==1)
              {
                strcpy(image,"minus.gif");
                strcpy(branch,"branch1");
              }
           }  
          ++openspans;
          ++labelcount;
          fprintf(fp0,"<span class=\"trigger\" onClick=\" showBranch('%s');\">\n",label);
          fprintf(fp0,"<table class=\"timer\"><tr>\n");
          fprintf(fp0,"<td class=\"%s\"><img src=\"%s\" id=\"I%s\" alt=\"\">%s</td>\n",
                  class1,image,label,it->tag);
          if(level1 == 0)  fprintf(fp0,"<td class=\"%s\">&nbsp;</td>",class2);       // "Processor" line
          else             fprintf(fp0,"<td class=\"%s\">%12.6f s</td>",class2,it->time);
          fprintf(fp0,"</tr></table>\n");
          fprintf(fp0,"</span>\n");
          fprintf(fp0,"<span class=\"%s\" id=\"%s\">\n",branch,label);
          ++it;
       }
      else
       {
          if(strcmp("TOTAL",it->tag) && strcmp("ALL THE REST",it->tag))
           {
              strcpy(image,"blank.gif");
           }
          fprintf(fp0,"<table class=\"timer\"><tr>\n");
          fprintf(fp0,"<td class=\"%s\"><img src=\"%s\" alt=\"\">%s</td>\n",
                  class1,image,it->tag);
          fprintf(fp0,"<td class=\"%s\">%12.6f s</td>",class2,it->time);
          fprintf(fp0,"</tr></table>\n");
          if(level1 > level2)      // Close span blocks
           {
              for(i=0; i < (level1-level2); ++i)
               {
                  fprintf(fp0,"</span>\n");
                  --openspans;
               }
           } 
          ++it;
       }
    }

   // Do the last one: It should be a "Total" register, but check
   if(!strcmp("TOTAL",it->tag))
    {
       fprintf(fp0,"<table class=\"timer\"><tr>\n");
       fprintf(fp0,"<td class=\"total\"><img src=\"sigma.gif\" alt=\"\">%s</td>\n",it->tag);
       fprintf(fp0,"<td class=\"calculated\">%12.6f s</td>",it->time);
       fprintf(fp0,"</tr></table>\n");
    }
   else
    {
       fprintf(fp0,"<table class=\"timer\"><tr>\n");
       fprintf(fp0,"<td class=\"tag\"><img src=\"blank.gif\" alt=\"\">%s</td>\n",it->tag);
       fprintf(fp0,"<td class=\"time\">%12.6f s</td>",it->time);
       fprintf(fp0,"</tr></table>\n");
    }
   while(openspans)
    {
       fprintf(fp0,"</span>\n");
       --openspans;
    }

   // Finally, close the file
   fprintf(fp0,"</body>\n");
   fprintf(fp0,"</html>\n");
   fclose(fp0);
 }

