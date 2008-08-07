// ------------------------------------------------------------------------------
// $Id$
// ------------------------------------------------------------------------------
// File cq_ut_general_cTimerList.hpp
// ------------------------------------------------------------------------------

// ****h* Conquest/utilities/cq_ut_general_cTimerList.hpp *
//
// NAME
//  cq_ut_general_cTimerList.hpp
// PURPOSE
//  Declaration of the cTimerList class
// USES
//  Standard Template Library: list
//  cq_ut_general_cTimerRegister.hpp
// AUTHOR
//  Antonio S. Torralba
// CREATION DATE
//  2008/07/28
// MODIFICATION HISTORY
//
// ****


#ifndef _CQ_UT_TIMER_LIST
#define _CQ_UT_TIMER_LIST

#include <list>

#include "cq_ut_general_cTimerRegister.hpp"

using namespace std;

// -----------------------------------------------------------------------------
// Class cTimerList
// -----------------------------------------------------------------------------

// ****f* cq_ut_general_cTimerList.hpp/cTimerList *
//
// NAME
//  cTimerList
// PURPOSE
//  Class declaration: For the moment built by composition, not inheritance
// INPUTS
//
// AUTHOR
//  Antonio S. Torralba
// CREATION DATE
//  2008/07/28
// MODIFICATION HISTORY
//
// ****

class cTimerList {
      list<cTimerRegister> l;             // The list itself

      static bool comparePerProcessor(cTimerRegister, cTimerRegister);
      int  orderByProcessor();
      void orderExecution(list<cTimerRegister> :: iterator,list<cTimerRegister> :: iterator,int);
        // This expects execution order already
      double assignSubtotals(list<cTimerRegister> :: iterator, list<cTimerRegister> :: iterator, int);
        // This function is the outer level, so it's the one to be called
      void calculateSubtotals(list<cTimerRegister> :: iterator,list<cTimerRegister> :: iterator);
      void fixList(list<cTimerRegister> :: iterator,list<cTimerRegister> :: iterator);

 public:
      cTimerList(const char *);
        // This should be done (more mantainable) using classes and configuration files
      void printHtml(const char *);
 };

#endif

