// ------------------------------------------------------------------------------
// $Id$
// ------------------------------------------------------------------------------
// File cq_ut_general_cTimerRegister.hpp
// ------------------------------------------------------------------------------

// ****h* Conquest/utilities/cq_ut_general_cTimerRegister.hpp *
//
// NAME
//  cq_ut_general_cTimerRegister.hpp
// PURPOSE
//  Declaration of the cTimerRegister class
// USES
//  cq_ut_general_cTimerRegister.hpp
// AUTHOR
//  Antonio S. Torralba
// CREATION DATE
//  2008/07/28
// MODIFICATION HISTORY
//
// ****

#ifndef _CQ_UT_GENERAL_TIMER_REGISTER
#define _CQ_UT_GENERAL_TIMER_REGISTER

#define MAX_TOK 500
#define MAX_LIN 2048

using namespace std;

// -----------------------------------------------------------------------------
// Class cTimerList
// -----------------------------------------------------------------------------

// ****f* cq_ut_general_cTimerRegister.hpp/cTimerRegister *
//
// NAME
//  cTimerRegister
// PURPOSE
//  Class declaration
// INPUTS
//
// AUTHOR
//  Antonio S. Torralba
// CREATION DATE
//  2008/07/28
// MODIFICATION HISTORY
//
// ****

class cTimerRegister {
    int tokenize(const char *, char **);

public:
    int processor;               // For the moment, easily accesible
    int level;
    char tag[MAX_LIN+1];
    double time;

    cTimerRegister(const char *);                    // From a line from Conquest
    cTimerRegister(int, int, const char *, double);  // From actual data
    bool operator==(int);                            // Required for finds by level
    bool operator==(long int);                       // Required for finds by processor
    bool operator!=(int);                            // Required for finds by level
    bool operator!=(long int);                       // Required for finds by processor
    void print();
 };

#endif
