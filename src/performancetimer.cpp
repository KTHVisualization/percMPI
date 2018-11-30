/*  Time-stamp: <@(#)PerformanceTimer.cpp   01/01/2010 4:08 PM   Tino Weinkauf>
 *
 *********************************************************************
 *  File    : PerformanceTimer.cpp
 *
 *  Project : Simple Mesh Viewer
 *
 *  Package : Utils
 *
 *  Company : NYU
 *
 *  Author  : Tino Weinkauf                           Date: 01/01/2010
 *
 *  Purpose : Implementation of methods for class PerformanceTimer
 *
 *********************************************************************
 * Version History:
 *
 * V 0.10  01/01/2010 16:08:24  TW : First Revision
 *
 **********************************************************************/

#include "performancetimer.h"

namespace perc {

PerformanceTimer::PerformanceTimer() {
#ifdef WIN32
    QueryPerformanceFrequency((LARGE_INTEGER*)&Frequency);
#else
    Frequency = CLOCKS_PER_SEC;
#endif

    Reset();
}

PerformanceTimer::~PerformanceTimer() {}

void PerformanceTimer::Reset() {
#ifdef WIN32
    QueryPerformanceCounter((LARGE_INTEGER*)&Start);
#else
    Start = clock();
#endif
}

float PerformanceTimer::ElapsedTime() const {
#ifdef WIN32
    __int64 End;
    QueryPerformanceCounter((LARGE_INTEGER*)&End);
#else
    clock_t End = clock();
#endif

    return (float)(End - Start) / Frequency;
}
}  // namespace perc