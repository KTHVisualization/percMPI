/*  Copied from KTH Inviwo
 *  Time-stamp: <@(#)PerformanceTimer.h   01/01/2010 4:08 PM   Tino Weinkauf>
 *
 *********************************************************************
 *  File    : PerformanceTimer.h
 *
 *  Project : Simple Mesh Viewer
 *
 *  Package : Utils
 *
 *  Company : NYU
 *
 *  Author  : Tino Weinkauf
 *
 *  Date    : 01/01/2010
 *
 *  Purpose : Declaration of class PerformanceTimer
 *
 *********************************************************************
 * Version History:
 *
 * V 0.10  01/01/2010 16:08:15  TW : First Revision
 *
 *********************************************************************
 */

#ifndef UTILS_PERFORMANCETIMER_H
#define UTILS_PERFORMANCETIMER_H

#ifdef WIN32
#include <windows.h>
#else
#include <time.h>
#endif

namespace perc {

/** Provides a simple timer to measure code performance.
@author Tino Weinkauf
*/
class PerformanceTimer {
    // Friends
    // Types
public:
    // Construction / Deconstruction
public:
    PerformanceTimer();
    virtual ~PerformanceTimer();

    // Methods
public:
    /// Returns elapsed time since construction of the object in seconds
    float ElapsedTime() const;

    /// Returns elapsed time since construction of the object in seconds and resets the clock.
    inline float ElapsedTimeAndReset() {
        const float Time = ElapsedTime();
        Reset();
        return Time;
    }

    /// Resets (or re-starts) the clock
    void Reset();

    // Attributes
protected:
#ifdef WIN32
    __int64 Start;
    __int64 Frequency;
#else
    clock_t Start;
    clock_t Frequency;
#endif
};

}  // namespace perc

#endif  // UTILS_PERFORMANCETIMER_H