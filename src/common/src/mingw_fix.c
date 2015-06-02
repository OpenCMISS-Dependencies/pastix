/*
 * Includes a hotfix where clock_gettime is not implemented on mingw 32bit environments.
 */
#if __MINGW32__
#include "common_pastix.h"
#include <time.h>
#include "windows.h"
#ifndef CLOCK_REALTIME
	#define CLOCK_REALTIME 1
#endif

    LARGE_INTEGER
    getFILETIMEoffset()
    {
        SYSTEMTIME s;
        FILETIME f;
        LARGE_INTEGER t;

        s.wYear = 1970;
        s.wMonth = 1;
        s.wDay = 1;
        s.wHour = 0;
        s.wMinute = 0;
        s.wSecond = 0;
        s.wMilliseconds = 0;
        SystemTimeToFileTime(&s, &f);
        t.QuadPart = f.dwHighDateTime;
        t.QuadPart <<= 32;
        t.QuadPart |= f.dwLowDateTime;
        return (t);
    }

    int
    clock_gettime(int X, struct timespec *tv)
    {
        LARGE_INTEGER           t;
        FILETIME            f;
        double                  microseconds;
        static LARGE_INTEGER    offset;
        static double           frequencyToMicroseconds;
        static int              initialized = 0;
        static BOOL             usePerformanceCounter = 0;

        if (!initialized) {
            LARGE_INTEGER performanceFrequency;
            initialized = 1;
            usePerformanceCounter = QueryPerformanceFrequency(&performanceFrequency);
            if (usePerformanceCounter) {
                QueryPerformanceCounter(&offset);
                frequencyToMicroseconds = (double)performanceFrequency.QuadPart / 1000000.;
            } else {
                offset = getFILETIMEoffset();
                frequencyToMicroseconds = 10.;
            }
        }
        if (usePerformanceCounter) QueryPerformanceCounter(&t);
        else {
            GetSystemTimeAsFileTime(&f);
            t.QuadPart = f.dwHighDateTime;
            t.QuadPart <<= 32;
            t.QuadPart |= f.dwLowDateTime;
        }

        t.QuadPart -= offset.QuadPart;
        microseconds = (double)t.QuadPart / frequencyToMicroseconds;
        t.QuadPart = microseconds;
        tv->tv_sec = t.QuadPart / 1000000;
        tv->tv_nsec = t.QuadPart % 1000000;
        return (0);
    }
#endif