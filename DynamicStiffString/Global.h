/*
  ==============================================================================

    Global.h
    Created: 10 Feb 2022 9:39:22am
    Author:  Silvin Willemsen

  ==============================================================================
*/

#pragma once
//#define RECORD

#include <random>

namespace Global {
#ifdef RECORD
    static int samplesToRecord = 500;

#endif
    static int margin = 10;
    static float NmaxChange = 1.0f / 20.0f;
    static float sig1min = 0.0002f;
    static float boundaryEllRad = 10.0f;
    // limiter
    static float limit(float val, float min, float max)
    {
        if (val < min)
        {
            val = min;
            return val;
        }
        else if (val > max)
        {
            val = max;
            return val;
        }
        return val;
    }

    static float f_random(float min, float max)
    {
        return min + static_cast <float> (rand()) / (static_cast <float> (RAND_MAX / (max - min)));
    }
}
