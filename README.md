
# Bela adaptation of the Dynamic Stiff String Simulation


## Introduction
This is a Bela project for simulating a stiff string in real time.

This project uses C++14 features. You will need to add

    CPPFLAGS=-std=c++14

to the "Make Parameters" section in Bela IDE.

## Acknowledgements
The simulation is based on the Juce plugin by Silvin Willemsen [github.com/SilvinWillemsen/RealTimeDynamic](https://github.com/SilvinWillemsen/RealTimeDynamic) *[1]*.

[1] Silvin Willemsen and Stefania Serafin, “REAL-TIME IMPLEMENTATION OF THE DYNAMIC STIFF STRING USING
FINITE-DIFFERENCE TIME-DOMAIN METHODS AND THE DYNAMIC GRID,”Proceedings of the 25th International Conference on Digital Audio Effects (DAFx20in22), Vienna, Austria, September 6-10, 2022.

## Final To Do's before Testing

In order of relevance:

### DigitalWrite LED handling. ESSENTIAL

### AnalogOut Parameter values. COOL





## DONE

### Fix Calibration. DOOOOONEEEE

### Introduce 1V/oct behaviour at Lenght Input. DONE

This is of the type:
- minimumLengthValue * 2 ^ (incomingvalue)

However, it will require Length to be mapped in reverse:
- from maximum length to minimum lenght (low pitch to high pitch, since that is how pitch voltage works).

### Fix handling of sigmas at De-clipping algorithm. DONE

### Introduce E handling stage NICE LITTLE DETAIL // UNNECESARY // RANGE DIMINISHED -> UNPERCEPTIBLE

need to introduce a de-clipping stage similar to sigmas since 

when L, rho, r and T are in extreme values (at the same time) 
and 
E <= a low value

it creates underruns.

So we need to introduce something like:

if (L >= w && rho >= x && r <= y && T <= z)
{
E lower range is a bigger number
}

Exact values need calibrated by ear.

It is optional since I could just increase the lower range of T (it only started happening when ranges were widened).

### Add Spray Function to loc. DOUNEEE. Also TRIGGER RANDOMIZER OMG

- while Button "(any button not in use)" = true
- loc = loc + spread * input[5]
- spread = random between -1 and 1

advantage of this is:

Since the CVinput for loc is lost, it is good to have an extra parameter that randomises loc slightly.

It is also a reference to the "spray parameter in Instruo's Arhbhar: explained in https://youtu.be/hw73DlxVWrI?t=500
