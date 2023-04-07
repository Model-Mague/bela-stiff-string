
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

### Fix handling of sigmas at De-clipping algorithm.
### DigitalWrite LED handling.
### AnalogOut Parameter values.
### Introduce 1V/oct behaviour at Lenght Input.

This is of the type:
- minimumLengthValue * 2 ^ (incomingvalue)

However, it will require Length to be mapped in reverse:
- from maximum length to minimum lenght (low pitch to high pitch, since that is how pitch voltage works).

#### Add Spray Function to loc.

- while Button "(any button not in use)" = true
- loc = loc + spread * input[5]
- spread = random between -1 and 1

advantage of this is:

Since the CVinput for loc is lost, it is good to have an extra parameter that randomises loc slightly.

It is also a reference to the "spray parameter in Instruo's Arhbhar: explained in https://youtu.be/hw73DlxVWrI?t=500

