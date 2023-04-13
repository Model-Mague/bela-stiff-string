
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

## Improvements.

### AnalogOut Parameter values.
### String Reverb.

If(AudioIn != 0)
    {audioinFlag = true;)

if(audioinFlag)
    exciteFlag = true and quickly false // constant excite
    raisedCosineFunction = AudioIn // Basically, if there is any audio coming in, this should substitute the raised cosine operation.
    
### ReadOut.

The points at which the variable "output" reads the information. These could be user definable to widen the stereo feeling.

### Several Others recommended by Instruo Staff

Will come up with full list after testing is finished.

## Model Extension.

### Exciter.

Substitute Raised Cosine function by

    - Struck (with a hammer - like in a piano)
    - Pluck (with a pick or finger - like in a guitar or harp)
    - Bow (with a violin bow)
    
These will require research and further searching but there are tons of papers about it from S. Bilbao, S. Willensen, etc...

### Bounds Info.

Bounds can be input in the model, as to say, limits. The same way that the sound is different at the edges, different several bounds
can be imposed on the system. This way you can create guitar frets.

### Several Strings/Poliphony.

On Numerical Sound Synthesis, S. Bilbao explains a chapter on how different strings interact between each other. This could be added into the model.


## The NESS project.

S. Bilboa studied several of these models under the NESS project at Edinburgh University. He managed to get all these working (not real time).
Anyhow, they have got very short videos that are very interesting explaining these things that can be taken as a reference.


- Guitar: https://www.youtube.com/watch?v=NTp7dRuld08 (pluck exciting mechanism + bounds + several strings)
- Bowed String: https://www.youtube.com/watch?v=fQMJm-YMXuQ (bowed exciting mechanism + bounds + several strings)




