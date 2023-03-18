/*
 ____  _____ _        _
| __ )| ____| |      / \
|  _ \|  _| | |     / _ \
| |_) | |___| |___ / ___ \
|____/|_____|_____/_/   \_\

In render() you'll see a nested for loop structure. You'll see this in all Bela projects.
The first for loop cycles through 'audioFrames', the second through 'audioChannels' (in this case left 0 and right 1).
------------------------------------

Dynamic Stiff String implementation based on DAFx 2022 paper submission https://dafx2020.mdw.ac.at/proceedings/papers/DAFx20in22_paper_11.pdf
"Real-Time Implementation of the Dynamic Stiff String using Finite-Difference Time-Domain Methods and the Dynamic Grid" by Silvin Willemsen and Stefania Serafin.

*/
#define _USE_MATH_DEFINES
#include <cmath>
#include <fstream>
#include <vector>
#include <iostream>
#include <vector>
#include <array>


// @HERE Uncomment this to run project on desktop
#define DESKTOP_BUILD

#ifndef DESKTOP_BUILD
#include <Bela.h>
#endif


#ifdef DESKTOP_BUILD
bool setup()
#else
bool setup(BelaContext* context, void* userData)
#endif
{
	return true;
}


#ifdef DESKTOP_BUILD
int main(int argc, char** argv)
#else
void render(BelaContext* context, void* userData)
#endif
{
#ifdef DESKTOP_BUILD
	setup();
#endif

#ifdef DESKTOP_BUILD
	for (unsigned int n = 0; n < 20; n++)
	{
#else
	for (unsigned int n = 0; n < context->audioFrames; n++)
	{
#endif

#ifndef DESKTOP_BUILD
		for (unsigned int channel = 0; channel < context->audioOutChannels; channel++)
		{
			audioWrite(context, n, channel, 0);
		}
#endif
	}
}

#ifndef DESKTOP_BUILD
void cleanup(BelaContext* context, void* userData)
{

}
#endif
