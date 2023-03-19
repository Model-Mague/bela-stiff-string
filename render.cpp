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

// @HERE Uncomment this to run project on desktop
//#define DESKTOP_BUILD

#ifndef DESKTOP_BUILD
#include <Bela.h>
#endif

#include "DynamicStiffString/DynamicStiffString.h"

#include <memory>
#include <iostream>


// Global variables
std::unique_ptr<DynamicStiffString> pDynamicStiffString;


#ifdef DESKTOP_BUILD
bool setup()
#else
bool setup(BelaContext* context, void* userData)
#endif
{
	DynamicStiffString::SimulationParameters parameters = {};
	parameters.L = 1.0f;
	parameters.rho = 7850.0f;
	parameters.r = 0.0005f;
	parameters.T = 300.0f;
	parameters.E = 2e11f;
	parameters.sigma0 = 1.0f;
	parameters.sigma1 = 0.005f;

#ifdef DESKTOP_BUILD
	float sampleRate = 44100.0;
#else
	float sampleRate = context->audioSampleRate;
#endif

	pDynamicStiffString = std::make_unique<DynamicStiffString>(parameters, 1.0f / sampleRate);
	pDynamicStiffString->excite();

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
		// @TODO: Handle parameter changes here
		// pDynamicStiffString->refreshParameter(paramId, value);

		pDynamicStiffString->refreshCoefficients();
		pDynamicStiffString->calculateScheme();
		pDynamicStiffString->updateStates();

		float output = (float)pDynamicStiffString->getOutput();

#ifdef DESKTOP_BUILD
		std::cout << output << std::endl;
#else
		for (unsigned int channel = 0; channel < context->audioOutChannels; channel++)
		{
			audioWrite(context, n, channel, Global::limit(output, -1.0, 1.0));
		}
#endif
	}
}

#ifndef DESKTOP_BUILD
void cleanup(BelaContext* context, void* userData)
{

}
#endif
