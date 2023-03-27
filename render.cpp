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

@TODO list:
1. Triggers to work (digital output)
2. Recalculate parameters based on analog reads from channels 0-7
3. Digital write to screen of currently-adjusted input; brighter light if input closer to max
4. Potentially all params might need to get locked after a button is pressed // possibly multiple modes based on button


*/
#define _USE_MATH_DEFINES

// The DESKTOP_BUILD preprocessor definition will be added if compiling with VS
#ifndef DESKTOP_BUILD
#include <Bela.h>
#include <libraries/Scope/Scope.h>
#else
#include "BelaMock.h"
#include <chrono>
#endif

#include "DynamicStiffString/DynamicStiffString.h"
#include "Simulation.h"

#include <memory>
#include <iostream>

// Global variables
Scope scope;
std::shared_ptr<DynamicStiffString> pDynamicStiffString;
std::unique_ptr<Simulation> pSimulation;


bool setup(BelaContext* context, void* userData)
{
	scope.setup(Simulation::analogInputs, context->audioSampleRate);

	DynamicStiffString::SimulationParameters parameters = {};
	parameters.L = 1.0f;
	parameters.rho = 7850.0f;
	parameters.r = 0.0005f;
	parameters.T = 300.0f;
	parameters.E = 2e11f;
	parameters.sigma0 = 1.0f;
	parameters.sigma1 = 0.005f;

	float inverseSampleRate = 1.0f / context->audioSampleRate;
	pDynamicStiffString = std::make_shared<DynamicStiffString>(parameters, inverseSampleRate);
	pSimulation = std::make_unique<Simulation>(context);

	return true;
}

void render(BelaContext* context, void* userData)
{
	for (unsigned int n = 0; n < context->audioFrames; n++)
	{
#ifdef DESKTOP_BUILD
		context->audioFramesElapsed++; // Bela does this automatically
#endif

		// 1. Read inputs
		pSimulation->readInputs(context, n);

		// 2. Update calculations if needed
		pSimulation->update(context, pDynamicStiffString);

		// 3. Write outputs
		pSimulation->writeOutputs(context, n);

		// 4. Write out audio
		float output = (float)pDynamicStiffString->getOutput();
		for (unsigned int channel = 0; channel < context->audioOutChannels; channel++)
		{
			audioWrite(context, n, channel, Global::limit(output, -1.0, 1.0));
		}

		scope.log(pSimulation->getAnalogIn());
	}
}

void cleanup(BelaContext* context, void* userData)
{

}

#ifdef DESKTOP_BUILD
int main(int argc, char** argv)
{
	std::unique_ptr<BelaContext> context = std::make_unique<BelaContext>();
	setup(context.get(), nullptr);

	// Max time we can spend producing a frame
	const float secondsPerFrame = 1000.f / context->audioSampleRate;
	const float microsecondsPerFrame = secondsPerFrame * 10000.f;

	auto lastUpdate = std::chrono::system_clock::now();

	// This loop is meant to run no faster than the real time constraint
	while (true)
	{
		auto currentTime = std::chrono::system_clock::now();
		auto deltaT = std::chrono::duration_cast<std::chrono::microseconds>(currentTime - lastUpdate).count();

		// At 44.1kHz we need to produce 1 frame every 0.22ms
		if (deltaT >= microsecondsPerFrame)
		{
			render(context.get(), nullptr);
			lastUpdate = currentTime;
		}
	}
	cleanup(context.get(), nullptr);
	return 0;
}
#endif
