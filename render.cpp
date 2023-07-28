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
#include "core/Simulation.h"

#include <memory>
#include <iostream>

// Global variables
Scope scope;
std::unique_ptr<Simulation> pSimulation;


bool setup(BelaContext* context, void* userData)
{
	scope.setup(Simulation::sAnalogInputCount, context->audioSampleRate);

	pSimulation = std::make_unique<Simulation>(context);

	//Initialize Compressor parameters.
	pSimulation->Compressor.setSampleRate(context->audioSampleRate);
	pSimulation->Compressor.setThresh(60);
	pSimulation->Compressor.setRatio(1);
	pSimulation->Compressor.setRelease(100);
	pSimulation->Compressor.setAttack(100);

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
		pSimulation->update(context);

		// 3. Write outputs
		pSimulation->writeOutputs(context, n);

		// 4. Write out audio
		pSimulation->writeAudio(context, n);
	}
}

void cleanup(BelaContext* context, void* userData)
{
	std::cout << pSimulation->getCalibrationResults();
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
