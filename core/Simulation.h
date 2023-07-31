#pragma once

#ifdef DESKTOP_BUILD
#include "BelaMock.h"
#else
#include "Bela.h"
#endif

#include "../components/AnalogInput.h"
#include "../components/Button.h"
#include "../components/LEDScreen.h"
#include "AudioBuffer.h"
#include "DynamicStiffString/DynamicStiffString.h"
#include "Parameters.h"

#include <array>
#include <map>
#include <memory>
#include <set>
#include <vector>
#include <fstream>

#include "Compressor/SimpleComp.h"

class Simulation {
public:
	static constexpr int sAnalogInputCount = 8;
	static constexpr short sDSSUpdateRate = 20; // DSS should be updated no more frequently than every 20 audio frames

	Simulation(BelaContext* context);
	void readInputs(BelaContext* context, int frame);
	void writeOutputs(BelaContext* context, int frame);
	void writeAudio(BelaContext* context, int frame);

	const float* getLfo() { return m_lfo; }
	const float* getPhase() { return m_phase; }
	const int getAnalogInputCount() { return sAnalogInputCount; }

	chunkware_simple::SimpleComp Compressor;

	//LOC & SPRAY VALUES
	float nonsprayloc; // pot 5 when button not pressed
	float sprayAmount; // pot 5 when button is pressed
	float sprayValue;  // randomised -0.5 -> 0.5 * sprayAmount
	float sprayedloc;  // loc position + sprayValue (sort of)


	void update(BelaContext* context); // Runs every audio frame

	std::string getCalibrationResults();

	//Diagnostics function to measure highest value before compression

	double maxValue = 0;

	void maxVal(double value)
	{
		if (value > maxValue)
		{
			maxValue = value;

			if (maxValues.is_open()) 
			{
				maxValues << maxValue << "\n";
			}
			else 
			std::cout << "Failed to open maxValues.txt" << std::endl;
		}
	}

	std::fstream maxValues;


private:
	// Screen
	LEDScreen m_screen;

	// DSS simulation
	std::unique_ptr<DynamicStiffString> m_pDynamicStiffString;

	// Maintains an internal state of all the simulation parameters
	// These may be adjusted internally and don't necessarily match up with inputs
	Parameters m_parameters;

	// Mapping of button type to button object
	std::map<Button::Type, Button> m_buttons;

	// When we encounter a change in inputs, we insert the channel # in here
	// Then this set is consumed in the update function
	std::set<ParameterName> m_parametersToUpdate;

	// Circular buffer for receiving audio input
	AudioBuffer m_audioBuffer;

	// Excitation functions. Take the current index and a total length and returns a float
	// Raised cosine
	std::function<float(int, int)> m_fnRaisedCos;
	// Sample based excitation
	std::function<float(int, int)> m_fnSampleExcitation;

	int m_audioFramesPerAnalogFrame;
	float m_inverseSampleRate;

	float m_lfo[sAnalogInputCount] = {};
	float m_phase[sAnalogInputCount] = {};

	float m_amplitude;
	float m_frequency;

	// Counter for ensuring no-more-than-every-20-frames update frequency
	short m_updateFrameCounter = 0;

	// Set if we're to actively damp the signal in the next update call
	bool clippingFlag = false;
	bool hasCorrectedFlag = false;
	bool stableFlag = true;
	float correctionValue; // Damping proportion

};
