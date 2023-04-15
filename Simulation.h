#pragma once

#ifdef DESKTOP_BUILD
#include "BelaMock.h"
#else
#include "Bela.h"
#endif

#include "AnalogInput.h"
#include "LEDScreen.h"
#include "DynamicStiffString/DynamicStiffString.h"

#include <array>
#include <map>
#include <memory>
#include <set>
#include <vector>

class Simulation {
public:
	static constexpr int sAnalogInputCount = 8;
	static constexpr short sDSSUpdateRate = 20; // DSS should be updated no more frequently than every 20 audio frames

	enum class Button : size_t {
		TRIGGER = 0,
		MODE = 1,
		SPRAY = 2,
		DOWN = 3
	};

	Simulation(BelaContext* context);
	void readInputs(BelaContext* context, int frame);
	void writeOutputs(BelaContext* context, int frame);
	void writeAudio(BelaContext* context, int frame);

	const float* getLfo() { return m_lfo; }
	const float* getPhase() { return m_phase; }
	const int getAnalogInputCount() { return sAnalogInputCount; }
	
	//SPRAY VALUES
	float sprayAmount; // pot 5 when button is pressed
	float sprayValue;  // randomised -0.5 -> 0.5 * sprayAmount
	float sprayedloc;  // loc position + sprayValue (sort of)
	
	
	void update(BelaContext* context); // Runs every audio frame

	bool isButtonReleased(const Button b) { return m_buttonPreviousState[(size_t)b] != 0 && m_buttonState[(size_t)b] == 0; }

	std::string getCalibrationResults();

private:
	// Screen
	LEDScreen m_screen;

	// DSS simulation
	std::unique_ptr<DynamicStiffString> m_pDynamicStiffString;

	// Maintains an internal state of all the simulation parameters
	// These may be adjusted internally and don't necessarily match up with inputs
	std::map<std::string, float> m_parameters;

	// Mapping of parameter name to parameter ID in DSS simulation
	std::map<std::string, int> m_parameterIdMap;

	// When we encounter a change in inputs, we insert the channel # in here
	// Then this set is consumed in the update function
	std::set<int> m_channelsToUpdate;

	// Mapping of label to AnalogInput (e.g. sigma0 -> AnalogInput for 6th channel) 
	std::map<std::string, int> m_labelToAnalogIn;

	int m_audioFramesPerAnalogFrame;
	float m_inverseSampleRate;

	float m_lfo[sAnalogInputCount] = {};
	float m_phase[sAnalogInputCount] = {};

	float m_amplitude;
	float m_frequency;

	// Vector of classes that allow us to read from an analog channels (0~7)
	std::vector<AnalogInput> m_analogInputs;

	// Counter for ensuring no-more-than-every-20-frames update frequency
	short m_updateFrameCounter = 0;
	
	// Set if we're to actively damp the signal in the next update call
	bool clippingFlag = false;
	float correctionValue; // Damping proportion

	int m_buttonPreviousState[4] = {}; // Last button state
	int m_buttonState[4] = {}; // Current button state
};
