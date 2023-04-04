#pragma once

#ifdef DESKTOP_BUILD
#include "BelaMock.h"
#else
#include "Bela.h"
#endif

#include "AnalogInput.h"
#include "DynamicStiffString/DynamicStiffString.h"

#include <array>
#include <memory>
#include <set>

class Simulation {
public:
	static constexpr int sAnalogInputCount = 8;
	static constexpr short sDSSUpdateRate = 20; // DSS should be updated no more frequently than every 20 audio frames

	enum class Button : size_t {
		TRIGGER = 0,
		MODE = 1,
		UP = 2,
		DOWN = 3
	};

	Simulation(BelaContext* context);
	void readInputs(BelaContext* context, int frame);
	void writeOutputs(BelaContext* context, int frame);
	void writeAudio(BelaContext* context, int frame);

	const float* getAnalogIn() { return m_rangeMappedInputs; }
	const float* getLfo() { return m_lfo; }
	const float* getPhase() { return m_phase; }
	const int getAnalogInputCount() { return sAnalogInputCount; }

	void update(BelaContext* context); // Runs every audio frame

	bool isButtonReleased(const Button b) { return m_buttonPreviousState[(size_t)b] != 0 && m_buttonState[(size_t)b] == 0; }

	std::string getCalibrationResults();

private:
	std::unique_ptr<DynamicStiffString> m_pDynamicStiffString;
	float m_excitationLoc;
	std::set<int> m_channelsToUpdate;

	int m_audioFramesPerAnalogFrame;
	float m_inverseSampleRate;

	float m_rangeMappedInputs[sAnalogInputCount] = {};
	float m_lfo[sAnalogInputCount] = {};
	float m_phase[sAnalogInputCount] = {};

	float m_amplitude;
	float m_frequency;

	std::vector<AnalogInput> m_analogInputs;

	short m_updateFrameCounter = 0;

	int m_buttonPreviousState[4] = {}; // Last button state
	int m_buttonState[4] = {}; // Current button state
	//uint64_t m_buttonLastActivated[4] = {}; // Timestamp of when the button was last activated
};
