#pragma once
#include "BelaMock.h"
#include "DynamicStiffString/DynamicStiffString.h"

class Simulation {
public:
	static constexpr int analogInputs = 8;

	enum class Button : size_t {
		TRIGGER = 0,
		MODE = 1,
		UP = 2,
		DOWN = 3
	};

	Simulation(BelaContext* context);
	void readInputs(BelaContext* context, int nFrame);
	void writeOutputs(BelaContext* context, int nFrame);
	const float* getAnalogIn() { return m_analogIn; }
	const float* getLfo() { return m_lfo; }
	const float* getPhase() { return m_phase; }
	const int getAnalogInputCount() { return analogInputs; }

	void update(BelaContext* context, std::shared_ptr<DynamicStiffString> pDynamicStiffString); // Runs every audio frame

	bool isButtonReleased(const Button b) { return m_buttonPreviousState[(size_t)b] != 0 && m_buttonState[(size_t)b] == 0; }

private:
	int m_audioFramesPerAnalogFrame;
	float m_inverseSampleRate;

	float m_analogIn[analogInputs] = {};
	float m_lfo[analogInputs] = {};
	float m_phase[analogInputs] = {};

	float m_amplitude;
	float m_frequency;

	unsigned long long frame;

	int m_buttonPreviousState[4] = {}; // Last button state
	int m_buttonState[4] = {}; // Current button state
	//uint64_t m_buttonLastActivated[4] = {}; // Timestamp of when the button was last activated
};
