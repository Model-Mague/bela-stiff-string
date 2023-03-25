#pragma once
#include "BelaMock.h"

class Simulation {
public:
	static constexpr int analogInputs = 8;

	Simulation(BelaContext* context);
	void readInputs(BelaContext* context, int nFrame);
	void writeOutputs(BelaContext* context, int nFrame);
	const float* getAnalogIn() { return m_analogIn; }
	const float* getLfo() { return m_lfo; }
	const float* getPhase() { return m_phase; }
	const int getAnalogInputCount() { return analogInputs; }

private:
	int m_audioFramesPerAnalogFrame;
	float m_inverseSampleRate;

	float m_analogIn[analogInputs] = {};
	float m_lfo[analogInputs] = {};
	float m_phase[analogInputs] = {};

	float m_amplitude;
	float m_frequency;
};
