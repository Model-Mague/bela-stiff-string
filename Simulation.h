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

	bool isTriggerPressed() { return m_buttonPressed[(size_t)Button::TRIGGER]; }
	bool isModePressed() { return m_buttonPressed[(size_t)Button::MODE]; }
	bool isUpPressed() { return m_buttonPressed[(size_t)Button::UP]; }
	bool isDownPressed() { return m_buttonPressed[(size_t)Button::DOWN]; }

private:
	int m_audioFramesPerAnalogFrame;
	float m_inverseSampleRate;

	float m_analogIn[analogInputs] = {};
	float m_lfo[analogInputs] = {};
	float m_phase[analogInputs] = {};

	float m_amplitude;
	float m_frequency;

	unsigned long long frame;

	bool m_buttonPressed[4] = {}; // Buttons pressed at this time
	uint64_t m_buttonLastActivated[4] = {}; // Timestamp of when the button was last activated
};
