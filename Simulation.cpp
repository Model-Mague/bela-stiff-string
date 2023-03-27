#include "Simulation.h"

#define _USE_MATH_DEFINES
#include <cmath>

#define INTERACTION_DELAY 1000 // frames

Simulation::Simulation(BelaContext* context) : m_amplitude(5.f), m_frequency(0.1f)
{
	m_inverseSampleRate = 1.0f / context->audioSampleRate;
	if (context->analogFrames)
		m_audioFramesPerAnalogFrame = context->audioFrames / context->analogFrames;

	memset((void*)&m_buttonPreviousState, 0, 4);
	memset((void*)&m_buttonState, 0, 4);
	memset((void*)&m_analogIn, 0, analogInputs * sizeof(float));
	memset((void*)&m_lfo, 0, analogInputs * sizeof(float));
	memset((void*)&m_phase, 0, analogInputs * sizeof(float));

	m_amplitude = 5;
	m_frequency = 0.1f;
}

void Simulation::update(BelaContext* context, std::shared_ptr<DynamicStiffString> pDynamicStiffString)
{
	// 1. Handle parameter changes
	// pDynamicStiffString->refreshParameter(paramId, value);

	// if (parameterChanged) { ...
	pDynamicStiffString->refreshCoefficients();
	pDynamicStiffString->calculateScheme();
	pDynamicStiffString->updateStates();
	// }

	// 2. Handle trigger button (should probably be last)
	if (isButtonReleased(Button::TRIGGER))
	{
		pDynamicStiffString->excite();
	}

	// 3. Update phase
	for (int channel = 0; channel < analogInputs; channel++)
	{
		m_phase[channel] += 2.0f * (float)M_PI * m_frequency * m_inverseSampleRate;
		if (m_phase[channel] > M_PI)
			m_phase[channel] -= 2.0f * (float)M_PI;
	}
}

void Simulation::readInputs(BelaContext* context, int nFrame)
{
	const int frame = nFrame / m_audioFramesPerAnalogFrame;
	for (int channel = 0; channel < analogInputs; channel++)
	{
		m_analogIn[channel] = analogRead(context, frame, channel);
	}

	// Buttons, left to right
	// 15, 14, 13, 12
	// Trigger, Mode, Up, Down
	static int buttonChannel[4] = { 15, 14, 13, 12 };
	for (const auto& button : { Button::TRIGGER, Button::MODE, Button::UP, Button::DOWN })
	{
		m_buttonPreviousState[(size_t)button] = m_buttonState[(size_t)button];
		m_buttonState[(size_t)button] = digitalRead(context, nFrame, buttonChannel[(size_t)button]);
	}
}

void Simulation::writeOutputs(BelaContext* context, int nFrame)
{
	const int frame = nFrame / m_audioFramesPerAnalogFrame;
	for (int channel = 0; channel < analogInputs; channel++)
	{
		analogWriteOnce(context, frame, channel, m_lfo[channel]);
	}
}
