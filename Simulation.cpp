#include "Simulation.h"

#define _USE_MATH_DEFINES
#include <cmath>

#define INTERACTION_DELAY 1000 // frames

Simulation::Simulation(BelaContext* context) : m_amplitude(5.f), m_frequency(0.1f), m_buttonPressed()
{
	m_inverseSampleRate = 1.0f / context->audioSampleRate;
	if (context->analogFrames)
		m_audioFramesPerAnalogFrame = context->audioFrames / context->analogFrames;

	memset((void*)&m_buttonPressed, 0, 4);
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
	if (isTriggerPressed())
	{
		const auto lastActivated = m_buttonLastActivated[(size_t)Button::TRIGGER];
		if (((context->audioFramesElapsed - lastActivated) > INTERACTION_DELAY) || (lastActivated == 0))
		{
			m_buttonLastActivated[(size_t)Button::TRIGGER] = context->audioFramesElapsed;

			// Trigger the stiff string
			pDynamicStiffString->excite();
		}
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
	m_buttonPressed[(size_t)Simulation::Button::TRIGGER] = (bool)digitalRead(context, nFrame, 15);
	m_buttonPressed[(size_t)Simulation::Button::MODE] = (bool)digitalRead(context, nFrame, 14);
	m_buttonPressed[(size_t)Simulation::Button::UP] = (bool)digitalRead(context, nFrame, 13);
	m_buttonPressed[(size_t)Simulation::Button::DOWN] = (bool)digitalRead(context, nFrame, 12);
}

void Simulation::writeOutputs(BelaContext* context, int nFrame)
{
	const int frame = nFrame / m_audioFramesPerAnalogFrame;
	for (int channel = 0; channel < analogInputs; channel++)
	{
		analogWriteOnce(context, frame, channel, m_lfo[channel]);
	}
}
