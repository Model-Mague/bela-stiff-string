#include "Simulation.h"

Simulation::Simulation(BelaContext* context) : m_amplitude(5.f), m_frequency(0.1f)
{
	m_inverseSampleRate = 1.0f / context->audioSampleRate;
	if (context->analogFrames)
		m_audioFramesPerAnalogFrame = context->audioFrames / context->analogFrames;

	memset((void*)&m_analogIn, 0, analogInputs * sizeof(float));
	memset((void*)&m_lfo, 0, analogInputs * sizeof(float));
	memset((void*)&m_phase, 0, analogInputs * sizeof(float));

	m_amplitude = 5;
	m_frequency = 0.1f;
}

void Simulation::readInputs(BelaContext* context, int nFrame)
{
	const int frame = nFrame / m_audioFramesPerAnalogFrame;
	for (int channel = 0; channel < analogInputs; channel++)
	{
		m_analogIn[channel] = analogRead(context, frame, channel);
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
