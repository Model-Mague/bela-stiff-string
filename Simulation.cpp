#include "Simulation.h"

#define _USE_MATH_DEFINES
#include <cmath>
#include <tuple>
#include <cstring>
#include <sstream>

#define INTERACTION_DELAY 1000 // frames


Simulation::Simulation(BelaContext* context) : m_amplitude(5.f), m_frequency(0.1f), m_updateParameters(true), m_excitationLoc(-1.f)
{
	m_inverseSampleRate = 1.0f / context->audioSampleRate;
	if (context->analogFrames)
		m_audioFramesPerAnalogFrame = context->audioFrames / context->analogFrames;

	// Setup Dynamic Stiff String
	DynamicStiffString::SimulationParameters parameters = {};
	parameters.L = 1.0f;
	parameters.rho = 7850.0f;
	parameters.r = 0.0005f;
	parameters.T = 300.0f;
	parameters.E = 200000000000.f;
	parameters.sigma0 = 1.0f;
	parameters.sigma1 = 0.005f;
	m_pDynamicStiffString = std::make_unique<DynamicStiffString>(parameters, m_inverseSampleRate);

	// Setup parameter ranges
	std::vector<std::pair<float, float>> parameterRanges;
	parameterRanges.reserve(8);
	parameterRanges.push_back({ 0.5f, 2.0f }); // L
	parameterRanges.push_back({ 3925.0f, 15700.0f }); // rho
	parameterRanges.push_back({ 0.00025f, 0.001f }); // r
	parameterRanges.push_back({ 150.f, 600.0f }); // T
	parameterRanges.push_back({ 0.f, 400000000000.f }); // E
	parameterRanges.push_back({ 0.f, 2.f }); // sigma0
	parameterRanges.push_back({ 0.0002f, 0.01f }); // sigma1
	parameterRanges.push_back({ 0.f, 1.f }); // loc (string excite position)

	// Setup analog inputs
	m_analogInputs.reserve(sAnalogInputCount);
	for (int i = 0; i < 8; i++)
	{
		m_analogInputs.push_back(AnalogInput(i, parameterRanges[i]));
	}

	m_amplitude = 5;
	m_frequency = 0.1f;
}

void Simulation::update(BelaContext* context)
{
	// 1. Handle parameter changes
	if (m_updateParameters)
	{
		//rt_printf("Updating parameters\n");
		m_pDynamicStiffString->refreshCoefficients();
		m_pDynamicStiffString->calculateScheme();
		m_pDynamicStiffString->updateStates();
		m_updateParameters = false;
	}

	// 2. Handle trigger button (should probably be last)
	if (isButtonReleased(Button::TRIGGER))
	{
		//rt_printf("Trigger button released\n");
		m_pDynamicStiffString->excite(m_excitationLoc);
	}

	// 3. Update phase
	for (int channel = 0; channel < sAnalogInputCount; channel++)
	{
		m_phase[channel] += 2.0f * (float)M_PI * m_frequency * m_inverseSampleRate;
		if (m_phase[channel] > M_PI)
			m_phase[channel] -= 2.0f * (float)M_PI;
	}
}

std::string Simulation::getCalibrationResults()
{
	std::stringstream ss;
	for (int i = 0; i < m_analogInputs.size(); i++)
	{
		float minValue, maxValue;
		std::tie(minValue, maxValue) = m_analogInputs[i].getValueRange();
		ss << "Pot " << i << " [" << minValue << ", " << maxValue << "]" << std::endl;
	}
	return ss.str();
}

void Simulation::readInputs(BelaContext* context, int frame)
{
	frame = frame / m_audioFramesPerAnalogFrame;
	for (int channel = 0; channel < sAnalogInputCount; channel++) // Start with channels 1-8; 0 is reserved for trigger
	{
		auto& analogIn = m_analogInputs[channel];
		m_analogIn[channel] = analogIn.read(context, frame);
		if (analogIn.hasChanged())
		{
			//rt_printf("New analog input at channel %d: %f\n", channel, m_analogIn[channel]);
			if (channel == 7)
			{
				m_excitationLoc = m_analogIn[channel];
			} else {
				m_pDynamicStiffString->refreshParameter(channel, m_analogIn[channel]);
				m_updateParameters = true;
			}

		}
	}

	// Buttons, left to right
	// 15, 14, 13, 12
	// Trigger, Mode, Up, Down
	static int buttonChannel[4] = { 15, 14, 13, 12 };
	for (const auto& button : { Button::TRIGGER, Button::MODE, Button::UP, Button::DOWN })
	{
		m_buttonPreviousState[(size_t)button] = m_buttonState[(size_t)button];
		m_buttonState[(size_t)button] = digitalRead(context, frame, buttonChannel[(size_t)button]);
	}
}

void Simulation::writeOutputs(BelaContext* context, int frame)
{
	frame = frame / m_audioFramesPerAnalogFrame;
	for (int channel = 0; channel < sAnalogInputCount; channel++)
	{
		analogWriteOnce(context, frame, channel, m_lfo[channel]);
	}
}

void Simulation::writeAudio(BelaContext* context, int frame)
{
	float output = Global::limit(m_pDynamicStiffString->getOutput(), -1.0, 1.0);
	for (unsigned int channel = 0; channel < context->audioOutChannels; channel++)
	{
		audioWrite(context, frame, channel, output);
	}
}
