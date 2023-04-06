#include "Simulation.h"

#define _USE_MATH_DEFINES
#include <cmath>
#include <tuple>
#include <cstring>
#include <sstream>


Simulation::Simulation(BelaContext* context) : m_excitationLoc(-1.f), m_amplitude(5.f), m_frequency(0.1f)
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
	parameterRanges.push_back({ 0.2f, 4.0f }); // L
	parameterRanges.push_back({ 1962.5f, 15700.0f }); // rho
	parameterRanges.push_back({ 0.0005f, 0.001f }); // r
	parameterRanges.push_back({ 150.f, 1200.0f }); // T
	parameterRanges.push_back({ 0.f, 400000000000.f }); // E
	parameterRanges.push_back({ 0.f, 2.f }); // sigma0
	parameterRanges.push_back({ 0.000f, 0.01f }); // sigma1
	parameterRanges.push_back({ 0.f, 1.f }); // loc (string excite position)

	// Setup analog inputs
	m_analogInputs.reserve(sAnalogInputCount);
	for (int i = 0; i < 8; i++)
	{
		m_analogInputs.push_back(AnalogInput(i, parameterRanges[i]));
		if (i != 7) m_channelsToUpdate.insert(i); // Force update to read initial values
	}

	m_amplitude = 5;
	m_frequency = 0.1f;
}

void Simulation::update(BelaContext* context)
{
	// 1. Handle trigger button (should probably be last)
	if (isButtonReleased(Button::TRIGGER))
	{
		m_pDynamicStiffString->excite(m_excitationLoc);
	}

	// 2. Process parameter changes
	if (!m_channelsToUpdate.empty() && (m_updateFrameCounter == 0))
	{
		// Process update queue
		for (int channel: m_channelsToUpdate)
		{
			rt_printf("Updating channel %d with value %f\n", channel, m_rangeMappedInputs[channel]);
			m_pDynamicStiffString->refreshParameter(channel, m_rangeMappedInputs[channel]);
		}
		m_channelsToUpdate.clear();
		m_updateFrameCounter = sDSSUpdateRate;
	}
	else
	{
		m_updateFrameCounter = std::max(0, m_updateFrameCounter - 1);
	}

	// 3. Update DSS simulation
	m_pDynamicStiffString->refreshCoefficients();
	m_pDynamicStiffString->calculateScheme();
	m_pDynamicStiffString->updateStates();

	// 4. Update LFO phase
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
	if (!(frame % m_audioFramesPerAnalogFrame))
	{
		const int analogFrame = frame / m_audioFramesPerAnalogFrame;
		for (int channel = 0; channel < sAnalogInputCount; channel++) // Start with channels 1-8; 0 is reserved for trigger
		{
			auto& analogIn = m_analogInputs[channel];
			analogIn.read(context, analogFrame);
			m_rangeMappedInputs[channel] = analogIn.getCurrentValueMapped();
			// We will always update channel 7 (excitation loc) as this param is only effective during excitation
			// So there is no risk of too frequent updates
			if (channel == 7)
			{
				m_excitationLoc = m_rangeMappedInputs[channel];
			}
			else if (analogIn.hasChanged())
			{
				//rt_printf("New analog input at channel %d: %f\n", channel, m_rangeMappedInputs[channel]);
				// Register channel as one that needs its value read and updated
				m_channelsToUpdate.insert(channel);
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
		
	/* Potential De-Clipping solution? How do I reference sigma? :)
	
	if (output >= 0.95f)
	{
		sigma0 = sigma0*10;
		sigma1 = sigma1*10;
	}
	
	if (output <= -0.95f)
	{
		sigma0 = sigma0*10;
		sigma1 = sigma1*10;
	}
	
	*/
	
	
}
