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

	// Order of params in DynamicDSS: L, rho, r, T, E, sigma0, sigma1
	// Do not change this
	std::vector<std::string> originalParameterOrder = { "L", "rho", "r", "T", "E", "sigma0", "sigma1" };
	for (int i = 0; i < originalParameterOrder.size(); i++)
	{
		const std::string& parameterName = originalParameterOrder[i];
		m_parameterIdMap[parameterName] = i;
	}

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
	std::map<std::string, std::pair<float, float>> parameterRanges;
	parameterRanges["L"] = { 0.2f, 4.0f };
	parameterRanges["rho"] = { 1962.5f, 15700.0f };
	parameterRanges["r"] = { 0.0005f, 0.001f };
	parameterRanges["T"] = { 150.f, 1200.0f };
	parameterRanges["E"] = { 0.f, 400000000000.f };
	parameterRanges["sigma0"] = { 0.000001f, 2.f };
	parameterRanges["sigma1"] = { 0.0002f, 0.01f };
	parameterRanges["loc"] = { 0.f, 1.f };

	// Order of inputs is L, rho, T, r, loc, E, sigma0, sigma1
	// Change this if you want to reorder inputs on the device
	const std::vector<std::string> parameterOrder = { "L", "rho", "T", "r", "loc", "E", "sigma0", "sigma1" };
	m_analogInputs.reserve(sAnalogInputCount);
	for (int i = 0; i < 8; i++)
	{
		const std::string& parameterName = parameterOrder[i];
		m_analogInputs.push_back(AnalogInput(parameterName, i, parameterRanges[parameterName]));
		m_labelToAnalogIn[parameterName] = i;
		if (parameterName != "loc") m_channelsToUpdate.insert(i); // Force update to read initial values
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
			const float mappedValue = m_analogInputs[channel].getCurrentValueMapped();
			rt_printf("Updating channel %d with value %f\n", channel, mappedValue);

			// Map the analog channel to intended parameter to parameter id in DSS simulation
			const auto& analogIn = m_analogInputs[channel];
			const int parameterId = m_parameterIdMap[analogIn.getLabel()];
			
			if(clippingFlag)
			{
				// sigma0
				int parameterId = m_parameterIdMap["sigma0"];
				float currentValue = m_analogInputs[m_labelToAnalogIn["sigma0"]].getCurrentValueMapped();
				m_pDynamicStiffString->refreshParameter(parameterId, currentValue * 10.f);
				// sigma1
				parameterId = m_parameterIdMap["sigma1"];
				currentValue = m_analogInputs[m_labelToAnalogIn["sigma1"]].getCurrentValueMapped();
				m_pDynamicStiffString->refreshParameter(parameterId, currentValue * 10.f);
			
				clippingFlag = false;
			}
			
			m_pDynamicStiffString->refreshParameter(parameterId, mappedValue);
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
			// We will always update channel 7 (excitation loc) as this param is only effective during excitation
			// So there is no risk of too frequent updates
			if (analogIn.getLabel() == "loc")
			{
				m_excitationLoc = analogIn.getCurrentValueMapped();
			}
			else if (analogIn.hasChanged())
			{
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

	float output = Global::limit(m_pDynamicStiffString->getOutput(), -5.f, 5.f);
	/* Audio range increased then mapped back
	 5x times more headroom (no clipping unless under loooots of extress) 
	 but will need external compression - will use compressor modules.
	 This would be the right place to look into compression algorithms.
	 Practically its own reasearch branch - > We can come back to this next month. */
	
	if ((output >= 3.f) || (output <= -3.f))
	{
		clippingFlag = true;
	}
	
	output = map(output, -5.f, 5.f, -1.f, 1.f);
	
	for (unsigned int channel = 0; channel < context->audioOutChannels; channel++)
	{
		audioWrite(context, frame, channel, output);
	}
	
}
