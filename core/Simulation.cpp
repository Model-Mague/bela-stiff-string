#include "Simulation.h"

#define _USE_MATH_DEFINES
#include <algorithm>
#include <cmath>
#include <tuple>
#include <cstring>
#include <sstream>
#include <random>


Simulation::Simulation(BelaContext* context) : m_amplitude(5.f), m_frequency(0.1f), m_screen(context), m_audioBuffer(10)
{
	m_inverseSampleRate = 1.0f / context->audioSampleRate;
	if (context->analogFrames)
		m_audioFramesPerAnalogFrame = context->audioFrames / context->analogFrames;

	// Setup buttons
	static const int buttonChannel[4] = { 15, 14, 13, 12 };
	static const Button::Type buttonTypes[4] = { Button::Type::TRIGGER, Button::Type::MODE, Button::Type::SPRAY, Button::Type::DOWN };
	for (int i = 0; i < 4; i++)
	{
		m_buttons[buttonTypes[i]] = Button(buttonChannel[i]);
	}

	m_amplitude = 5;
	m_frequency = 0.1f;

	m_pDynamicStiffString = std::make_unique<DynamicStiffString>(m_parameters.getDSSParameters(), m_inverseSampleRate);

	for (auto& parameter: m_parameters.getParameters())
	{
		const auto name = parameter.second.getName();
		if (name != ParameterName::loc) m_parametersToUpdate.insert(name); // Force update to read initial values
	}

	// Setup excitation functions
	m_fnRaisedCos = [](int index, int size) {
		return 0.5f * (1 - cos(2.0f * static_cast<float>(M_PI) * index / (size - 1.0f)));
	};
	m_fnSampleExcitation = [this](int index, int size) {
		return m_audioBuffer.get(index);
	};
}


void Simulation::update(BelaContext* context)
{
	// 1. Handle trigger button (should probably be last)
	if (m_buttons[Button::Type::TRIGGER].isReleased())
	{
		static std::mt19937 gen(0); 
		static std::uniform_real_distribution<float> dis(-0.5, 0.5);

		sprayValue = dis(gen) * sprayAmount;
		sprayedloc += sprayValue;

		sprayedloc = [](float location) 
		{ 
			if (location < 0) return -location;
			else if (location > 1) return 1 - (location - 1);
			else return location; 
		} (sprayedloc);
		
		rt_printf("sprayValue is %f\n", sprayValue);

		auto fnExcitation = m_audioBuffer.containsSilence() ? m_fnRaisedCos : m_fnSampleExcitation;
		m_pDynamicStiffString->excite(m_parameters.getParameter(ParameterName::loc).getValue(), fnExcitation);
	}

	// 2. Process parameter changes
	if (m_updateFrameCounter == 0)
	{
		// Process update queue
		for (const auto parameterName: m_parametersToUpdate)
		{
			auto& parameter = m_parameters.getParameter(parameterName);

			float mappedValue = parameter.getAnalogInput()->getCurrentValueMapped();

			if (parameterName == ParameterName::L) // 1 Volt per Octave
			{
				while (mappedValue > 3.f)		   // Three Octaves Available (Length 4m to 0.5m)
				{								   // While loop converts any number over 3
					mappedValue -= 3.f;			   // Back to the 0-3 range
				}
				mappedValue = 4.f * powf(2, -mappedValue);
			}

			rt_printf("Updating channel %d with value %f\n", parameter.getChannel(), mappedValue);
	
			// Save state and send change to DSS
			parameter.setValue(mappedValue);
			m_pDynamicStiffString->refreshParameter(parameter.getId(), mappedValue);
			
			// Update screen
			rt_printf("sending channel %d to the LEDScreen with value %f\n", parameter.getChannel(), mappedValue);
			m_screen.setBrightness(parameter.getChannel(), parameter.getAnalogInput()->unmapValue(mappedValue)); // Passes the parameter's value to the correct channel
		}
		
		if (clippingFlag == true)
		{
			// sigma0
			auto& sigma0 = m_parameters.getParameter(ParameterName::sigma0);
			float updatedValue = std::min(sigma0.getValue() * correctionValue, 2.f);
			sigma0.setValue(updatedValue);
			m_pDynamicStiffString->refreshParameter(sigma0.getId(), updatedValue);
			rt_printf("Updating sigma0 with value %f\n", updatedValue);

			// Update screen
			m_screen.setBrightness(7, sigma0.getAnalogInput()->unmapValue(updatedValue)); // passes the new value to the LEDScreen

			// sigma1
			auto& sigma1 = m_parameters.getParameter(ParameterName::sigma1);
			updatedValue = std::min(sigma1.getValue() * correctionValue, 0.01f);
			sigma1.setValue(updatedValue);
			m_pDynamicStiffString->refreshParameter(sigma1.getId(), updatedValue);
			rt_printf("Updating sigma1 with value %f\n", updatedValue);

			// Update screen
			m_screen.setBrightness(8, sigma1.getAnalogInput()->unmapValue(updatedValue)); // passes the new value to the LEDScreen
			
			clippingFlag = false;	
		}
		
		m_parametersToUpdate.clear();
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

	// 4. Update screen
	m_screen.update(context);

	// 5. Update LFO phase
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
	for (auto& parameter: m_parameters.getParameters())
	{
		float minValue, maxValue;
		std::tie(minValue, maxValue) = parameter.second.getRange();
		ss << "Pot " << parameter.second.getChannel() << " [" << minValue << ", " << maxValue << "]" << std::endl;
	}
	return ss.str();
}

void Simulation::readInputs(BelaContext* context, int frame)
{
	// First, read in a frame of audio from the input
	m_audioBuffer.put(audioRead(context, frame, 0));

	if (!(frame % m_audioFramesPerAnalogFrame))
	{
		const int analogFrame = frame / m_audioFramesPerAnalogFrame;
		for (auto kvp: m_parameters.getParameters()) // Start with channels 1-8; 0 is reserved for trigger
		{
			auto& parameter = kvp.second;
			auto analogIn = kvp.second.getAnalogInput();

			analogIn->read(context, analogFrame);
			// We will always update channel 7 (excitation loc) as this param is only effective during excitation
			// So there is no risk of too frequent updates
			if (parameter.getName() == ParameterName::loc)
			{
				// Handle Spray button
				if (m_buttons[Button::Type::SPRAY].isPressed())
				{
					sprayAmount = analogIn->getCurrentValueMapped();
					rt_printf("Spray Amount is %f\n", sprayAmount);
				}
				
				sprayedloc = analogIn->getCurrentValueMapped() + sprayValue;
				
				m_parameters.getParameter(ParameterName::loc).setValue(sprayedloc);
				m_screen.setBrightness(parameter.getChannel(), analogIn->getCurrentValue()); // needs mapped
			}
			else if (analogIn->hasChanged())
			{
				// Register channel as one that needs its value read and updated
				m_parametersToUpdate.insert(parameter.getName());
			}
		}
	}

	for (auto& kvp : m_buttons)
	{
		const int read = digitalRead(context, frame, kvp.second.getChannel());
		kvp.second.updateState(read);
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
	
	if ((output >= 2.5f) || (output <= - 2.5f))
	{
		clippingFlag = true;		
		correctionValue = 1 + (powf(output, 10.f) / powf(10, 6.f)); // Very steep exponential function
	}
	
	output = map(output, -5.f, 5.f, -1.f, 1.f);
	
	for (unsigned int channel = 0; channel < context->audioOutChannels; channel++)
	{
		audioWrite(context, frame, channel, output);
	}
}
