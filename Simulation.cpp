#include "Simulation.h"

#include "LEDScreen.h"

#define _USE_MATH_DEFINES
#include <cmath>
#include <tuple>
#include <cstring>
#include <sstream>


Simulation::Simulation(BelaContext* context) : m_amplitude(5.f), m_frequency(0.1f), m_screen(context)
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

	// Order of inputs is L, rho, T, r, loc, E, sigma0, sigma1
	// Change this if you want to reorder inputs on the device
	m_analogInputs.reserve(Parameters::sNames.size());
	for (int i = 0; i < Parameters::sNames.size(); i++)
	{
		const std::string& parameterName = Parameters::sNames[i];
		m_analogInputs.push_back(AnalogInput(parameterName, i, m_parameters.getParameter(parameterName).getRange()));
		m_labelToAnalogIn[parameterName] = i;
		if (parameterName != "loc") m_channelsToUpdate.insert(i); // Force update to read initial values
	}
}


void Simulation::update(BelaContext* context)
{
	// 1. Handle trigger button (should probably be last)
	if (m_buttons[Button::Type::TRIGGER].isReleased())
	{
		sprayValue = ((rand() % 100)/100.f); // Since rand() only delivers integers, the range 0 - 100 is divided to be 0.00 - 1.00
		
		if (sprayValue > 50) // Since rand() only delivers positive numbers, any number > 50 is - 100 -> making the range -0.50 -> 0.50
		{
			sprayValue = sprayValue - 100;
		}
		sprayValue = sprayValue * sprayAmount;
		
				
		if (sprayedloc > 1.f || sprayedloc < 0.f) // if over 1 or below 0, sprayValue reversed.
		{
			sprayValue =  - sprayValue * 2;
			sprayedloc = sprayedloc + sprayValue;
		}

		rt_printf("sprayValue is %f\n", sprayValue);
		
		m_pDynamicStiffString->excite(m_parameters.getParameter("loc").getValue());
	}

	// 2. Process parameter changes
	if (m_updateFrameCounter == 0)
	{
		// Process update queue
		for (int channel: m_channelsToUpdate)
		{
			const auto& analogIn = m_analogInputs[channel];
			float mappedValue = m_analogInputs[channel].getCurrentValueMapped();
			if (channel == 0) // Length Handled differently -> 1V/oct 
			{
				mappedValue = 0.5f * powf(2, map(Global::limit(mappedValue, 0.f, 1.33f), 0.5f, 1.33f, 3.f, 0.f)); // <- magic: input limited to 0-3V (which are the number of octaves made available by changing the Length in our range, then mapped to oposite values, then made exponent. It is very messy, sort of a desperate measure tbh.
			}

			rt_printf("Updating channel %d with value %f\n", channel, mappedValue);
	
			// Map the analog channel to intended parameter to parameter id in DSS simulation
			const std::string& paramName = analogIn.getLabel();
			auto& parameter = m_parameters.getParameter(paramName);

			// Save state and send change to DSS
			parameter.setValue(mappedValue);
			m_pDynamicStiffString->refreshParameter(parameter.getId(), mappedValue);
			
			// Update screen
			rt_printf("sending channel %d to the LEDScreen with value %f\n", channel, mappedValue);
			m_screen.setBrightness(channel, analogIn.unmapValue(mappedValue)); // Passes the parameter's value to the correct channel
		}
		
		if (clippingFlag == true)
		{
			// sigma0
			auto& sigma0 = m_parameters.getParameter("sigma0");
			float updatedValue = std::min(sigma0.getValue() * correctionValue, 2.f);
			sigma0.setValue(updatedValue);
			m_pDynamicStiffString->refreshParameter(sigma0.getId(), updatedValue);
			rt_printf("Updating sigma0 with value %f\n", updatedValue);

			// Update screen
			const auto& analogInSigma0 = m_analogInputs[m_labelToAnalogIn["sigma0"]];
			m_screen.setBrightness(7, analogInSigma0.unmapValue(updatedValue)); // passes the new value to the LEDScreen

			// sigma1
			auto& sigma1 = m_parameters.getParameter("sigma1");
			updatedValue = std::min(sigma1.getValue() * correctionValue, 0.01f);
			sigma1.setValue(updatedValue);
			m_pDynamicStiffString->refreshParameter(sigma1.getId(), updatedValue);
			rt_printf("Updating sigma1 with value %f\n", updatedValue);

			// Update screen
			const auto& analogInSigma1 = m_analogInputs[m_labelToAnalogIn["sigma1"]];
			m_screen.setBrightness(8, analogInSigma1.unmapValue(updatedValue)); // passes the new value to the LEDScreen
			
			clippingFlag = false;	
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
				// Handle Spray button
				if (m_buttons[Button::Type::SPRAY].isPressed())
				{
					sprayAmount = analogIn.getCurrentValueMapped();
					rt_printf("Spray Amount is %f\n", sprayAmount);
				}
				
				sprayedloc = analogIn.getCurrentValueMapped() + sprayValue;
				
				m_parameters.getParameter("loc").setValue(sprayedloc);
				m_screen.setBrightness(channel, analogIn.getCurrentValue()); // needs mapped
			}
			else if (analogIn.hasChanged())
			{
				// Register channel as one that needs its value read and updated
				m_channelsToUpdate.insert(channel);
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
	/* Audio range increased then mapped back
	 5x times more headroom (no clipping unless under loooots of extress) 
	 but will need external compression - will use compressor modules.
	 This would be the right place to look into compression algorithms.
	 Practically its own reasearch branch - > We can come back to this next month. */
	
	if ((output >= 2.5f) || (output <= - 2.5f))
	{
		clippingFlag = true;
		correctionValue = 1.001f;
		
		if ((output >= 3.5f) || (output <= - 3.5f))
		{
			correctionValue = 1.01f;
		}
		if ((output >= 4.f) || (output <= - 4.5f))
		{
			correctionValue = 1.1f;
		}
		if ((output >= 4.5f) || (output <= - 4.5f))
		{
			correctionValue = 5.f;
		}
		if ((output >= 4.9f) || (output <= - 4.9f))
		{
			correctionValue = 10.f;
		}
	}
	
	output = map(output, -5.f, 5.f, -1.f, 1.f);
	
	for (unsigned int channel = 0; channel < context->audioOutChannels; channel++)
	{
		audioWrite(context, frame, channel, output);
	}
}
