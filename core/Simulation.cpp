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

	for (auto& parameter : m_parameters.getParameters())
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
	if (m_buttons[Button::Type::TRIGGER].isPressed())
	{
		sprayValue = Global::f_random(-1, 1) * sprayAmount;
		sprayedloc = nonsprayloc + sprayValue;
		sprayedloc = (sprayedloc < 0) ? -sprayedloc : (sprayedloc > 1) ? (1 - (sprayedloc - 1)) : sprayedloc;

		auto fnExcitation = m_fnRaisedCos;

		m_parameters.getParameter(ParameterName::loc).setValue(sprayedloc);

		m_screen.setBrightness(m_parameters.getParameter(ParameterName::loc).getChannel(), sprayedloc);
		m_pDynamicStiffString->excite(m_parameters.getParameter(ParameterName::loc).getValue(), fnExcitation);
	}

	// 2. Handle Audio Excitation (Audio Excitation should be triggered unrelated to Trigger Button)
	if (!m_audioBuffer.containsSilence())
	{
		auto fnExcitation = m_fnSampleExcitation;
		m_pDynamicStiffString->excite(m_parameters.getParameter(ParameterName::loc).getValue(), fnExcitation);
	}



	// 2. Process parameter changes
	if (m_updateFrameCounter == 0)
	{
		// Process update queue
		for (const auto parameterName : m_parametersToUpdate)
		{
			auto& parameter = m_parameters.getParameter(parameterName);

			float mappedValue;

			if (parameterName == ParameterName::L) // 1 Volt per Octave
			{
				float lValue_inVolts = parameter.getAnalogInput()->getCurrentValueinVolts();

				while (lValue_inVolts > 3.f)		   // Length from 2m to 0.5m - 3 Octaves
				{								       // While loop converts any number over 3
					lValue_inVolts -= 3.f;			   // Back to the 0-3 range (but not including our first 3)
				}

				mappedValue = 2.f * powf(2, -lValue_inVolts);
			}

			else
				mappedValue = parameter.getAnalogInput()->getCurrentValueMapped();


			// Save state and send change to DSS
			parameter.setValue(mappedValue);
			m_pDynamicStiffString->refreshParameter(parameter.getId(), mappedValue);

			// Update screen
			//rt_printf("sending channel %d to the LEDScreen with value %f\n", parameter.getChannel(), mappedValue);
			m_screen.setBrightness(parameter.getChannel(), parameter.getAnalogInput()->unmapValue(mappedValue)); // Passes the parameter's value to the correct channel
		}

		if (hasCorrectedFlag) //Only After Correction and Once Signal has stabilized
		{
			// compares pots with increased updatedValues, then decreases values at speed depending on distance

			auto& sigma1 = m_parameters.getParameter(ParameterName::sigma1);
			float sigma1pot = sigma1.getAnalogInput()->getCurrentValue(); // min value is 0.0008f
			float sigma1cor = sigma1.getValue(); // max     value is 1
			float sigma1dif = sigma1pot - sigma1cor; // max value is 0.9992
			float sigma1coef = sigma1dif * 0.0002;

			sigma1.setValue(sigma1cor + sigma1coef);
			m_pDynamicStiffString->refreshParameter(sigma1.getId(), sigma1cor - sigma1coef);

			if (sigma1dif == 0)
				hasCorrectedFlag = false;
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
	for (auto& parameter : m_parameters.getParameters())
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
	m_audioBuffer.put(map(audioRead(context, frame, 0), -1.f, 1.f, -0.5f, 0.5f));

	if (!(frame % m_audioFramesPerAnalogFrame))
	{
		const int analogFrame = frame / m_audioFramesPerAnalogFrame;
		for (auto& kvp : m_parameters.getParameters()) // Start with channels 1-8; 0 is reserved for trigger
		{
			auto& parameter = kvp.second;
			auto analogIn = kvp.second.getAnalogInput();

			analogIn->read(context, analogFrame);

			// We will always update L (length). as this param needs to be as accurate as possible
			// since it affects pitch.

			if (parameter.getName() == ParameterName::L)
				m_parametersToUpdate.insert(parameter.getName());

			// We will always update channel 7 (excitation loc) as this param is only effective during excitation
			// So there is no risk of too frequent updates
			if (parameter.getName() == ParameterName::loc)
			{

				// Handle Spray button
				if (m_buttons[Button::Type::SPRAY].isHeld())
				{
					sprayAmount = analogIn->getCurrentValueMapped();
					m_screen.setBrightness(m_parameters.getParameter(ParameterName::loc).getChannel(), sprayAmount);
					//rt_printf("Spray Amount is %f\n", sprayAmount);
				}

				nonsprayloc = analogIn->getCurrentValueMapped();

				if (m_buttons[Button::Type::SPRAY].isReleased() || analogIn->hasChanged())
					m_screen.setBrightness(parameter.getChannel(), nonsprayloc);

				//rt_printf("nonsprayLoc is %f\n", nonsprayloc);
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

	double output = m_pDynamicStiffString->getOutput();

	if (output >= 10000000.f)
	{
		auto& sigma1 = m_parameters.getParameter(ParameterName::sigma1);
		sigma1.setValue(1.f);
		m_pDynamicStiffString->refreshParameter(sigma1.getId(), 1.f);
		m_screen.setBrightness(sigma1.getChannel(), 1.f);
		hasCorrectedFlag = true;
	}

	Compressor.process(output, output);

	for (unsigned int channel = 0; channel < context->audioOutChannels; channel++)
	{
		audioWrite(context, frame, channel, output);
	}
}
