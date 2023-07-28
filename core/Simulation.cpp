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

		//rt_printf("sprayValue is %f\n", sprayValue);
		//rt_printf("sprayedloc is %f\n", sprayedloc);

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

				while (lValue_inVolts > 4.f)		   // Four Octaves Available (Length 2m to 0.125m)
				{								       // While loop converts any number over 4
					lValue_inVolts -= 4.f;			   // Back to the 0-4 range
				}

				mappedValue = 2.f * powf(2, -lValue_inVolts);
			}

			else
				mappedValue = parameter.getAnalogInput()->getCurrentValueMapped();


			//rt_printf("Updating channel %d with value %f\n", parameter.getChannel(), mappedValue);

			// Save state and send change to DSS
			parameter.setValue(mappedValue);
			m_pDynamicStiffString->refreshParameter(parameter.getId(), mappedValue);

			// Update screen
			//rt_printf("sending channel %d to the LEDScreen with value %f\n", parameter.getChannel(), mappedValue);
			m_screen.setBrightness(parameter.getChannel(), parameter.getAnalogInput()->unmapValue(mappedValue)); // Passes the parameter's value to the correct channel
		}

		/*if (clippingFlag == true)
		{
			// sigma0
			/* auto& sigma0 = m_parameters.getParameter(ParameterName::sigma0);
			float updatedValue = std::min(sigma0.getValue() * correctionValue, 2.f);
			sigma0.setValue(updatedValue);
			m_pDynamicStiffString->refreshParameter(sigma0.getId(), updatedValue);
			//rt_printf("Updating sigma0 with value %f\n", updatedValue);

			// Update screen
			m_screen.setBrightness(sigma0.getChannel(), updatedValue); // passes the new value to the LEDScreen


			// sigma1
			auto& sigma1 = m_parameters.getParameter(ParameterName::sigma1);
			float updatedValue = std::min(sigma1.getValue() * correctionValue, 1.f);
			sigma1.setValue(updatedValue);
			m_pDynamicStiffString->refreshParameter(sigma1.getId(), updatedValue);
			//rt_printf("Updating sigma1 with value %f\n", updatedValue);

			// Update screen
			m_screen.setBrightness(sigma1.getChannel(), updatedValue); // passes the new value to the LEDScreen

			clippingFlag = false;
			hasCorrectedFlag = true;
		}

		if (hasCorrectedFlag && stableFlag) //Only After Correction and Once Signal has stabilized
		{
			// compares pots with increased updatedValues, then decreases values at speed depending on distance

			/*auto& sigma0 = m_parameters.getParameter(ParameterName::sigma0);
			float sigma0pot = sigma0.getAnalogInput()->getCurrentValue(); // min value is 0
			float sigma0cor = sigma0.getValue(); // max value is 2
			float sigma0dif = sigma0pot - sigma0cor; // max value is 2
			float sigma0coef = sigma0dif * 0.0002;


			sigma0.setValue(sigma0cor + sigma0coef);
			m_pDynamicStiffString->refreshParameter(sigma0.getId(), sigma0cor - sigma0coef);


			auto& sigma1 = m_parameters.getParameter(ParameterName::sigma1);
			float sigma1pot = sigma1.getAnalogInput()->getCurrentValue(); // min value is 0.0008f
			float sigma1cor = sigma1.getValue(); // max     value is 1
			float sigma1dif = sigma1pot - sigma1cor; // max value is 0.9992
			float sigma1coef = sigma1dif * 0.0002;

			sigma1.setValue(sigma1cor + sigma1coef);
			m_pDynamicStiffString->refreshParameter(sigma1.getId(), sigma1cor - sigma1coef);

			if(/*sigma0dif == 0 && sigma1dif == 0)
			hasCorrectedFlag = false;
		}
		*/
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
	float l_Range = 10.f; // loudness range
	double output = Global::limit(m_pDynamicStiffString->getOutput(), -l_Range, l_Range);

	before_compression.push_back(std::abs(output));

	Compressor.process(output, output);

	after_compression.push_back(std::abs(output));

	/*if ((output >= l_Range) || (output <= -l_Range))
	{
		float positive_output = output >= 0 ? output : -output;

		clippingFlag = true;
		stableFlag = false;
		correctionValue = powf(1.1, positive_output - (l_Range - 0.8 * l_Range));
	}

	else
		stableFlag = true;
	*/
	
	output = map(output, -l_Range, l_Range, -1.f, 1.f);

	for (unsigned int channel = 0; channel < context->audioOutChannels; channel++)
	{
		audioWrite(context, frame, channel, output);
	}

	if (before_compression.size() == context->audioFrames)
	{
		result_before = std::max_element(before_compression.begin(), before_compression.end());
		rt_printf("b: \f", *result_before);

		result_after = std::max_element(after_compression.begin(), after_compression.end());
		rt_printf("a: \f", *result_after);

		before_compression.clear();
		after_compression.clear();	
	}
}
