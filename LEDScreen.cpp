#include "LEDScreen.h"

LEDScreen::LEDScreen(BelaContext* context) : m_phase(0.f)
{
	for (unsigned int i = 0; i < sCellCount; i++)
	{
		pinMode(context, 0, sLedPins[i], OUTPUT);
	}

	m_invSampleRate = 1 / context->digitalSampleRate;
}

void LEDScreen::setBrightness(const int cellIndex, const float brightness)
{
	m_brightness[cellIndex] = brightness;
}

void LEDScreen::update(BelaContext* context)
{
	for (unsigned int n = 0; n < context->digitalFrames; ++n)
	{
		computePhase();

		for (int i = 0; i < sCellCount; i++)
		{
			digitalWriteOnce(context, n, sLedPins[i], squarePWM(m_brightness[i]));
		}
	}
}

void LEDScreen::computePhase()
{
	m_phase += 440 * m_invSampleRate;
	if (m_phase >= 1)
		m_phase -= 1;
}

int LEDScreen::squarePWM(float width) // map Width to be 0 <-> 2.0f* (float)M_PI
{
	if (m_phase < width)
	{
		return 1;
	}
	else
	{
		return 0;
	}
}
