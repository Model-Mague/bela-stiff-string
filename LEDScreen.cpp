#include "LEDScreen.h"

LEDScreen::LEDScreen(BelaContext* context) : m_phase(0.f)
{
	m_ledPins = { 6, 7, 10, 2, 3, 0, 1, 4, 5, 8 }; // Bela Pepper Pin Numbering
	for (unsigned int i = 0; i < sCellCount; i++)
	{
		pinMode(context, 0, m_ledPins[i], OUTPUT);
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

		for (int i = 0; i < 2; i++) // 0-1
		{
			digitalWriteOnce(context, n, m_ledPins[i], squarePWM(m_brightness[0]));
		}
		for (int i = 2; i < 8; i++) // 2-7
		{
			digitalWriteOnce(context, n, m_ledPins[i], squarePWM(m_brightness[i - 1]));
		}
		for (int i = 8; i < sCellCount; i++) // 8-9
		{
			digitalWriteOnce(context, n, m_ledPins[i], squarePWM(m_brightness[7]));
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
