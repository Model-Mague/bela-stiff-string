#include "AnalogInput.h"

#include <cmath>

AnalogInput::AnalogInput(const int channel, const std::pair<float, float>& valueRange)
	: m_channel(channel), m_valueRange(valueRange), m_hasChanged(false), m_currentValue(0.f), m_maxValue(-1000.f), m_minValue(1000.f)
{
}

float AnalogInput::read(BelaContext* context, const int frame)
{
	float read = analogRead(context, frame, m_channel);
	if (fabs(read - m_currentValue) > 0.005f)
	{
		m_hasChanged = true;
		m_currentValue = read;
	}
	else
	{
		m_hasChanged = false;
	}

	if (m_calibrate)
	{
		if (read > m_maxValue)
			m_maxValue = read;
		if (read < m_minValue)
			m_minValue = read;
	}

	return map(read, 0.f, 1.f, m_valueRange.first, m_valueRange.second);
}
