#include "AnalogInput.h"

#include <cmath>

AnalogInput::AnalogInput(const int channel, const std::pair<float, float>& valueRange, const float readThreshold)
	: m_channel(channel), m_valueRange(valueRange), m_readThreshold(readThreshold),
	m_hasChanged(false), m_currentValue(0.f), m_maxValue(-1000.f), m_minValue(1000.f)
{
}

void AnalogInput::read(BelaContext* context, const int frame)
{
	float read = analogRead(context, frame, m_channel);
	if (fabs(read - m_currentValue) > m_readThreshold)
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
}

float AnalogInput::getCurrentValueMapped() const
{
	return map(m_currentValue, 0.f, 1.f, m_valueRange.first, m_valueRange.second);
}
