#include "AnalogInput.h"


AnalogInput::AnalogInput(const int channel, const std::pair<float, float>& valueRange)
	: m_channel(channel), m_currentValue(0.f), m_previousValue(0.f), m_maxValue(-1000.f), m_minValue(1000.f), m_valueRange(valueRange)
{
}

float AnalogInput::read(BelaContext* context, const int frame)
{
	m_previousValue = m_currentValue;
	m_currentValue = analogRead(context, frame, m_channel);
	if (m_calibrate)
	{
		if (m_currentValue > m_maxValue)
			m_maxValue = m_currentValue;
		if (m_currentValue < m_minValue)
			m_minValue = m_currentValue;
	}
	return map(m_currentValue, 0.f, 1.f, m_valueRange.first, m_valueRange.second);
}
