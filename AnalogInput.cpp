#include "AnalogInput.h"


AnalogInput::AnalogInput(const int channel) : m_channel(channel), m_currentValue(0.f), m_previousValue(0.f), m_maxValue(-1000.f), m_minValue(1000.f)
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
	return m_currentValue;
}
