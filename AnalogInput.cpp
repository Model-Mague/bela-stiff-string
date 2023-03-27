#include "AnalogInput.h"


AnalogInput::AnalogInput(const int channel) : m_channel(channel), m_currentValue(0.f), m_previousValue(0.f)
{

}

float AnalogInput::read(BelaContext* context, const int frame)
{
	m_previousValue = m_currentValue;
	m_currentValue = analogRead(context, frame, m_channel);
	return m_currentValue;
}
