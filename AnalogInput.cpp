#include "AnalogInput.h"
#include "DynamicStiffString/Global.h"
#include <cmath>
#include <vector>

// Now incoming values are limited by imposed limits -> we loose an unperceptible amount of tweaking but guarantees the range limits are reachable
static std::vector<float> sUpperLimitValue = { 0.89f, 0.89f, 0.87f, 0.89f, 0.92f, 0.92f, 0.91f, 0.92f };
static std::vector<float> sLowerLimitValue = { 0.01f, 0.01f, 0.01f, 0.01f, 0.01f, 0.01f, 0.01f, 0.01f };

AnalogInput::AnalogInput(const std::string& label, const int channel, const std::pair<float, float>& valueRange, const float readThreshold)
	: m_label(label), m_channel(channel), m_valueRange(valueRange), m_readThreshold(readThreshold),
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
		/* <-- I planned to update the values as I went but I have realized exactitude does not pay off.
		if(read > sUpperLimitValue[m_channel])
			sUpperLimitValue[m_channel] = read;
			//rt_printf("updated channel %d upper limit with value %f\n" , m_channel, sUpperLimitValue[m_channel]);
		if(read < sLowerLimitValue[m_channel])
			sLowerLimitValue[m_channel] = read;
			//rt_printf("updated channel %d lower limit with value %f\n" , m_channel, sLowerLimitValue[m_channel]);
		*/	
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
	return mapValue(m_currentValue);
} 

float AnalogInput::mapValue(const float value) const
{
	return map(Global::limit(value, sLowerLimitValue[m_channel], sUpperLimitValue[m_channel]), sLowerLimitValue[m_channel], sUpperLimitValue[m_channel], m_valueRange.first, m_valueRange.second);
}

float AnalogInput::unmapValue(const float value) const
{
	return map(value, m_valueRange.first, m_valueRange.second, sLowerLimitValue[m_channel], sUpperLimitValue[m_channel]);

}
