#include "AnalogInput.h"
#include "DynamicStiffString/Global.h"

#include <cmath>
#include <vector>

// Now incoming values are limited by imposed limits -> we loose an unperceptible amount of tweaking but guarantees the range limits are reachable
static std::vector<float> sUpperLimitValue = { 0.89f, 0.89f, 0.87f, 0.89f, 0.92f, 0.92f, 0.91f, 0.92f };
static std::vector<float> sLowerLimitValue = { 0.01f, 0.01f, 0.01f, 0.01f, 0.01f, 0.01f, 0.01f, 0.01f };

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
	return mapValue(m_currentValue);
}

float AnalogInput::getCurrentValueinVolts() const
{
	return maptoVolts(m_currentValue);
}

// Maps (the *resistor tolerance calibrated* digitalized 0-1) to the parameter's range within the DSS engine.

float AnalogInput::mapValue(const float value) const
{
	return map(Global::limit(value, sLowerLimitValue[m_channel], sUpperLimitValue[m_channel]), sLowerLimitValue[m_channel], sUpperLimitValue[m_channel], m_valueRange.first, m_valueRange.second);
}

// Unmaps to the 0-1 digitalized Analog Input Range *resistor tolerance calibrated

float AnalogInput::unmapValue(const float value) const
{
	return map(value, m_valueRange.first, m_valueRange.second, sLowerLimitValue[m_channel], sUpperLimitValue[m_channel]);

}

// Maps (the digitalized 0-1) to the real analog value in volts 0-10V * resistor tolerance calibrated

float AnalogInput::maptoVolts(const float value) const
{
	return Global::limit(value, sLowerLimitValue[m_channel], sUpperLimitValue[m_channel]) * 10.f;
}
