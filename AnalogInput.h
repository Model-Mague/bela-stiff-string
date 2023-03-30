#pragma once

#ifdef DESKTOP_BUILD
#include "BelaMock.h"
#else
#include "Bela.h"
#endif

#include <utility>

class AnalogInput {
public:
	AnalogInput(const int channel, const std::pair<float, float>& valueRange);

	float read(BelaContext* context, const int frame);
	float getCurrentValue() { return m_currentValue; }
	float hasChanged() { return m_currentValue != m_previousValue; }

	std::pair<float, float> getValueRange() { return std::make_pair(m_minValue, m_maxValue); }

private:
	int m_channel;
	const std::pair<float, float> m_valueRange;

	float m_currentValue;
	float m_previousValue;

	// Calibration
	bool m_calibrate = 0;
	float m_maxValue;
	float m_minValue;
};
