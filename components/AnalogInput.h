#pragma once

#ifdef DESKTOP_BUILD
#include "BelaMock.h"
#include "vector"
#else
#include "Bela.h"
#endif

#include <string>
#include <utility>

class AnalogInput {
public:
	AnalogInput(const int channel, const std::pair<float, float>& valueRange, const float readThreshold = 0.005f);

	void read(BelaContext* context, const int frame);

	float getCurrentValue() const { return m_currentValue; }
	float getCurrentValueMapped() const;
	float getCurrentValueinVolts() const;

	float mapValue(const float value) const;
	float maptoVolts(const float value) const;
	float unmapValue(const float value) const;

	bool hasChanged() { return m_hasChanged; }

	std::pair<float, float> getValueRange() { return std::make_pair(m_minValue, m_maxValue); }

private:
	int m_channel;
	const std::pair<float, float> m_valueRange;
	float m_readThreshold;

	bool m_hasChanged;
	float m_currentValue;
	
	// Calibration
	bool m_calibrate = 1;
	float m_maxValue;
	float m_minValue;
};
