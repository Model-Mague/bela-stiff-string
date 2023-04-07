#pragma once

#ifdef DESKTOP_BUILD
#include "BelaMock.h"
#else
#include "Bela.h"
#endif

#include <string>
#include <utility>

class AnalogInput {
public:
	AnalogInput(const std::string& label, const int channel, const std::pair<float, float>& valueRange, const float readThreshold = 0.005f);

	void read(BelaContext* context, const int frame);

	float getCurrentValue() const { return m_currentValue; }
	float getCurrentValueMapped() const;

	bool hasChanged() { return m_hasChanged; }

	const std::string& getLabel() const { return m_label; }
	std::pair<float, float> getValueRange() { return std::make_pair(m_minValue, m_maxValue); }

private:
	const std::string m_label;
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
