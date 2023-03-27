#pragma once

#ifdef DESKTOP_BUILD
#include "BelaMock.h"
#else
#include "Bela.h"
#endif

class AnalogInput {
public:
	AnalogInput(const int channel);

	float read(BelaContext* context, const int frame);
	float getCurrentValue() { return m_currentValue; }
	float hasChanged() { return m_currentValue != m_previousValue; }

private:
	int m_channel;
	float m_currentValue;
	float m_previousValue;
};
