#pragma once

#ifdef DESKTOP_BUILD
#include "BelaMock.h"
#else
#include "Bela.h"
#endif

#include <vector>

class LEDScreen {
public:
	static constexpr int sCellCount = 10;

	LEDScreen(BelaContext* context);
	void setBrightness(const int cellIndex, const float brightness);
	void update(BelaContext* context);

private:
	void computePhase();
	int squarePWM(float width);

	std::vector<int> m_ledPins;
	float m_phase;
	float m_invSampleRate;
	float m_brightness[sCellCount] = {};


};