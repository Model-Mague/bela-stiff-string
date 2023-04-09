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
	static constexpr int sLedPins[sCellCount] = {6, 7, 10, 2, 3, 0, 1, 4, 5, 8}; // Bela Pepper Pin Numbering

	LEDScreen(BelaContext* context);
	void setBrightness(const int cellIndex, const float brightness);
	void update(BelaContext* context);

private:
	void computePhase();
	int squarePWM(float width);

	float m_phase;
	float m_invSampleRate;
	float m_brightness[sCellCount] = {};


};