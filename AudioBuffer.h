#pragma once
#include <vector>

// Circular buffer
class AudioBuffer {
public:
	AudioBuffer(const int size);

	void put(float value);
	bool containsSilence() const;

private:
	unsigned int m_index;
	int m_silenceCounter;
	std::vector<float> m_buffer;
};