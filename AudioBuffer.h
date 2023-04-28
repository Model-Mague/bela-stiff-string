#pragma once
#include <vector>

// Circular buffer
class AudioBuffer {
public:
	AudioBuffer(const int size);

	void put(float value);
	float get(const int index);
	bool containsSilence() const;

private:
	unsigned int m_head;
	int m_silenceCounter;
	std::vector<float> m_buffer;
};