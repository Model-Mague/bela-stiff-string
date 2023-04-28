#include "AudioBuffer.h"
#include <cmath>

AudioBuffer::AudioBuffer(const int size)
	: m_index(0), m_silenceCounter(0)
{
	m_buffer.resize(size, 0.f);
}

void AudioBuffer::put(float value)
{
	const auto size = m_buffer.size();
	m_buffer[m_index++ % size] = value;
	if (value != 0.f) m_silenceCounter = static_cast<int>(size);
	else m_silenceCounter = std::max(0, m_silenceCounter - 1);
}

bool AudioBuffer::containsSilence() const
{
	return m_silenceCounter == 0;
}
