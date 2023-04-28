#include "AudioBuffer.h"
#include <cmath>

AudioBuffer::AudioBuffer(const int size)
	: m_head(0), m_silenceCounter(0)
{
	m_buffer.resize(size, 0.f);
}

void AudioBuffer::put(float value)
{
	const auto size = m_buffer.size();
	m_buffer[m_head++ % size] = value;
	if (value != 0.f)
	{
		m_silenceCounter = static_cast<int>(size);
	}
	else 
	{
		m_silenceCounter = std::max(0, m_silenceCounter - 1);
	}
}

float AudioBuffer::get(const int index)
{
	// If the head is at, say, index 2 then the oldest entry is at index 3
	// If the oldest entry is queried (index = 0) then we pick the 3+1=4th element 
	return m_buffer[(index + m_head + 1) % m_buffer.size()];
}

bool AudioBuffer::containsSilence() const
{
	return m_silenceCounter == 0;
}
