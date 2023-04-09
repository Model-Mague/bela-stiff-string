#include "AnalogInput.h"
//#include "Global.h"
#include <cmath>
#include <vector>

static float limit(float val, float min, float max) // Global.h could not be found? Ended up re-defining
    {
        if (val < min)
        {
            val = min;
            return val;
        }
        else if (val > max)
        {
            val = max;
            return val;
        }
        return val;
    }


std::vector<float> upperlimitValue = {0.89, 0.89, 0.87, 0.89, 0.92, 0.92, 0.91, 0.92};
std::vector<float> lowerlimitValue = {0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01};

AnalogInput::AnalogInput(const std::string& label, const int channel, const std::pair<float, float>& valueRange, const float readThreshold)
	: m_label(label), m_channel(channel), m_valueRange(valueRange), m_readThreshold(readThreshold),
	m_hasChanged(false), m_currentValue(0.f), m_maxValue(-1000.f), m_minValue(1000.f)
{
}

void AnalogInput::read(BelaContext* context, const int frame)
{
	float read = analogRead(context, frame, m_channel);
	if (fabs(read - m_currentValue) > m_readThreshold)
	{
		m_hasChanged = true;
		m_currentValue = read;
		/* <-- I planned to update the values as I went but I have realized exactitude does not pay off.
		if(read > upperlimitValue[m_channel])
			upperlimitValue[m_channel] = read;
			//rt_printf("updated channel %d upper limit with value %f\n" , m_channel, upperlimitValue[m_channel]);
		if(read < lowerlimitValue[m_channel])
			lowerlimitValue[m_channel] = read;
			//rt_printf("updated channel %d lower limit with value %f\n" , m_channel, lowerlimitValue[m_channel]);
		*/	
	}
	else
	{
		m_hasChanged = false;
	}

	if (m_calibrate)
	{
		if (read > m_maxValue)
			m_maxValue = read;
		if (read < m_minValue)
			m_minValue = read;
	}
	
}

float AnalogInput::getCurrentValueMapped() const
{
	return map(limit(m_currentValue, lowerlimitValue[m_channel], upperlimitValue[m_channel]), lowerlimitValue[m_channel], upperlimitValue[m_channel], m_valueRange.first, m_valueRange.second);
} // Now incoming values are limited by imposed limits -> we loose an unperceptible amount of tweaking but guarantees the range limits are reachable
