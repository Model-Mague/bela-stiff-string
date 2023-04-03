#pragma once
#define _USE_MATH_DEFINES
#include <cstdint>
#include <cstdarg>
#include <iostream>
#include <sstream>

#define AUDIO_FRAMES 16
#define ANALOG_FRAMES 8

struct BelaContext {
	float  audioSampleRate = 44100.f;
	uint32_t  audioFrames = AUDIO_FRAMES;
	uint32_t  analogFrames = ANALOG_FRAMES;
	uint32_t digitalFrames = 16;
	const uint32_t  audioInChannels = 8;
	const uint32_t audioOutChannels = 2;
	const uint32_t  analogInChannels = 8;
	const uint32_t  analogOutChannels = 8;
	float audioIn[AUDIO_FRAMES * 2] = {}; // 2 channels
	float audioOut[AUDIO_FRAMES * 2] = {}; // 4 frames * 2 channels

	// Analog I/O runs at half the sample rate of audio I/O
	float analogIn[ANALOG_FRAMES * 8] = {}; // 8 channels
	float analogOut[ANALOG_FRAMES * 8] = {}; // 8 channels

	/// \brief Number of elapsed audio frames since the start of rendering.
	///
	/// This holds the total number of audio frames as of the beginning of the current block. To
	/// find the current number of analog or digital frames elapsed, multiply by the ratio of the
	/// sample rates (e.g. half the number of analog frames will have elapsed if the analog sample
	/// rate is 22050).
	uint64_t audioFramesElapsed = 0;
};

/**
 * \brief Read an analog input, specifying the frame number (when to read) and the channel.
 *
 * This function returns the value of an analog input, at the time indicated by \c frame.
 * The returned value ranges from 0 to 1, corresponding to a voltage range of 0 to 4.096V.
 *
 * \param context The I/O data structure which is passed by Bela to render().
 * \param frame Which frame (i.e. what time) to read the analog input. Valid values range
 * from 0 to (context->analogFrames - 1).
 * \param channel Which analog input to read. Valid values are between 0 and
 * (context->analogInChannels - 1), typically 0 to 7 by default.
 * \return Value of the analog input, range 0 to 1.
 */
static inline float analogRead(BelaContext* context, int frame, int channel) { return rand() > (RAND_MAX - 10); }

/**
 * \brief Write an audio output, specifying the frame number (when to write) and the channel.
 *
 * This function sets the value of an audio output, at the time indicated by \c frame. Valid
 * values are between -1 and 1.
 *
 * \param context The I/O data structure which is passed by Bela to render().
 * \param frame Which frame (i.e. what time) to write the audio output. Valid values range
 * from 0 to (context->audioFrames - 1).
 * \param channel Which analog output to write. Valid values are between 0 and
 * (context->audioChannels - 1), typically 0 to 1 by default.
 * \param value Value to write to the output, range -1 to 1.
 */
static inline void audioWrite(BelaContext* context, int frame, int channel, float value) {}

/**
 * \brief Write an analog output, specifying the frame number (when to write) and the channel.
 *
 * This function sets the value of an analog output, at the time indicated by \c frame. Valid
 * values are between 0 and 1, corresponding to the range 0 to 5V.
 *
 * Unlike analogWrite(), the value written will affect \b only the frame specified, with
 * future values unchanged. This is faster than analogWrite() so is better suited
 * to applications where every frame will be written to a different value. If
 * BELA_FLAG_ANALOG_OUTPUTS_PERSIST is not set within context->flags, then
 * analogWriteOnce() and analogWrite() are equivalent.
 *
 * \param context The I/O data structure which is passed by Bela to render().
 * \param frame Which frame (i.e. what time) to write the analog output. Valid values range
 * from 0 to (context->analogFrames - 1).
 * \param channel Which analog output to write. Valid values are between 0 and
 * (context->analogOutChannels - 1), typically 0 to 7 by default.
 * \param value Value to write to the output, range 0 to 1.
 */
static inline void analogWriteOnce(BelaContext* context, int frame, int channel, float value) {}

/**
 * \brief Read a digital input, specifying the frame number (when to read) and the pin.
 *
 * This function returns the value of a digital input, at the time indicated by \c frame.
 * The value is 0 if the pin is low, and nonzero if the pin is high (3.3V).
 *
 * \param context The I/O data structure which is passed by Bela to render().
 * \param frame Which frame (i.e. what time) to read the digital input. Valid values range
 * from 0 to (context->digitalFrames - 1).
 * \param channel Which digital input to read. 16 pins across the headers are
 * available. Check your board diagram to know where they are on the specific
 * board you have.
 * \return Value of the digital input.
 */
static inline int digitalRead(BelaContext* context, int frame, int channel) { return rand() > (RAND_MAX-1); }

/**
 * \brief Write a digital output, specifying the frame number (when to write) and the pin.
 *
 * This function sets the value of a digital output, at the time indicated by \c frame.
 * A value of 0 sets the pin low; any other value sets the pin high (3.3V).
 *
 * Unlike digitalWrite(), the value written will affect \b only the frame specified, with
 * future values unchanged. This is faster than digitalWrite() so is better suited
 * to applications where every frame will be written to a different value.
 *
 * \param context The I/O data structure which is passed by Bela to render().
 * \param frame Which frame (i.e. what time) to write the digital output. Valid values range
 * from 0 to (context->digitalFrames - 1).
 * \param channel Which digital output to write. 16 pins across the headers are
 * available. Check your board diagram to know where they are on the specific
 * board you have.
 * \param value Value to write to the output.
 */
static inline void digitalWriteOnce(BelaContext* context, int frame, int channel, int value) {}


class Scope {
public:
	Scope() : m_nChannels(0) {}
	void setup(unsigned int numChannels, float sampleRate) { m_nChannels = numChannels; }
	void log(const float* values)
	{
		std::stringstream ss;
		for (unsigned int i = 0; i < m_nChannels; i++) {
			ss << i << ": " << values[i];
			if (i != (m_nChannels - 1)) ss << " | ";
		}
		std::cout << ss.str() << std::endl;
	}

	void log(double chn1, ...)
	{
		std::stringstream ss;

		va_list args;
		va_start(args, chn1);
		for (unsigned int i = 0; i < m_nChannels; i++) {
			ss << i << ": " << (float)va_arg(args, double);
			if (i != (m_nChannels - 1)) ss << " | ";
		}
		va_end(args);

		std::cout << ss.str() << std::endl;
	}
private:
	unsigned int m_nChannels;
};

// map()
//
// Scale an input value from one range to another. Works like its Wiring language equivalent.
// x is the value to scale; in_min and in_max are the input range; out_min and out_max
// are the output range.

static inline float map(float x, float in_min, float in_max, float out_min, float out_max)
{
	return (x - in_min) * (out_max - out_min) / (in_max - in_min) + out_min;
}

static int rt_printf(const char* format, ...)
{
	va_list arg;
	int done;

	va_start(arg, format);
	done = vfprintf(stdout, format, arg);
	va_end(arg);

	return done;
}
