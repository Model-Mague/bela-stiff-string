#include <Bela.h>
#include <cmath>
#include <libraries/math_neon/math_neon.h>
#include <vector>

float phase_;
float Width;
float out;
float invSampleRate_;
std::vector<int> ledPins = {6, 7, 10, 2, 3, 0, 1, 4, 5, 8}; // Bela Pepper Pin Numbering
std::vector<float> analogInput(8);

void computePhase()
{
	phase_ += 440 * invSampleRate_;
	if(phase_ >= 1)
		phase_ -= 1;
}

float squarePWM(float Width) // map Width to be 0 <-> 2.0f* (float)M_PI
{
	if(phase_ < Width)
	{
		out = 1;
	}
	else
	{
		out = 0;
	}
	return out;
}

bool setup(BelaContext *context, void *userData)
{
	for(unsigned int i = 0; i < 10; i++)
	{
		pinMode(context, 0, ledPins[i], OUTPUT);
	}
	
	invSampleRate_ = 1 / context->digitalSampleRate;

	return true;
}


void render(BelaContext *context, void *userData)
{
		for(unsigned int n = 0; n < context->analogFrames; n++) 
		{
			for(unsigned int m = 0; m < 8; m++)
			{
				analogInput[m] = analogRead(context, n, m);
			}
		
		}
		
		for(unsigned int n = 0; n < context->digitalFrames; n++)
		{
			computePhase();
			
			for(unsigned int m = 0; m < 2; m++)
			{
				digitalWriteOnce(context, n, ledPins[m], squarePWM(analogInput[0]));
			}
			
			for(unsigned int m = 2; m < 8; m++)
			{
				digitalWriteOnce(context, n, ledPins[m], squarePWM(analogInput[m - 1]));
			}
			
			for(unsigned int m = 8; m < 10; m++)
			{
				digitalWriteOnce(context, n, ledPins[m], squarePWM(analogInput[7]));
			}
			
		}
	
}

void cleanup(BelaContext *context, void *userData)
{
	
}
