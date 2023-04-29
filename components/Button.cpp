#include "Button.h"

void Button::updateState(int newValue)
{
	m_previousState = m_state;
	m_state = static_cast<Button::State>(newValue);
}
