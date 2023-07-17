#pragma once

#include <memory>

class Button {
public:
	enum class State : int {
		OFF = 0,
		ON = 1
	};
	enum class Type : size_t {
		TRIGGER = 0,
		MODE = 1,
		SPRAY = 2,
		DOWN = 3
	};

	Button(const int channel = -1): m_state(State::OFF), m_previousState(State::OFF), m_channel(channel) {}
	bool isReleased() const { return m_previousState == State::ON && m_state == State::OFF; }
	bool isPressed() const { return m_state == State::ON; }
	int getChannel() const { return m_channel; };
	void updateState(int newValue);

private:
	State m_state;
	State m_previousState;
	int m_channel;
};
