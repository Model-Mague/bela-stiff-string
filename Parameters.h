#pragma once
#include "DynamicStiffString/DynamicStiffString.h"

#include <map>
#include <string>

class Parameter {
public:
	Parameter() : m_value(0.f), m_range({0.f, 0.f}), m_id(-1) {}
	Parameter(const int id, const float value, const std::pair<float, float>& range) : m_id(id), m_value(value), m_range(range) {}
	void setValue(const float value) { m_value = value; }
	float getValue() const { return m_value; }
	const std::pair<float, float>& getRange() const { return m_range; }
	int getId() const { return m_id; }

private:
	float m_value;
	std::pair<float, float> m_range;
	int m_id; // ID to match it in DSS parameters
};

class Parameters {
public:
	/////////////////////////////////////////////////////////////////////////////////////////////////
	// Order of inputs is L, rho, T, r, loc, E, sigma0, sigma1
	// Change the order if you want to reorder inputs on the device
	/////////////////////////////////////////////////////////////////////////////////////////////////
	enum class Name : uint8_t {
		L,
		rho,
		T,
		r,
		loc,
		E,
		sigma0,
		sigma1
	};
	
	/////////////////////////////////////////////////////////////////////////////////////////////////
	using ParameterMap = std::map<Parameters::Name, Parameter>;
	/////////////////////////////////////////////////////////////////////////////////////////////////
	
	Parameters();
	Parameter& getParameter(const Parameters::Name name) { return m_parameters.find(name)->second; }
	DynamicStiffString::SimulationParameters getDSSParameters() const;
	ParameterMap getParameters() { return m_parameters; }

private:
	// Maintains an internal state of all the simulation parameters
	// These may be adjusted internally and don't necessarily match up with inputs
	ParameterMap m_parameters;
};