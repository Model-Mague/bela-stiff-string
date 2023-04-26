#pragma once
#include "DynamicStiffString/DynamicStiffString.h"
#include "AnalogInput.h"

#include <map>
#include <string>


/////////////////////////////////////////////////////////////////////////////////////////////////
// Order of inputs is L, rho, T, r, loc, E, sigma0, sigma1
// Change the order if you want to reorder inputs on the device
/////////////////////////////////////////////////////////////////////////////////////////////////
enum class ParameterName : uint8_t {
	L,
	rho,
	T,
	r,
	loc,
	E,
	sigma0,
	sigma1
};

class Parameter {
public:
	Parameter(const ParameterName name, const float value, const std::pair<float, float>& range);

	float getValue() const { return m_value; }
	void setValue(const float value) { m_value = value; }

	ParameterName getName() const { return m_name; }
	const std::pair<float, float>& getRange() const { return m_range; }
	int getId() const { return m_id; }
	std::shared_ptr<AnalogInput> getAnalogInput() { return m_analogInput; }
	int getChannel() const { return static_cast<int>(m_name); }

private:
	ParameterName m_name;
	float m_value;
	std::pair<float, float> m_range;
	int m_id; // ID to match it in DSS parameters
	std::shared_ptr<AnalogInput> m_analogInput;
};

class Parameters {
public:
	/////////////////////////////////////////////////////////////////////////////////////////////////
	using ParameterMap = std::map<ParameterName, Parameter>;
	/////////////////////////////////////////////////////////////////////////////////////////////////
	
	Parameters();
	Parameter& getParameter(const ParameterName name) { return m_parameters.find(name)->second; }
	DynamicStiffString::SimulationParameters getDSSParameters() const;
	ParameterMap getParameters() { return m_parameters; }

private:
	// Maintains an internal state of all the simulation parameters
	// These may be adjusted internally and don't necessarily match up with inputs
	ParameterMap m_parameters;
};
