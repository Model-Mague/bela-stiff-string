#pragma once
#include "DynamicStiffString/DynamicStiffString.h"
#include "../components/AnalogInput.h"

#include <map>
#include <string>
#include <memory>


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

enum class ParameterBehaviour : uint8_t {
	None,
	Pitch,
	Correction,
	Spray
};

class Parameter {
public:

	Parameter(const ParameterName name, 
		const float value, 
		const std::pair<float, float>& range, 
		const std::pair<ParameterBehaviour, float>& behaviour);

	float getValue() const { return m_value; }
	void setValue(const float value) { m_value = value; }
	
	ParameterName getName() const { return m_name; }
	ParameterBehaviour getBehaviour() const { return m_behaviour; }

	const std::pair<float, float>& getRange() const { return m_range; }
	int getId() const { return m_id; }
	std::shared_ptr<AnalogInput> getAnalogInput() { return m_analogInput; }
	int getChannel() const { return static_cast<int>(m_name); }

	// Pitch-Behaviour Only

	float getpitchRatio() const { return m_pitchRatio; }
	void setpitchRatio(const float ratio) { m_pitchRatio = ratio; }

	void calcOctaves();
	float getOctaves() const { return m_octaves; }

	float Volt_perOctave();

	bool activate1VMode() { m_1Vactive = true; }
	bool deactivate1Vmode() { m_1Vactive = false; }
	bool is1Vmodeactive() const { return m_1Vactive; }

private:
	ParameterName m_name;
	ParameterBehaviour m_behaviour;

	float m_value;
	std::pair<float, float> m_range;
	int m_id; // ID to match it in DSS parameters
	std::shared_ptr<AnalogInput> m_analogInput;

	// Pitch-behaviour Only

	float m_pitchRatio;
	float m_octaves;
	bool m_1Vactive;
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
