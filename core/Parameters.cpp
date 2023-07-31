#include "Parameters.h"

#include <algorithm>
#include <cassert>


Parameters::Parameters()
{
	// Convenience function for addding a parameter to map
	auto fnCreateParameter = [&](const ParameterName name, const float initValue, const std::pair<float, float>& range) {
		m_parameters.emplace(std::make_pair(name, Parameter(name, initValue, range)));
	};

	// Setup internal state, cloning the DSS values
	fnCreateParameter(ParameterName::L, 1.0f, { 0.125f, 2.0f });
	fnCreateParameter(ParameterName::rho, 7850.0f, { 15700.0f, 3925.f }); // intentionally reversed range
	fnCreateParameter(ParameterName::r, 0.0005f, { 0.001f, 0.00025f }); // intentionally reversed range
	fnCreateParameter(ParameterName::T, 300.0f, { 75.f, 600.0f });
	fnCreateParameter(ParameterName::E, 2e11, { 1e9, 4e13 });
	fnCreateParameter(ParameterName::sigma0, 1.0f, { 0.f, 2.f, });
	fnCreateParameter(ParameterName::sigma1, 1.f, { 0.0008f, 1.f });
	fnCreateParameter(ParameterName::loc, -1.f, { 0.f, 1.f });
}

DynamicStiffString::SimulationParameters Parameters::getDSSParameters() const
{
	// Setup Dynamic Stiff String
	DynamicStiffString::SimulationParameters parameters = {};
	parameters.L = m_parameters.find(ParameterName::L)->second.getValue();
	parameters.rho = m_parameters.find(ParameterName::rho)->second.getValue();
	parameters.r = m_parameters.find(ParameterName::r)->second.getValue();
	parameters.T = m_parameters.find(ParameterName::T)->second.getValue();
	parameters.E = m_parameters.find(ParameterName::E)->second.getValue();
	parameters.sigma0 = m_parameters.find(ParameterName::sigma0)->second.getValue();
	parameters.sigma1 = m_parameters.find(ParameterName::sigma1)->second.getValue();
	return parameters;
}

Parameter::Parameter(const ParameterName name, const float value, const std::pair<float, float>& range)
	: m_name(name), m_value(value), m_range(range)
{
	// Order of params in DynamicDSS: L, rho, r, T, E, sigma0, sigma1
	// Do not change this
	const std::vector<ParameterName> originalParameterOrder = {
		ParameterName::L,
		ParameterName::rho,
		ParameterName::r,
		ParameterName::T,
		ParameterName::E,
		ParameterName::sigma0,
		ParameterName::sigma1
	};

	// Lookup function
	auto fnGetParamId = [&originalParameterOrder](ParameterName name) {
		const auto it = std::find(originalParameterOrder.begin(), originalParameterOrder.end(), name);
		return it == originalParameterOrder.end() ? -1 : static_cast<int>(std::distance(originalParameterOrder.begin(), it));
	};
	m_id = fnGetParamId(name);

	const int channel = static_cast<int>(name);
	m_analogInput = std::make_shared<AnalogInput>(channel, range);
}
