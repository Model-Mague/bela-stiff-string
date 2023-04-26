#include "Parameters.h"
#include <cassert>


Parameters::Parameters()
{
	// Order of params in DynamicDSS: L, rho, r, T, E, sigma0, sigma1
	// Do not change this
	const std::vector<Parameters::Name> originalParameterOrder = { Name::L, Name::rho, Name::r, Name::T, Name::E, Name::sigma0, Name::sigma1 };

	// Lookup function
	auto fnGetParamId = [&originalParameterOrder](Parameters::Name name) {
		const auto it = std::find(originalParameterOrder.begin(), originalParameterOrder.end(), name);
		return static_cast<int>(it - originalParameterOrder.begin());
	};

	// Setup internal state, cloning the DSS values
	m_parameters[Name::L] = Parameter(fnGetParamId(Name::L), 1.0f, {0.5f, 4.0f});
	m_parameters[Name::rho] = Parameter(fnGetParamId(Name::rho), 7850.0f, { 15700.0f, 1962.5f}); // intentionally reversed range
	m_parameters[Name::r] = Parameter(fnGetParamId(Name::r), 0.0005f, { 0.0005f, 0.001f }); // intentionally reversed range
	m_parameters[Name::T] = Parameter(fnGetParamId(Name::T), 300.0f, { 150.f, 1200.0f });
	
	// This strange limit on the left hand dodges block-dropping without being perceptible
	m_parameters[Name::E] = Parameter(fnGetParamId(Name::E), 200000000000.f, { 5000000000.f, 400000000000.f });

	// Quadriplied Right side range
	m_parameters[Name::sigma0] = Parameter(fnGetParamId(Name::sigma0), 1.0f, { 0.f, 8.f, });

	// For that clean P I Z Z I C A T O 
	m_parameters[Name::sigma1] = Parameter(fnGetParamId(Name::sigma1), 0.005f, { 0.0002f, 0.04f });
	m_parameters[Name::loc] = Parameter(-1, -1.f, { 0.f, 1.f });
}

DynamicStiffString::SimulationParameters Parameters::getDSSParameters() const
{
	// Setup Dynamic Stiff String
	DynamicStiffString::SimulationParameters parameters = {};
	parameters.L = m_parameters.find(Name::L)->second.getValue();
	parameters.rho = m_parameters.find(Name::rho)->second.getValue();
	parameters.r = m_parameters.find(Name::r)->second.getValue();
	parameters.T = m_parameters.find(Name::T)->second.getValue();
	parameters.E = m_parameters.find(Name::E)->second.getValue();
	parameters.sigma0 = m_parameters.find(Name::sigma0)->second.getValue();
	parameters.sigma1 = m_parameters.find(Name::sigma1)->second.getValue();
	return parameters;
}
