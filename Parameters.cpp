#include "Parameters.h"
#include <cassert>

// Order of inputs is L, rho, T, r, loc, E, sigma0, sigma1
// Change this if you want to reorder inputs on the device
std::vector<std::string> Parameters::sNames = { "L", "rho", "T", "r", "loc", "E", "sigma0", "sigma1" };

Parameters::Parameters()
{
	// Order of params in DynamicDSS: L, rho, r, T, E, sigma0, sigma1
	// Do not change this
	std::vector<std::string> originalParameterOrder = { "L", "rho", "r", "T", "E", "sigma0", "sigma1" };
	for (int i = 0; i < originalParameterOrder.size(); i++)
	{
		const std::string& parameterName = originalParameterOrder[i];
		m_parameterIdMap[parameterName] = i;
	}

	// Setup internal state, cloning the DSS values
	m_parameters["L"] = Parameter(m_parameterIdMap["L"], 1.0f, { 0.5f, 4.0f });
	m_parameters["rho"] = Parameter(m_parameterIdMap["rho"], 7850.0f, { 15700.0f, 1962.5f}); // intentionally reversed range
	m_parameters["r"] = Parameter(m_parameterIdMap["r"], 0.0005f, { 0.0005f, 0.001f }); // intentionally reversed range
	m_parameters["T"] = Parameter(m_parameterIdMap["T"], 300.0f, { 150.f, 1200.0f });
	
	// This strange limit on the left hand dodges block-dropping without being perceptible
	m_parameters["E"] = Parameter(m_parameterIdMap["E"], 200000000000.f, { 5000000000.f, 400000000000.f });

	// Quadriplied Right side range
	m_parameters["sigma0"] = Parameter(m_parameterIdMap["sigma0"], 1.0f, { 0.f, 8.f, });

	// For that clean P I Z Z I C A T O 
	m_parameters["sigma1"] = Parameter(m_parameterIdMap["sigma1"], 0.005f, { 0.0002f, 0.04f });
	m_parameters["loc"] = Parameter(-1, -1.f, { 0.f, 1.f });
}

DynamicStiffString::SimulationParameters Parameters::getDSSParameters() const
{
	// Setup Dynamic Stiff String
	DynamicStiffString::SimulationParameters parameters = {};
	parameters.L = m_parameters.find("L")->second.getValue();
	parameters.rho = m_parameters.find("rho")->second.getValue();
	parameters.r = m_parameters.find("r")->second.getValue();
	parameters.T = m_parameters.find("T")->second.getValue();
	parameters.E = m_parameters.find("E")->second.getValue();
	parameters.sigma0 = m_parameters.find("sigma0")->second.getValue();
	parameters.sigma1 = m_parameters.find("sigma1")->second.getValue();
	return parameters;
}
