#pragma once

#include <string>
#include <fstream>
#include <map>
#include <cstdlib>

template<typename T> 
class CSVWriter {

public:
	static CSVWriter& getInstance()
	{
		static CSVWriter instance;
		return instance;
	}

	void addNamedRow(const std::string& rowName) { m_rows.emplace(std::make_pair(rowName, std::vector<T>())); }
	std::vector<T>& getNamedRow(const std::string& rowName) { return m_rows[rowName]; }
	void addValue(const std::string& rowName, const T value) { m_rows[rowName].push_back(value); }

	void write(const std::string& filename)
	{
		std::ofstream out(filename);
		for (const auto& row : m_rows)
		{
			out << row.first;
			for (const auto val : row.second)
			{
				out << "," << val;
			}
			out << std::endl;
		}
	}

private:
	CSVWriter() {}
	std::map<std::string, std::vector<T>> m_rows;

};