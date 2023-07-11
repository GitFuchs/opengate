
#ifndef GateAMDMActor_lookUpTable_h
#define GateAMDMActor_lookUpTable_h 1

#include <map>
#include <vector>
#include "G4UImanager.hh"

class lookUpTable
{
private:
    // int m_rows;
    // int m_columns;

    // empty map container
    // std::map<int, std::vector<std::vector<double>>> **m_delta_e;
    // std::map<int, std::vector<std::vector<double>>> **m_yd;
    std::vector<std::vector<double>> data;
    int key_columns;
    int charge_column;
    // int a_mass_column;
    // int let_column;
    int energy_column;
    int min_number_of_columns;

    // map key is charge, value is 2D vector energy/LET with values

public:
    // lookUpTable(); // Default Constructor - Null or Empty Matrix
    // lookUpTable(const int rows, const int columns);                 // Empty Matrix With Defined Size
    // lookUpTable(const int rows, const int columns, const T **data); // Defined Matrix

    // std::map<int, std::vector<std::vector<double>>> **getData() const;
    lookUpTable();
    ~lookUpTable();
    void printLookUpTable();
    void readLookUpTable(const std::string &strFilename);
    void determineRowsAndColumnsFromFile(const std::string &strFilename, int &rows, int &columns);
    std::vector<double> interpolateVectors(double targetValue, int column_to_use, std::vector<double> vector1, std::vector<double> vector2);

    double interpolateValue(double target, double x_1, double x_2, double f_1, double f_2);
    // bool findEntriesByLet(int charge, double let, std::vector<double> &result);
    bool findEntriesByEnergy(int charge, double energy, std::vector<double> &result);
    bool verifyLUT();

    // bool findEntries(int charge, int atomic_mass, double let, std::vector<double> &result);

}; // lookUpTable

#endif