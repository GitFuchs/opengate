#include <algorithm>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>

#include "GateAMDMActor_lookUpTable.hh"

lookUpTable::lookUpTable() {
  // std::cout << "initialized lookUpTable" << std::endl;
  // how many columns are used as key
  key_columns = 2;
  charge_column = 0;
  // a_mass_column = -1;
  energy_column = 1;
  // let_column = -3;
  // there should be 20 data columns
  min_number_of_columns = key_columns + 20;
}
lookUpTable::~lookUpTable() {
  // std::cout << "deleted HistoManager" << std::endl;
}

void lookUpTable::readLookUpTable(const std::string &strFilename) {
  // int rows = 0;
  // int columns = 0;
  // determineRowsAndColumnsFromFile(strFilename, rows, columns);
  // std::cout << "found columns: " << columns << " found rows: " << rows <<
  // std::endl;

  // std::map<int, std::vector<std::vector<double>>> **m_delta_e;
  // std::map<int, std::vector<std::vector<double>>> **m_yd;
  // std::vector<std::vector<double>> **data;

  // Opening the file
  std::ifstream inputfile(strFilename);
  // check if file is open
  if (!inputfile.is_open()) {
    std::cout << "Error opening file: " << strFilename << std::endl;
    exit(-1);
  }

  // int i = 0;
  int j = 0;
  std::string line;
  // std::cout << "reading data from file" << std::endl;
  // std::cout << "number of vector before reading: " << data.size() <<
  // std::endl;

  while (std::getline(inputfile, line)) {
    if (line.size() && line[0] == '#')
      continue;
    std::istringstream iss(line);
    // iss << line; // send the line to the stringstream object...
    double value;
    // Define the inner vector
    std::vector<double> vectorRow;

    while (iss >> value) {
      vectorRow.push_back(value);
      j++;
    }
    data.push_back(vectorRow);
  }
  verifyLUT();
  // std::cout << data << std::endl;
}

bool lookUpTable::verifyLUT() {
  // check the number of rows
  if (data.size() <= 0) {
    std::cout << "ERROR: no data in look up table" << std::endl;
    return false;
  }
  // check if all rows have the same number of columns
  int columns = data[0].size();
  for (long unsigned int i = 0; i < data.size(); i++) {
    if (data[i].size() != columns) {
      std::cout << "ERROR: row " << i << " has " << data[i].size()
                << " columns, but should have " << columns << std::endl;
      return false;
    }
  }
  // check if the first column is sorted
  for (long unsigned int i = 1; i < data.size(); i++) {
    if (data[i][0] < data[i - 1][0]) {
      std::cout << "ERROR: row " << i << " has a smaller value than row "
                << i - 1 << std::endl;
      return false;
    }
  }
  // check if enough columns are present
  if (columns < min_number_of_columns) {
    std::cout << "ERROR: not enough columns in look up table. Number of "
                 "columns found: "
              << columns << std::endl;
    return false;
  }

  return true;
}

void lookUpTable::printLookUpTable() {
  std::cout << "printing look up table" << std::endl;
  std::cout << "number of rows: " << data.size() << std::endl;
  if (data.size() <= 0) {
    std::cout << "look up table is empty" << std::endl;
    return;
  } else {
    std::cout << "number of columns: " << data[0].size() << std::endl;
    for (long unsigned int i = 0; i < data.size(); i++) {
      for (long unsigned int j = 0; j < data[i].size(); j++)
        std::cout << data[i][j] << " ";
      std::cout << std::endl;
    }
  }
}

/* bool lookUpTable::findEntriesByLet(int charge, double let,
std::vector<double> &result)
{
    // returns 0 if successful
    // printLookUpTable();
    // index 0 contains charge
    // index 1 contains atomic mass
    // index 2 contains energy
    // index 3 contains LET

    // key_columns = 4;
    // charge_column = 0;
    // a_mass_column = 1;
    // energy_column = 2;
    // let_column = 3;

    // int key_columns = 4;
    long unsigned int i = 0;
    // std::vector<double>
    //     slice;
    // std::vector<double> slice(10, 0);
    result.clear();
    // if arriving here, return zero vector
    // std::cout << "zero entry" << std::endl;

    // Declaring an iterator
    // find range of iterators in which the charge is correct
    auto lower_charge = std::lower_bound(data.begin(), data.end(), charge,
[this](auto &vec, int value) { return vec[charge_column] < value; });

    auto upper_charge = std::upper_bound(data.begin(), data.end(), charge,
[this](int value, auto &vec) { return value < vec[charge_column]; });
    // if charge not found exit.
    if (lower_charge == data.end())
    {
        for (long unsigned int j = 0; j < data[0].size() - key_columns; j++)
        {
            result.push_back(0.);
            // std::cout << "j: " << j << "to acchieve " << data[i].size() <<
"size: " << (int)result.size() << "capacity: " << (int)result.capacity() <<
std::endl; return false;
        }
    }
    // int charge_tmp = *(lower_charge->begin() + 0);
    // double mass_tmp = *(lower_charge->begin() + 1);
    // double energy_tmp = *(lower_charge->begin() + 2);
    // double let_tmp = *(lower_charge->begin() + 3);
    // std::cout << "charge selection: value at postion of lower iterator is: "
    //           << " charge: " << charge_tmp
    //           << " mass: " << mass_tmp
    //           << " energy: " << energy_tmp
    //           << " let: " << let_tmp << std::endl;
    // charge_tmp = *(upper_charge->begin() + 0);
    // mass_tmp = *(upper_charge->begin() + 1);
    // energy_tmp = *(upper_charge->begin() + 2);
    // let_tmp = *(upper_charge->begin() + 3);
    // std::cout << "charge selection: value at postion of upper iterator is: "
    //           << " charge: " << charge_tmp
    //           << " mass: " << mass_tmp
    //           << " energy: " << energy_tmp
    //           << " let: " << let_tmp << std::endl;

    // check if no suitable charge was found, in this case lower und
upper_charge should be the same, e.g. at end. if (lower_charge == upper_charge)
    {
        std::cout << "WARNING: no suitable charge entry found for charge: " <<
charge << std::endl;
        // exit(-1);
        for (long unsigned int j = 0; j < data[0].size() - key_columns; j++)
        {
            result.push_back(0.);
            // std::cout << "j: " << j << "to acchieve " << data[i].size() <<
"size: " << (int)result.size() << "capacity: " << (int)result.capacity() <<
std::endl;
        }
        return false;
    }

    // ##############################
    // find LET
    // find range of iterators in which the LET is correct
    auto lower_let = std::lower_bound(lower_charge, upper_charge, let,
[this](auto &vec, double value) { return vec[let_column] < value; });
//     // auto upper_let = std::upper_bound(lower_charge, upper_charge, let,
[this](double value, auto &vec)
//     //                                   { return value < vec[let_column];
});
//
//     // std::cout << "target let: " << let << std::endl;
//     // auto charge_tmp = *(lower_let->begin() + 0);
//     // auto mass_tmp = *(lower_let->begin() + 1);
//     // auto energy_tmp = *(lower_let->begin() + 2);
//     // auto let_tmp = *(lower_let->begin() + 3);
//     // std::cout << "let selection: value at postion of lower iterator is: "
//     //           << " charge: " << charge_tmp
//     //           << " mass: " << mass_tmp
//     //           << " energy: " << energy_tmp
//     //           << " let: " << let_tmp << std::endl;
//     // charge_tmp = *(lower_let->begin() + charge_column);
//     // mass_tmp = *(lower_let->begin() + a_mass_column);
//     // energy_tmp = *(lower_let->begin() + 2);
//     // let_tmp = *(lower_let->begin() + let_column);
//     // std::cout << "let selection: value at postion of upper iterator is: "
//     //           << " charge: " << charge_tmp
//     //           << " mass: " << mass_tmp
//     //           << " energy: " << energy_tmp
//     //           << " let: " << let_tmp
//     //           << " target let: " << let << std::endl;
//     // exit(-1);

    // check if end of lut was reached, if so go to last entry which is in list
    if (lower_let >= data.end())
    {
        lower_let = std::prev(lower_let);
    }

    // if precisely matching, return value entries
    if (*(lower_let->begin() + let_column) == let)
    {
        // convert iterator to index
        i = std::distance(data.begin(), lower_let);
        std::copy(data[i].begin() + key_columns, data[i].end(),
std::back_inserter(result)); return true;
    }
    // // if let too low and not in list, return first one
    else if ((*(lower_charge->begin() + let_column) > let) && (lower_let ==
lower_charge))
    {
        // convert iterator to index
        i = std::distance(data.begin(), lower_charge);
        std::copy(data[i].begin() + key_columns, data[i].end(),
std::back_inserter(result)); return true;
    }
    // // if let too high and not in list, return last one
    else if ((*(std::prev(upper_charge)->begin() + let_column) < let) &&
(lower_let == upper_charge))
    {
        // convert iterator to index
        i = std::distance(data.begin(), std::prev(upper_charge));
        std::copy(data[i].begin() + key_columns, data[i].end(),
std::back_inserter(result)); return true;
    }
    // if not matching, interpolate
    else
    {
        // convert iterator to index
        // i = std::distance(data.begin(), upper_let);
        if (lower_let > lower_charge)
        {
            // std::cout << "interpolating target let: " << let
            //           << " first let: " << *((lower_let - 1)->begin() +
let_column)
            //           << " second let: " << *((lower_let)->begin() +
let_column)
            //           << " indexed let: " << *(lower_let->begin() +
let_column)
            //           << std::endl;

            // std::cout << "in, so range valid" << std::endl;
            result = interpolateVectors(let, let_column, *std::prev(lower_let),
*(lower_let)); return true;
        }
    }

    // just in case, if nothing was found
    // std::cout << "ERROR: should not end up here in findEntries" << std::endl;
    // exit(-1);
    for (long unsigned int j = 0; j < data[0].size() - key_columns; j++)
    {
        result.push_back(0.);
        // std::cout << "j: " << j << "to acchieve " << data[i].size() << "size:
" << (int)result.size() << "capacity: " << (int)result.capacity() << std::endl;
    }
    return false;
}
*/

bool lookUpTable::findEntriesByEnergy(int charge, double energy,
                                      std::vector<double> &result) {
  // returns 0 if successful

  // key_columns = 4;
  // charge_column = 0;
  // a_mass_column = 1;
  // energy_column = 2;
  // let_column = 3;

  // int key_columns = 4;
  long unsigned int i = 0;
  result.clear();

  // Declaring an iterator
  // find range of iterators in which the charge is correct
  auto lower_charge = std::lower_bound(
      data.begin(), data.end(), charge,
      [this](auto &vec, int value) { return vec[charge_column] < value; });

  auto upper_charge = std::upper_bound(
      data.begin(), data.end(), charge,
      [this](int value, auto &vec) { return value < vec[charge_column]; });
  // if charge not found exit.
  if (lower_charge == data.end()) {
    std::cout << "INFO: no suitable charge entry found for charge: " << charge
              << std::endl;
    for (long unsigned int j = 0; j < data[0].size() - key_columns; j++) {
      result.push_back(0.);
      // std::cout << "j: " << j << "to acchieve " << data[i].size() << "size: "
      // << (int)result.size() << "capacity: " << (int)result.capacity() <<
      // std::endl;
      return false;
    }
  }
  // int charge_tmp = *(lower_charge->begin() + 0);
  // double mass_tmp = *(lower_charge->begin() + 1);
  // double energy_tmp = *(lower_charge->begin() + 2);
  // double let_tmp = *(lower_charge->begin() + 3);
  // std::cout << "charge selection: value at postion of lower iterator is: "
  //           << " charge: " << charge_tmp
  //           << " mass: " << mass_tmp
  //           << " energy: " << energy_tmp
  //           << " let: " << let_tmp << std::endl;
  // charge_tmp = *(upper_charge->begin() + 0);
  // mass_tmp = *(upper_charge->begin() + 1);
  // energy_tmp = *(upper_charge->begin() + 2);
  // let_tmp = *(upper_charge->begin() + 3);
  // std::cout << "charge selection: value at postion of upper iterator is: "
  //           << " charge: " << charge_tmp
  //           << " mass: " << mass_tmp
  //           << " energy: " << energy_tmp
  //           << " let: " << let_tmp << std::endl;

  // check if no suitable charge was found, in this case lower und upper_charge
  // should be the same, e.g. at end.
  if (lower_charge == upper_charge) {
    // std::cout << "INFO: no suitable charge entry found for charge: " <<
    // charge << std::endl;
    for (long unsigned int j = 0; j < data[0].size() - key_columns; j++) {
      result.push_back(0.);
      // std::cout << "j: " << j << "to acchieve " << data[i].size() << "size: "
      // << (int)result.size() << "capacity: " << (int)result.capacity() <<
      // std::endl;
    }
    return false;
  }

  // ##############################
  // find range of iterators in which the charge is correct
  auto lower_energy = std::lower_bound(
      lower_charge, upper_charge, energy,
      [this](auto &vec, double value) { return vec[energy_column] < value; });

  // std::cout << "target energy: " << energy << std::endl;
  // auto charge_tmp = *(lower_charge->begin() + charge_column);
  // auto energy_tmp = *(lower_charge->begin() + energy_column);
  // std::cout << "energy selection: value at postion of lower iterator is: "
  //           << " charge: " << charge_tmp
  //           << " energy: " << energy_tmp << std::endl;
  // charge_tmp = *(lower_charge->begin() + charge_column);
  // mass_tmp = *(lower_charge->begin() + a_mass_column);
  // energy_tmp = *(lower_charge->begin() + 2);
  // std::cout << "energy selection: value at postion of upper iterator is: "
  //           << " charge: " << charge_tmp
  //           << " mass: " << mass_tmp
  //           << " energy: " << energy_tmp
  //           << " target energy: " << energy << std::endl;

  // safety check if end of lut was exceeded, if so go to last entry which is in
  // LUT
  if (lower_energy > data.end()) {
    lower_energy = data.end();
    std::cout
        << "WARNING: Look-up table: reseting index. This should not happen."
        << std::endl;
  }

  // if precisely matching, return value entries
  if (*(lower_energy->begin() + energy_column) == energy) {
    // convert iterator to index
    i = std::distance(data.begin(), lower_energy);
    std::copy(data[i].begin() + key_columns, data[i].end(),
              std::back_inserter(result));
    // std::cout << "precise match target energy: " << energy
    //           << " previous energy: " << *(std::prev(lower_energy)->begin() +
    //           energy_column)
    //           << " indexed energy: " << *(lower_energy->begin() +
    //           energy_column)
    //           << " returned energy: " << *(lower_energy->begin() +
    //           energy_column)
    //           << std::endl;
    return true;
  }
  // // if energy too low and not in list, return first one
  else if ((*(lower_charge->begin() + energy_column) > energy) &&
           (lower_energy == lower_charge)) {
    // convert iterator to index
    i = std::distance(data.begin(), lower_charge);
    std::copy(data[i].begin() + key_columns, data[i].end(),
              std::back_inserter(result));
    // std::cout << "TOO LOW target energy: " << energy
    //           //   << " previous energy: " <<
    //           *(std::prev(lower_energy)->begin() + energy_column)
    //           << " indexed energy: " << *(lower_energy->begin() +
    //           energy_column)
    //           << " returned energy: " << *(lower_energy->begin() +
    //           energy_column)
    //           << std::endl;
    return true;
  }
  // // if energy too high and not in list, return last one
  else if ((*(std::prev(upper_charge)->begin() + energy_column) < energy) &&
           (lower_energy == upper_charge)) {
    // convert iterator to index
    i = std::distance(data.begin(), std::prev(upper_charge));
    std::copy(data[i].begin() + key_columns, data[i].end(),
              std::back_inserter(result));
    // std::cout << "TOO HIGH target energy: " << energy
    //           << " previous energy: " << *(std::prev(lower_energy)->begin() +
    //           energy_column)
    //           << " indexed energy: " << *(lower_energy->begin() +
    //           energy_column)
    //           << " returned energy: " << *(upper_charge->begin() +
    //           energy_column)
    //           << std::endl;
    return true;
  }
  // if not matching, interpolate
  else {
    if (lower_energy > lower_charge) {
      // std::cout << "interpolating target energy: " << energy
      //           << " previous energy: " << *(std::prev(lower_energy)->begin()
      //           + energy_column)
      //           << " indexed energy: " << *(lower_energy->begin() +
      //           energy_column)
      //           << std::endl;

      // std::cout << "in, so range valid" << std::endl;
      result = interpolateVectors(energy, energy_column,
                                  *std::prev(lower_energy), *(lower_energy));
      return true;
    }
  }

  // if no suitable entry was found, return zero and false
  for (long unsigned int j = 0; j < data[0].size() - key_columns; j++) {
    result.push_back(0.);
    // std::cout << "j: " << j << "to acchieve " << data[i].size() << "size: "
    // << (int)result.size() << "capacity: " << (int)result.capacity() <<
    // std::endl;
  }
  std::cout
      << "NOT FOUND target energy: "
      << energy
      //   << " previous energy: " << *(std::prev(lower_energy)->begin() +
      //   energy_column)
      //   << " indexed energy: " << *(lower_energy->begin() + energy_column)
      //   << " returned energy: " << *(upper_charge->begin() + energy_column)
      << std::endl;
  return false;
}

std::vector<double>
lookUpTable::interpolateVectors(double targetValue, int column_to_use,
                                std::vector<double> vector1,
                                std::vector<double> vector2) {
  // passed vectors contain 3 leading columns
  //  index 0 contains charge
  // index 1 contains energy
  // index 2 contains LET

  std::vector<double> result;
  double resultValue;
  // loop over whole vector
  for (long unsigned int i = key_columns; i < vector1.size(); i++) {
    resultValue =
        interpolateValue(targetValue, vector1[column_to_use],
                         vector2[column_to_use], vector1[i], vector2[i]);
    result.push_back(resultValue);
  }
  return result;
}

double lookUpTable::interpolateValue(double target, double x_1, double x_2,
                                     double f_1, double f_2) {
  double result;

  result = f_1 + (f_2 - f_1) / (x_2 - x_1) * (target - x_1);
  // std::cout << "interpol target " << target << " x_1 " << x_1 << " x_2 " <<
  // x_2 << " f_1 " << f_1 << " f_2 " << f_2 << " result " << result <<
  // std::endl;
  if ((result <= f_2 && result >= f_1) || (result >= f_2 && result <= f_1)) {
  } else {
    std::cout << "+++++++++++++++++++Interpolation "
                 "WARNING++++++++++++++++++++++++++++++++++++"
              << std::endl;
    std::cout << "in interpolateValue: target " << target << " x_1 " << x_1
              << " x_2 " << x_2 << " f_1 " << f_1 << " f_2 " << f_2
              << " result " << result << std::endl;
  }
  return result;
}

void lookUpTable::determineRowsAndColumnsFromFile(
    const std::string &strFilename, int &rows, int &columns) {
  rows = 0;
  columns = 0;

  // Opening the file
  std::ifstream inputfile(strFilename);
  // check if file is open
  if (!inputfile.is_open()) {
    std::cout << "Error opening file: " << strFilename << " exiting."
              << std::endl;
    exit(-1);
  }

  // count columns in file
  std::string line;
  int max_columns = 0;
  int columns_in_line = 0;

  while (std::getline(inputfile, line)) {
    columns_in_line = 0;
    // std::getline(inputfile, line); // read the first line of your file to
    // string

    std::istringstream iss(line);
    // iss << line; // send the line to the stringstream object...

    double value;
    while (iss >> value)
      columns_in_line++; // while there's something in the line, increase the
                         // number of columns

    if (max_columns < columns_in_line)
      max_columns = columns_in_line;
    rows++; // increase number of rows by one
  }
  // close file
  inputfile.close();
  // // rewind inputfile
  // inputfile.clear();
  // inputfile.seekg(0);
  // columns = max_columns;
  // // std::cout << "found columns: " << columns_in_line << " found rows: " <<
  // rows << std::endl;

  // // // count rows in file
  // // rows = std::count(std::istreambuf_iterator<char>(inputfile),
  // //                   std::istreambuf_iterator<char>(), '\n');
  // // std::rewind(inputfile);
}
