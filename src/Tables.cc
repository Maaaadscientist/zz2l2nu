#include <Tables.h>

#include <cmath>
#include <fstream>

namespace utils
{

bool table::record::belongTo(double pt, double eta) const
{
    return (pt < ptHi && pt >= ptLow) && (eta < etaHi && eta >= etaLow);
}

bool table::record::equalTo(int num) const
{
    return (std::abs(num - etaLow) < 0.5); // etaLow means the first value
}

table::table(const std::string &filename)
{
    std::ifstream file(filename);
    if (!file) {
        throw std::invalid_argument("File " + filename + " doesn't exist");
    }
    double data[7];
    while (file) {
        for (int i(0); i < 7; i++) {
            file >> data[i];
        }
        _recd.push_back(record{data[2], data[3], data[0], data[1], data[4], data[5], data[6]});
    }
}

double table::getEfficiency(double pt, double eta) const
{
    double hiPtBin = 0;
    for (unsigned int i = 0; i != _recd.size(); i++) {
        // if finds the proper bin, then return the efficiency
        if ((_recd[i]).belongTo(pt, eta)) return _recd[i].effi;
        // else store the average pt of the current bin efficency but do not return and try the next
        // bin
        if ((_recd[i]).belongTo(0.5 * (_recd[i].ptHi + _recd[i].ptLow), eta))
            hiPtBin = _recd[i].effi;
    }
    return hiPtBin;
}

double table::getEfficiencyLow(double pt, double eta) const
{
    double hiPtBin = 0;
    for (unsigned int i = 0; i != _recd.size(); i++) {
        if ((_recd[i]).belongTo(pt, eta)) return _recd[i].effi - _recd[i].effiErrorLow;
        if ((_recd[i]).belongTo(350, eta)) hiPtBin = _recd[i].effi;
    }
    return hiPtBin;
}

double table::getEfficiencyHigh(double pt, double eta) const
{
    double hiPtBin = 0;
    for (unsigned int i = 0; i != _recd.size(); i++) {
        if ((_recd[i]).belongTo(pt, eta)) return _recd[i].effi + _recd[i].effiErrorHigh;
        if ((_recd[i]).belongTo(350, eta)) hiPtBin = _recd[i].effi;
    }
    return hiPtBin;
}
} // namespace util
