
#include "VirVubFitter/VirHelper.hh"

// << operator for vector<double>
std::ostream& operator<< (std::ostream& os, const std::vector<double>& v)
{
  for (int i(0); i < v.size(); ++i) {
    if (i > 0) os << " ";
    os << v[i];
  }
  return os;
}

