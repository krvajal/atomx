#pragma once

#include <istream>
#include <map>
#include <optional>
#include <string>

namespace atomx::io {

// A Fortran-style namelist group, e.g.
//   &input
//      zed=10., config='1s2 2s2 2p6', rel=0, dft='KLI',
//   /
// Keys are case-insensitive. Values are kept as text and converted on access.
class Namelist {
 public:
  // Reads the group named `group` (without '&'); throws if absent.
  static Namelist Read(std::istream& in, const std::string& group);
  static Namelist FromString(const std::string& text, const std::string& group);

  bool Has(const std::string& key) const;
  std::string GetString(const std::string& key, const std::string& fallback) const;
  double GetDouble(const std::string& key, double fallback) const;
  int GetInt(const std::string& key, int fallback) const;
  // Keys present in the input that were never read (likely typos).
  std::map<std::string, std::string> Unused() const;

 private:
  std::optional<std::string> Raw(const std::string& key) const;

  std::map<std::string, std::string> values_;
  mutable std::map<std::string, bool> used_;
};

}  // namespace atomx::io
