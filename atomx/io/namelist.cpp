#include "atomx/io/namelist.hpp"

#include <algorithm>
#include <cctype>
#include <sstream>

#include "atomx/core/types.hpp"

namespace atomx::io {

namespace {

std::string Lower(std::string s) {
  std::transform(s.begin(), s.end(), s.begin(),
                 [](unsigned char c) { return std::tolower(c); });
  return s;
}

std::string Trim(const std::string& s) {
  const auto b = s.find_first_not_of(" \t\r\n");
  if (b == std::string::npos) return "";
  const auto e = s.find_last_not_of(" \t\r\n");
  return s.substr(b, e - b + 1);
}

// Fortran real literals may use d/D exponents.
double ParseReal(std::string v) {
  std::replace(v.begin(), v.end(), 'd', 'e');
  std::replace(v.begin(), v.end(), 'D', 'e');
  return std::stod(v);
}

}  // namespace

Namelist Namelist::FromString(const std::string& text,
                              const std::string& group) {
  std::istringstream in(text);
  return Read(in, group);
}

Namelist Namelist::Read(std::istream& in, const std::string& group) {
  // Collect the body between "&group" and the closing "/" outside quotes.
  std::string all((std::istreambuf_iterator<char>(in)),
                  std::istreambuf_iterator<char>());
  const std::string lower = Lower(all);
  const auto start = lower.find("&" + Lower(group));
  if (start == std::string::npos) {
    throw Error("namelist &" + group + " not found");
  }
  std::size_t i = start + group.size() + 1;
  std::string body;
  char quote = 0;
  for (; i < all.size(); ++i) {
    const char c = all[i];
    if (quote) {
      if (c == quote) quote = 0;
      body += c;
      continue;
    }
    if (c == '\'' || c == '"') {
      quote = c;
      body += c;
      continue;
    }
    if (c == '!') {  // comment to end of line
      while (i < all.size() && all[i] != '\n') ++i;
      body += '\n';
      continue;
    }
    if (c == '/') break;
    body += c;
  }
  if (i >= all.size()) throw Error("namelist &" + group + " is not closed");

  // Split into key = value items on top-level commas and newlines.
  Namelist nl;
  std::vector<std::string> items;
  std::string cur;
  quote = 0;
  for (const char c : body) {
    if (quote) {
      if (c == quote) quote = 0;
      cur += c;
    } else if (c == '\'' || c == '"') {
      quote = c;
      cur += c;
    } else if (c == ',' || c == '\n') {
      items.push_back(cur);
      cur.clear();
    } else {
      cur += c;
    }
  }
  items.push_back(cur);
  for (const auto& raw : items) {
    const std::string item = Trim(raw);
    if (item.empty()) continue;
    const auto eq = item.find('=');
    if (eq == std::string::npos) {
      throw Error("cannot parse namelist item '" + item + "'");
    }
    std::string key = Lower(Trim(item.substr(0, eq)));
    std::string value = Trim(item.substr(eq + 1));
    if (value.size() >= 2 && (value.front() == '\'' || value.front() == '"') &&
        value.back() == value.front()) {
      value = value.substr(1, value.size() - 2);
    }
    nl.values_[key] = value;
  }
  return nl;
}

std::optional<std::string> Namelist::Raw(const std::string& key) const {
  const auto it = values_.find(Lower(key));
  if (it == values_.end()) return std::nullopt;
  used_[it->first] = true;
  return it->second;
}

bool Namelist::Has(const std::string& key) const {
  return values_.count(Lower(key)) > 0;
}

std::string Namelist::GetString(const std::string& key,
                                const std::string& fallback) const {
  return Raw(key).value_or(fallback);
}

double Namelist::GetDouble(const std::string& key, double fallback) const {
  const auto v = Raw(key);
  if (!v) return fallback;
  try {
    return ParseReal(*v);
  } catch (const std::exception&) {
    throw Error("'" + key + "' is not a number: " + *v);
  }
}

int Namelist::GetInt(const std::string& key, int fallback) const {
  const auto v = Raw(key);
  if (!v) return fallback;
  try {
    return std::stoi(*v);
  } catch (const std::exception&) {
    throw Error("'" + key + "' is not an integer: " + *v);
  }
}

std::map<std::string, std::string> Namelist::Unused() const {
  std::map<std::string, std::string> out;
  for (const auto& [k, v] : values_) {
    if (!used_.count(k)) out[k] = v;
  }
  return out;
}

}  // namespace atomx::io
