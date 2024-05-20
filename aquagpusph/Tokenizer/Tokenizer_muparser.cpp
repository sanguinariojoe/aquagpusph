/*
 *  This file is part of AQUAgpusph, a free CFD program based on SPH.
 *  Copyright (C) 2012  Jose Luis Cercos Pita <jl.cercos@upm.es>
 *
 *  AQUAgpusph is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  AQUAgpusph is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with AQUAgpusph.  If not, see <http://www.gnu.org/licenses/>.
 */

/** @file
 * @brief Math expression evaluator.
 * (See Aqua::Tokenizer for details)
 */

// #include <matheval.h>

#include "Tokenizer.hpp"
#include <cstdint>

using namespace std;

namespace Aqua {

std::mutex Tokenizer_muparser::mutex;

double
mod_operator(double v, double w)
{
#ifdef MAX
#undef MAX
#endif
#define MAX(a, b) ((a > b) ? a : b)
	return (int)v % MAX(1, (int)w);
#undef MAX
}

double
not_operator(double v)
{
	return v == 0.0;
}


Tokenizer_muparser::Tokenizer_muparser()
{
	struct lconv* lc;
	char* s;

	// Set the decimal-point character (which is depending on the locale)
	s = setlocale(LC_NUMERIC, NULL);
	if (strcmp(s, "C")) {
		std::ostringstream msg;
		msg << "\"" << s << "\" numeric locale found" << std::endl;
		LOG(L_INFO, msg.str());
		LOG0(L_DEBUG, "\tIt is replaced by \"C\"\n");
		setlocale(LC_NUMERIC, "C");
	}
	lc = localeconv();
	s = lc->decimal_point;
	if (strcmp(s, ".")) {
		std::ostringstream msg;
		msg << "\"" << s << "\" decimal point character found" << std::endl;
		LOG(L_WARNING, msg.str());
		LOG0(L_DEBUG, "\tIt is replaced by \".\"\n");
		lc->decimal_point = ".";
	}
	s = lc->thousands_sep;
	if (strcmp(s, "")) {
		std::ostringstream msg;
		msg << "\"" << s << "\" thousands separator character found"
		    << std::endl;
		LOG(L_WARNING, msg.str());
		LOG0(L_DEBUG, "\tIt is removed\n");
		lc->thousands_sep = "";
	}

	// Register a modulus operator
	p.DefineOprtChars("%");
	p.DefineOprt("%", mod_operator, mu::prINFIX);
	p.DefineInfixOprt("!", not_operator, 0);
	q.DefineOprtChars("%");
	q.DefineOprt("%", mod_operator, mu::prINFIX);
	q.DefineInfixOprt("!", not_operator, 0);
}

Tokenizer_muparser::~Tokenizer_muparser()
{
	p.ClearConst();
}

std::vector<std::string>
Tokenizer_muparser::exprVariables(const std::string eq)
{
	const std::lock_guard<std::mutex> lock(this->mutex);
	std::vector<std::string> vars;
	q.SetExpr(eq);
	auto variables = q.GetUsedVar();
	for (auto it = variables.begin(); it != variables.end(); ++it)
		vars.push_back(it->first);
	return vars;
}

bool
Tokenizer_muparser::isVariable(const std::string name)
{
	mu::valmap_type cmap = p.GetConst();
	if (cmap.size()) {
		mu::valmap_type::const_iterator item = cmap.begin();
		for (; item != cmap.end(); ++item) {
			if (!name.compare(item->first.c_str())) {
				return true;
			}
		}
	}
	return false;
}

template<>
double
Tokenizer_muparser::cast(const double& val)
{
	return val;
}

template<>
float
Tokenizer_muparser::cast(const double& val)
{
	return (float)val;
}

template<>
int32_t
Tokenizer_muparser::cast(const double& val)
{
	return Aqua::round<int32_t>(val);
}

template<>
int64_t
Tokenizer_muparser::cast(const double& val)
{
	return Aqua::round<int64_t>(val);
}

template<>
uint32_t
Tokenizer_muparser::cast(const double& val)
{
	return Aqua::round<uint32_t>(val);
}

template<>
uint64_t
Tokenizer_muparser::cast(const double& val)
{
	return Aqua::round<uint64_t>(val);
}

template<>
double
Tokenizer_muparser::cast(const int64_t& val)
{
	return (double)val;
}

template<>
float
Tokenizer_muparser::cast(const int64_t& val)
{
	return (float)val;
}

template<>
int32_t
Tokenizer_muparser::cast(const int64_t& val)
{
	return (int32_t)val;
}

template<>
int64_t
Tokenizer_muparser::cast(const int64_t& val)
{
	return (int64_t)val;
}

template<>
uint32_t
Tokenizer_muparser::cast(const int64_t& val)
{
	return (uint32_t)val;
}

template<>
uint64_t
Tokenizer_muparser::cast(const int64_t& val)
{
	return (uint64_t)val;
}

} // Aqua::
