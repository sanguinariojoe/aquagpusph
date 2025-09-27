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

#include <sstream>
#include <cstdint>
#include "Tokenizer.hpp"
#include "aquagpusph/AuxiliarMethods.hpp"

using namespace std;

namespace Aqua {

std::mutex Tokenizer_exprtk::mutex;

Tokenizer_exprtk::Tokenizer_exprtk()
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

	// Register the default constants
	using namespace exprtk::details::numeric::details;
	vars.add_constants();
	const typename number_type<tokenizer_t>::type num_type;
	static const tokenizer_t local_e = const_e_impl<tokenizer_t>(num_type);
	vars.add_constant("e", local_e);

}

Tokenizer_exprtk::~Tokenizer_exprtk()
{
	vars.clear();
}

std::vector<std::string>
Tokenizer_exprtk::exprVariables(const std::string eq)
{
	std::vector<std::string> vars_list;

	auto subeqs = split_formulae(eq);
	if (subeqs.size() > 1) {
		for (auto subeq : subeqs) {
			auto subvars = exprVariables(subeq);
			vars_list.insert(vars_list.end(), subvars.begin(), subvars.end());
		}
		return vars_list;
	}

	const std::lock_guard<std::mutex> lock(this->mutex);

	exprtk::parser<tokenizer_t> parser;
	exprtk::expression<tokenizer_t> expr;

	// Create a copy of the symbols list that we will extend until no more
	// undefined symbols are received
	exprtk::symbol_table<tokenizer_t> vars_tmp(vars);
	parser.dec().collect_variables() = true;
	while(true) {
		expr.register_symbol_table(vars_tmp);
		if (!parser.compile(eq, expr))
		{
			if (hasPrefix(parser.error(), "ERR234")) {
				auto var_name = replaceAllCopy(parser.error(),
				                               "ERR234 - Undefined symbol: '",
				                               "");
				replaceAll(var_name, "'", "");
				trim(var_name);
				vars_tmp.create_variable(var_name, 0.f);
				continue;
			}
			LOG(L_ERROR, std::string("Error parsing \"") + eq + "\":\n");
			LOG0(L_DEBUG, "\n");
			LOG0(L_DEBUG, parser.error() + "\n");
			LOG0(L_DEBUG, "\n");
			throw std::runtime_error("Invalid expression");
		}
		break;
	}

	std::deque<
		exprtk::parser<tokenizer_t>::dependent_entity_collector::symbol_t>
			symbol_list;
	parser.dec().symbols(symbol_list);
	for (auto symbol : symbol_list)
	{
		vars_list.push_back(symbol.first);
	}
	return vars_list;
}

} // Aqua::
