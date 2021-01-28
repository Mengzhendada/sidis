#include "sidis/extra/exception.hpp"

using namespace sidis;
using namespace constant;

MassThresholdOutOfRange::MassThresholdOutOfRange(
	Real Mth,
	Real Mth_min,
	Real Mth_max) :
	_what(
		"Mass threshold " + std::to_string(Mth)
		+ " GeV is not within valid kinematic range " + std::to_string(Mth_min)
		+ " GeV to " + std::to_string(Mth_max) + " GeV"),
	Mth(Mth),
	Mth_min(Mth_min),
	Mth_max(Mth_max) { }

ComEnergyOutOfRange::ComEnergyOutOfRange(
	Real E_b,
	Real E_b_min) :
	_what(
		"Beam energy " + std::to_string(E_b)
		+ " GeV is not above threshold " + std::to_string(E_b_min)
		+ " GeV"),
	E_b(E_b),
	E_b_min(E_b_min) { }

TargetMismatch::TargetMismatch(
	Nucleus target,
	Nucleus target_expected) :
	_what(
		std::string("Target nucleus '") + name(target)
		+ "' is not the same as expected nucleus '" + name(target_expected)
		+ "'"),
	target(target),
	target_expected(target_expected) { }

FlavorOutOfRange::FlavorOutOfRange(unsigned flavor) :
	_what(
		"Flavor index " + std::to_string(flavor)
		+ " is out of range for TMD"),
	flavor(flavor) { }

HadronOutOfRange::HadronOutOfRange(Hadron hadron) :
	_what(
		std::string("Hadron '") + name(hadron)
		+ "' is not supported by FF"),
	hadron(hadron) { }

DataFileNotFound::DataFileNotFound(std::string file_name) :
	_what(
		"Could not open data file '" + file_name + "'"),
	file_name(file_name) { }

DataFileParseError::DataFileParseError(std::string file_name) :
	_what(
		"Could not parse data file '" + file_name + "'"),
	file_name(file_name) { }

