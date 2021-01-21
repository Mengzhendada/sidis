#include <iomanip>
#include <ios>
#include <iostream>
#include <random>
#include <stdexcept>
#include <string>

#include <sidis/sidis.hpp>
#include <sidis/extra/bounds.hpp>

using namespace sidis;
using namespace sidis::constant;
using namespace sidis::kin;
using namespace sidis::math;

// This program returns a random valid point in phase space.
int main(int argc, char** argv) {
	Real M_th = MASS_P + MASS_PI_0;
	Lepton beam = Lepton::E;
	Nucleus target = Nucleus::P;
	Hadron hadron = Hadron::PI_P;

	// Read input parameters from command line.
	Real beam_energy;
	bool radiative;
	try {
		if (argc != 3) {
			throw std::invalid_argument(
				"Unexpected number of command line arguments");
		}
		beam_energy = std::stold(argv[1]);
		std::string radiative_str = argv[2];
		if (radiative_str == "true"
				|| radiative_str == "on"
				|| radiative_str == "rad") {
			radiative = true;
		} else if (radiative_str == "false"
				|| radiative_str == "off"
				|| radiative_str == "nrad") {
			radiative = false;
		} else {
			throw std::out_of_range(
				"Must select radiative (rad) or "
				"non-radiative (nrad) phase space point");
		}
	} catch (std::exception const& e) {
		std::cerr << "Error: " << e.what() << std::endl;
		std::cout << "Usage: "
			<< "random_phase_space <E_b> <rad,nrad>"
			<< std::endl;
		return 1;
	}

	// Repeatedly choose a random point within phase space until we get one that
	// is kinematically valid.
	Initial initial_state(target, beam, beam_energy);
	Kinematics kin;
	KinematicsRad kin_rad;
	std::random_device rd;
	std::mt19937 rng(rd());
	std::uniform_real_distribution<Real> dist(0., 1.);
	do {
		Real x = x_bounds(initial_state).lerp(dist(rng));
		Real y = y_bounds(initial_state, x).lerp(dist(rng));
		Real z = z_bounds(initial_state, hadron, M_th, x, y).lerp(dist(rng));
		Real ph_t_sq = ph_t_sq_bounds(initial_state, hadron, M_th, x, y, z).lerp(dist(rng));
		Real phi_h = Bounds(-PI, PI).lerp(dist(rng));
		Real phi = Bounds(-PI, PI).lerp(dist(rng));
		PhaseSpace phase_space { x, y, z, ph_t_sq, phi_h, phi };
		kin = Kinematics(initial_state, phase_space, hadron, M_th);

		Real tau = tau_bounds(kin).lerp(dist(rng));
		Real phi_k = Bounds(-PI, PI).lerp(dist(rng));
		Real R = R_bounds(kin, tau, phi_k).lerp(dist(rng));
		kin_rad = KinematicsRad(kin, tau, phi_k, R);
	} while ((!radiative && !valid(kin)) || (radiative && !valid(kin_rad)));

	std::cout << std::scientific << std::setprecision(16);

	std::cout << "x     = " << kin.x << std::endl;
	std::cout << "y     = " << kin.y << std::endl;
	std::cout << "z     = " << kin.z << std::endl;
	std::cout << "ph_t² = " << kin.ph_t_sq << std::endl;
	std::cout << "φ_h   = " << kin.phi_h << std::endl;
	std::cout << "φ     = " << kin.phi << std::endl;
	if (radiative) {
		std::cout << "τ     = " << kin_rad.tau << std::endl;
		std::cout << "φ_k   = " << kin_rad.phi_k << std::endl;
		std::cout << "R     = " << kin_rad.R << std::endl;
	}

	return 0;
}
