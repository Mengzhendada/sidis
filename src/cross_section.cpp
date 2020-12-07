#include "sidis/cross_section.hpp"

#include <cmath>
#include <limits>
#include <string>
#include <stdexcept>

#include <cubature.hpp>

#include "sidis/constant.hpp"
#include "sidis/frame.hpp"
#include "sidis/hadronic_coeff.hpp"
#include "sidis/kinematics.hpp"
#include "sidis/leptonic_coeff.hpp"
#include "sidis/structure_function.hpp"
#include "sidis/extra/bounds.hpp"
#include "sidis/extra/math.hpp"
#include "sidis/extra/integrate.hpp"
#include "sidis/extra/transform.hpp"
#include "sidis/extra/vector.hpp"

using namespace sidis;
using namespace sidis::xs;
using namespace sidis::constant;
using namespace sidis::had;
using namespace sidis::integ;
using namespace sidis::kin;
using namespace sidis::lep;
using namespace sidis::math;
using namespace sidis::sf;

Real const xs::SMALL_R_REL = std::cbrt(2. * std::numeric_limits<Real>::epsilon());

namespace {

Real const DELTA_R_REL = std::sqrt(2. * std::numeric_limits<Real>::epsilon());

// TODO: Move this into the header.
struct UnexpectedNucleusException : public std::runtime_error {
	Nucleus expected;
	Nucleus provided;
	UnexpectedNucleusException(Nucleus expected, Nucleus provided) :
		std::runtime_error(
			std::string("Target nucleus is ") + name(expected)
			+ " but provided structure functions are for " + name(provided)),
		expected(expected),
		provided(provided) { }
};

Real delta_vert_rad_0(Kinematics kin) {
	// Equation [1.3].
	Real Q_m_sq = kin.Q_sq + 2.*sq(kin.m);
	Real S_prime = kin.S - kin.Q_sq - kin.V_1;
	Real X_prime = kin.X + kin.Q_sq - kin.V_2;
	Real lambda_m = kin.Q_sq*(kin.Q_sq + 4.*sq(kin.m));
	Real lambda_S_prime = sq(S_prime) - 4.*sq(kin.m)*kin.mx_sq;
	Real lambda_X_prime = sq(X_prime) - 4.*sq(kin.m)*kin.mx_sq;
	Real lambda_m_sqrt = std::sqrt(lambda_m);
	Real lambda_S_prime_sqrt = std::sqrt(lambda_S_prime);
	Real lambda_X_prime_sqrt = std::sqrt(lambda_X_prime);

	// Differences of the form `√λ/|S| - 1`.
	Real diff_m = sqrt1p_1m((4.*sq(kin.m))/kin.Q_sq);
	Real diff_S_prime = sqrt1p_1m(-(4.*sq(kin.m)*kin.mx_sq)/sq(S_prime));
	Real diff_X_prime = sqrt1p_1m(-(4.*sq(kin.m)*kin.mx_sq)/sq(X_prime));
	Real sum_m = 2. + diff_m;
	Real sum_S_prime = 2. + diff_S_prime;
	Real sum_X_prime = 2. + diff_X_prime;

	// Equation [1.C10].
	Real L_m = 1./lambda_m_sqrt*std::log(sum_m/diff_m);
	Real L_S_prime = 1./lambda_S_prime_sqrt*std::log(-sum_S_prime/diff_S_prime);
	Real L_X_prime = 1./lambda_X_prime_sqrt*std::log(-sum_X_prime/diff_X_prime);

	// Equation [1.40].
	Real rho = 1./lambda_m_sqrt*(
		(Q_m_sq + lambda_m_sqrt)*S_prime
		- 2.*sq(kin.m)*X_prime);
	Real S_phi = Q_m_sq/lambda_m_sqrt*(
		lambda_S_prime*L_S_prime/4. - lambda_X_prime*L_X_prime/4.
		+ dilog(1. - 1./(sum_S_prime*S_prime)*rho)
		+ dilog(1. - (sum_S_prime*S_prime)/(4.*sq(kin.m)*kin.mx_sq)*rho)
		- dilog(1. - (sum_X_prime*X_prime)/(kin.mx_sq*sq(sum_m))*rho)
		- dilog(1. - (4.*sq(kin.m))/(sum_X_prime*X_prime*sq(sum_m))*rho));

	// Equation [1.52].
	Real delta = 0.5*S_prime*L_S_prime + 0.5*X_prime*L_X_prime + S_phi - 2.
		+ (1.5*kin.Q_sq + 4.*sq(kin.m))*L_m
		- Q_m_sq/lambda_m_sqrt*(
			0.5*lambda_m*sq(L_m)
			+ 2.*dilog((2.*lambda_m_sqrt)/(kin.Q_sq + lambda_m_sqrt))
			- 0.5*sq(PI));
	return delta;
}

}

Real xs::born(Real lambda_e, Vec3 eta, Kinematics kin, Model const& model) {
	if (model.target != kin.target) {
		throw UnexpectedNucleusException(kin.target, model.target);
	}
	Born b(kin);
	SfXX sf = model.sf(kin.hadron, kin.x, kin.z, kin.Q_sq, kin.ph_t_sq);
	LepBornXX lep(kin);
	HadXX had(kin, sf);

	Real uu = born_uu_base(b, lep, had);
	Vec3 up(
		born_ut1_base(b, lep, had),
		born_ut2_base(b, lep, had),
		born_ul_base(b, lep, had));
	Real lu = born_lu_base(b, lep, had);
	Vec3 lp(
		born_lt1_base(b, lep, had),
		born_lt2_base(b, lep, had),
		born_ll_base(b, lep, had));
	return uu + dot(eta, up) + lambda_e * (lu + dot(eta, lp));
}

Real xs::amm(Real lambda_e, Vec3 eta, Kinematics kin, Model const& model) {
	if (model.target != kin.target) {
		throw UnexpectedNucleusException(kin.target, model.target);
	}
	Amm b(kin);
	SfXX sf = model.sf(kin.hadron, kin.x, kin.z, kin.Q_sq, kin.ph_t_sq);
	LepAmmXX lep(kin);
	HadXX had(kin, sf);

	Real uu = amm_uu_base(b, lep, had);
	Vec3 up(
		amm_ut1_base(b, lep, had),
		amm_ut2_base(b, lep, had),
		amm_ul_base(b, lep, had));
	Real lu = amm_lu_base(b, lep, had);
	Vec3 lp(
		amm_lt1_base(b, lep, had),
		amm_lt2_base(b, lep, had),
		amm_ll_base(b, lep, had));
	return uu + dot(eta, up) + lambda_e * (lu + dot(eta, lp));
}

Real xs::nrad_ir(Real lambda_e, Vec3 eta, Kinematics kin, Model const& model, Real k0_cut) {
	if (model.target != kin.target) {
		throw UnexpectedNucleusException(kin.target, model.target);
	}
	NRadIR b(kin, k0_cut);
	SfXX sf = model.sf(kin.hadron, kin.x, kin.z, kin.Q_sq, kin.ph_t_sq);
	LepBornXX lep_born(kin);
	LepAmmXX lep_amm(kin);
	HadXX had(kin, sf);

	Real uu = nrad_ir_uu_base(b, lep_born, lep_amm, had);
	Vec3 up(
		nrad_ir_ut1_base(b, lep_born, lep_amm, had),
		nrad_ir_ut2_base(b, lep_born, lep_amm, had),
		nrad_ir_ul_base(b, lep_born, lep_amm, had));
	Real lu = nrad_ir_lu_base(b, lep_born, lep_amm, had);
	Vec3 lp(
		nrad_ir_lt1_base(b, lep_born, lep_amm, had),
		nrad_ir_lt2_base(b, lep_born, lep_amm, had),
		nrad_ir_ll_base(b, lep_born, lep_amm, had));
	return uu + dot(eta, up) + lambda_e * (lu + dot(eta, lp));
}

Real xs::nrad(Real lambda_e, Vec3 eta, Kinematics kin, Model const& model, Real k0_cut) {
	// The soft part of the radiative cross-section (below `k0_cut`) is bundled
	// into the return value here.
	Real xs_nrad_ir = nrad_ir(lambda_e, eta, kin, model, k0_cut);
	Real xs_rad_f = rad_f_integ(lambda_e, eta, kin, model, k0_cut);
	return xs_nrad_ir + xs_rad_f;
}

Real xs::rad(Real lambda_e, Vec3 eta, KinematicsRad kin, Model const& model) {
	if (model.target != kin.target) {
		throw UnexpectedNucleusException(kin.target, model.target);
	}
	Rad b(kin);
	Transform3 shift_rot = frame::hadron_from_shift(kin);
	Sf shift_sf = model.sf(
		kin.hadron,
		kin.shift_x, kin.shift_z, kin.shift_Q_sq, kin.shift_ph_t_sq);
	LepRadXX lep(kin);
	HadXX shift_had(kin.project_shift(), shift_sf);
	Real uu = rad_uu_base(b, lep, shift_had);
	Vec3 up = rad_up_base(b, lep, shift_had, shift_rot);
	Real lu = rad_lu_base(b, lep, shift_had);
	Vec3 lp = rad_lp_base(b, lep, shift_had, shift_rot);
	return uu + dot(eta, up) + lambda_e * (lu + dot(eta, lp));
}

Real xs::rad_f(Real lambda_e, Vec3 eta, KinematicsRad kin, Model const& model) {
	if (model.target != kin.target) {
		throw UnexpectedNucleusException(kin.target, model.target);
	}
	Rad b(kin);
	Transform3 shift_rot = frame::hadron_from_shift(kin);
	Sf sf = model.sf(
		kin.hadron,
		kin.x, kin.z, kin.Q_sq, kin.ph_t_sq);
	Sf shift_sf = model.sf(
		kin.hadron,
		kin.shift_x, kin.shift_z, kin.shift_Q_sq, kin.shift_ph_t_sq);
	LepRadXX lep(kin);
	HadXX shift_had(kin.project_shift(), shift_sf);
	HadXX had(kin.project(), sf);

	// Compute hard cross-section.
	Real uu_h = rad_f_hard_uu_base(b, lep, shift_had);
	Vec3 up_h = rad_f_hard_up_base(b, lep, shift_had, shift_rot);
	Real lu_h = rad_f_hard_lu_base(b, lep, shift_had);
	Vec3 lp_h = rad_f_hard_lp_base(b, lep, shift_had, shift_rot);
	Real xs_h = uu_h + dot(eta, up_h) + lambda_e * (lu_h + dot(eta, lp_h));

	// Compute soft cross-section.
	Real xs_s;
	if (std::abs(kin.R) < std::abs(kin.R_max) * SMALL_R_REL) {
		xs_s = rad_f_soft_0(lambda_e, eta, kin, model);
	} else {
		Real uu_s = rad_f_soft_uu_base(b, lep, had, shift_had);
		Vec3 up_s = rad_f_soft_up_base(b, lep, had, shift_had, shift_rot);
		Real lu_s = rad_f_soft_lu_base(b, lep, had, shift_had);
		Vec3 lp_s = rad_f_soft_lp_base(b, lep, had, shift_had, shift_rot);
		xs_s = uu_s + dot(eta, up_s) + lambda_e * (lu_s + dot(eta, lp_s));
	}
	return xs_s + xs_h;
}

Real xs::rad_f_integ(Real lambda_e, math::Vec3 eta, kin::Kinematics kin, sf::Model const& model, Real k0_cut) {
	if (model.target != kin.target) {
		throw UnexpectedNucleusException(kin.target, model.target);
	}
	Sf sf = model.sf(
		kin.hadron,
		kin.x, kin.z, kin.Q_sq, kin.ph_t_sq);
	HadXX had(kin, sf);
	cubature::EstErr<Real> xs_integ = cubature::cubature<3>(
		[&](cubature::Point<3, Real> x) {
			Bounds tau_b = tau_bounds(kin);
			Real tau = tau_b.lerp(x[0]);

			Bounds phi_k_b(0., 2. * PI);
			Real phi_k = phi_k_b.lerp(x[1]);

			Bounds R_b = R_bounds_soft(kin, tau, phi_k, k0_cut);
			Real R = R_b.lerp(x[2]);

			Real jacobian = tau_b.size() * phi_k_b.size() * R_b.size();

			KinematicsRad kin_rad(kin, tau, phi_k, R);
			Rad b(kin_rad);
			Transform3 shift_rot = frame::hadron_from_shift(kin_rad);
			Sf shift_sf = model.sf(
				kin.hadron,
				kin_rad.shift_x, kin_rad.shift_z, kin_rad.shift_Q_sq, kin_rad.shift_ph_t_sq);
			LepRadXX lep(kin_rad);
			HadXX shift_had(kin_rad.project_shift(), shift_sf);

			// Compute hard cross-section.
			Real uu_h = rad_f_hard_uu_base(b, lep, shift_had);
			Vec3 up_h = rad_f_hard_up_base(b, lep, shift_had, shift_rot);
			Real lu_h = rad_f_hard_lu_base(b, lep, shift_had);
			Vec3 lp_h = rad_f_hard_lp_base(b, lep, shift_had, shift_rot);
			Real xs_h = uu_h + dot(eta, up_h) + lambda_e * (lu_h + dot(eta, lp_h));

			// Compute soft cross-section.
			Real xs_s;
			if (std::abs(kin_rad.R) < std::abs(kin_rad.R_max) * SMALL_R_REL) {
				xs_s = rad_f_soft_0(lambda_e, eta, kin_rad, model);
			} else {
				Real uu_s = rad_f_soft_uu_base(b, lep, had, shift_had);
				Vec3 up_s = rad_f_soft_up_base(b, lep, had, shift_had, shift_rot);
				Real lu_s = rad_f_soft_lu_base(b, lep, had, shift_had);
				Vec3 lp_s = rad_f_soft_lp_base(b, lep, had, shift_had, shift_rot);
				xs_s = uu_s + dot(eta, up_s) + lambda_e * (lu_s + dot(eta, lp_s));
			}
			Real xs = xs_s + xs_h;

			if (std::isnan(xs)) {
				// If the result is `NaN`, it most likely means we went out of
				// the allowed region for the structure function grids. In that
				// case, just return zero.
				// TODO: Handle this case in a more correct way.
				return 0.;
			} else {
				return jacobian * xs;
			}
		},
		cubature::Point<3, Real>{ 0., 0., 0. },
		cubature::Point<3, Real>{ 1., 1., 1. },
		100000, 0., 1e-6);
	return xs_integ.val;
}

Real xs::rad_integ(Real lambda_e, math::Vec3 eta, kin::Kinematics kin, sf::Model const& model, Real k0_cut) {
	if (model.target != kin.target) {
		throw UnexpectedNucleusException(kin.target, model.target);
	}
	cubature::EstErr<Real> xs_integ = cubature::cubature<3>(
		[&](cubature::Point<3, Real> x) {
			Bounds tau_b = tau_bounds(kin);
			Real tau = tau_b.lerp(x[0]);

			Bounds phi_k_b(0., 2. * PI);
			Real phi_k = phi_k_b.lerp(x[1]);

			Bounds R_b = R_bounds_hard(kin, tau, phi_k, k0_cut);
			Real R = R_b.lerp(x[2]);

			Real jacobian = tau_b.size() * phi_k_b.size() * R_b.size();

			KinematicsRad kin_rad(kin, tau, phi_k, R);
			Rad b(kin_rad);
			Transform3 shift_rot = frame::hadron_from_shift(kin_rad);
			Sf shift_sf = model.sf(
				kin.hadron,
				kin_rad.shift_x, kin_rad.shift_z, kin_rad.shift_Q_sq, kin_rad.shift_ph_t_sq);
			LepRadXX lep(kin_rad);
			HadXX shift_had(kin_rad.project_shift(), shift_sf);

			Real uu = rad_uu_base(b, lep, shift_had);
			Vec3 up = rad_up_base(b, lep, shift_had, shift_rot);
			Real lu = rad_lu_base(b, lep, shift_had);
			Vec3 lp = rad_lp_base(b, lep, shift_had, shift_rot);
			Real xs = uu + dot(eta, up) + lambda_e * (lu + dot(eta, lp));
			if (std::isnan(xs)) {
				// TODO: Handle this case more correctly.
				return 0.;
			} else {
				return jacobian * xs;
			}
		},
		cubature::Point<3, Real>{ 0., 0., 0. },
		cubature::Point<3, Real>{ 1., 1., 1. },
		100000, 0., 1e-6);
	return xs_integ.val;
}

Real xs::rad_f_hard(Real lambda_e, Vec3 eta, KinematicsRad kin, Model const& model) {
	if (model.target != kin.target) {
		throw UnexpectedNucleusException(kin.target, model.target);
	}
	Rad b(kin);
	Transform3 shift_rot = frame::hadron_from_shift(kin);
	Sf shift_sf = model.sf(
		kin.hadron,
		kin.shift_x, kin.shift_z, kin.shift_Q_sq, kin.shift_ph_t_sq);
	LepRadXX lep(kin);
	HadXX shift_had(kin.project_shift(), shift_sf);

	Real uu = rad_f_hard_uu_base(b, lep, shift_had);
	Vec3 up = rad_f_hard_up_base(b, lep, shift_had, shift_rot);
	Real lu = rad_f_hard_lu_base(b, lep, shift_had);
	Vec3 lp = rad_f_hard_lp_base(b, lep, shift_had, shift_rot);
	return uu + dot(eta, up) + lambda_e * (lu + dot(eta, lp));
}

Real xs::rad_f_soft(Real lambda_e, Vec3 eta, KinematicsRad kin, Model const& model) {
	if (model.target != kin.target) {
		throw UnexpectedNucleusException(kin.target, model.target);
	}
	if (std::abs(kin.R) < std::abs(kin.R_max) * SMALL_R_REL) {
		return rad_f_soft_0(lambda_e, eta, kin, model);
	} else {
		Rad b(kin);
		Transform3 shift_rot = frame::hadron_from_shift(kin);
		Sf sf = model.sf(
			kin.hadron,
			kin.x, kin.z, kin.Q_sq, kin.ph_t_sq);
		Sf shift_sf = model.sf(
			kin.hadron,
			kin.shift_x, kin.shift_z, kin.shift_Q_sq, kin.shift_ph_t_sq);
		LepRadXX lep(kin);
		HadXX had(kin.project(), sf);
		HadXX shift_had(kin.project_shift(), shift_sf);

		Real uu = rad_f_soft_uu_base(b, lep, had, shift_had);
		Vec3 up = rad_f_soft_up_base(b, lep, had, shift_had, shift_rot);
		Real lu = rad_f_soft_lu_base(b, lep, had, shift_had);
		Vec3 lp = rad_f_soft_lp_base(b, lep, had, shift_had, shift_rot);
		return uu + dot(eta, up) + lambda_e * (lu + dot(eta, lp));
	}
}

Real xs::rad_f_soft_0(Real lambda_e, Vec3 eta, KinematicsRad kin, Model const& model) {
	if (model.target != kin.target) {
		throw UnexpectedNucleusException(kin.target, model.target);
	}
	// For small `R`, the evaluation of the soft contribution to the infrared-
	// divergent-free part of the radiative cross-section becomes inaccurate due
	// to catastrophic cancellation. At the same time, this is the region that
	// contributes the most to the cross-section, so it's important to get
	// right. We use linear extrapolation to extend the radiative cross-section
	// to this region. The first and second derivatives of the structure
	// functions are used for this.
	Real delta_R = kin.R_max * DELTA_R_REL;
	Real R_rel = kin.R / delta_R;
	Real R_1 = delta_R;
	Real R_2 = 2. * delta_R;
	Rad b(kin);
	KinematicsRad kin_1(kin.project(), kin.tau, kin.phi_k, R_1);
	KinematicsRad kin_2(kin.project(), kin.tau, kin.phi_k, R_2);
	Transform3 shift_rot_1 = frame::hadron_from_shift(kin_1);
	Transform3 shift_rot_2 = frame::hadron_from_shift(kin_2);
	Sf sf_0 = model.sf(
		kin_1.hadron,
		kin_1.x, kin_1.z, kin_1.Q_sq, kin_1.ph_t_sq);
	Sf sf_1 = model.sf(
		kin_1.hadron,
		kin_1.shift_x, kin_1.shift_z, kin_1.shift_Q_sq, kin_1.shift_ph_t_sq);
	Sf sf_2 = model.sf(
		kin_2.hadron,
		kin_2.shift_x, kin_2.shift_z, kin_2.shift_Q_sq, kin_2.shift_ph_t_sq);
	LepRadXX lep(kin);
	HadXX had_0(kin.project(), sf_0);
	HadXX had_1(kin_1.project_shift(), sf_1);
	HadXX had_2(kin_2.project_shift(), sf_2);
	Real uu_0 = rad_f_soft_uu_base_R0(b, lep, had_0);
	Vec3 up_0(
		rad_f_soft_ul_base_R0(b, lep, had_0),
		rad_f_soft_ut1_base_R0(b, lep, had_0),
		rad_f_soft_ut2_base_R0(b, lep, had_0));
	Real lu_0 = rad_f_soft_lu_base_R0(b, lep, had_0);
	Vec3 lp_0(
		rad_f_soft_ll_base_R0(b, lep, had_0),
		rad_f_soft_lt1_base_R0(b, lep, had_0),
		rad_f_soft_lt2_base_R0(b, lep, had_0));
	Real uu_1 = rad_f_soft_uu_base_R(b, lep, had_1);
	Vec3 up_1 = rad_f_soft_up_base_R(b, lep, had_1, shift_rot_1);
	Real lu_1 = rad_f_soft_lu_base_R(b, lep, had_1);
	Vec3 lp_1 = rad_f_soft_lp_base_R(b, lep, had_1, shift_rot_1);
	Real uu_2 = rad_f_soft_uu_base_R(b, lep, had_2);
	Vec3 up_2 = rad_f_soft_up_base_R(b, lep, had_2, shift_rot_2);
	Real lu_2 = rad_f_soft_lu_base_R(b, lep, had_2);
	Vec3 lp_2 = rad_f_soft_lp_base_R(b, lep, had_2, shift_rot_2);
	Real xs_R0 = uu_0 + dot(eta, up_0) + lambda_e * (lu_0 + dot(eta, lp_0));
	Real xs_R1 = uu_1 + dot(eta, up_1) + lambda_e * (lu_1 + dot(eta, lp_1));
	Real xs_R2 = uu_2 + dot(eta, up_2) + lambda_e * (lu_2 + dot(eta, lp_2));

	// Compute first and second derivatives, and use the result as part of a
	// linear extrapolation for the cross-section.
	return 1. / delta_R * (
		(xs_R1 - xs_R0)
		+ 0.5 * R_rel * (xs_R2 - 2. * xs_R1 + xs_R0));
}

// Radiative corrections to Born cross-section.
Real xs::delta_vert_rad_ir(Kinematics kin, Real k0_cut) {
	// Paragraph following equation [1.C17].
	Real k0_max = (kin.mx_sq - sq(kin.M_th))/(2.*kin.mx);
	if (!(k0_cut > 0.)) {
		return -INFINITY;
	}
	Real Q_m_sq = kin.Q_sq + 2.*sq(kin.m);
	Real diff_m = sqrt1p_1m((4.*sq(kin.m))/kin.Q_sq);
	Real sum_m = 2. + diff_m;
	Real lambda_m = kin.Q_sq*(kin.Q_sq + 4.*sq(kin.m));
	Real lambda_m_sqrt = std::sqrt(lambda_m);
	Real L_m = 1./lambda_m_sqrt*std::log(sum_m/diff_m);
	Real delta_0 = delta_vert_rad_0(kin);
	// This comes from subtracting `delta_H` (equation [1.38]) from `delta_VR`
	// (equation [1.52]).
	Real delta_shift = 2.*(Q_m_sq*L_m - 1.)*std::log(
		k0_cut < k0_max ?
		(2.*k0_cut)/kin.m :
		(kin.mx_sq - sq(kin.M_th))/(kin.m*kin.mx));
	return delta_0 + delta_shift;
}
Real xs::delta_rad_ir_hard(Kinematics kin, Real k0_cut) {
	Real k0_max = (kin.mx_sq - sq(kin.M_th))/(2.*kin.mx);
	if (!(k0_cut > 0.)) {
		return INFINITY;
	} else if (!(k0_cut < k0_max)) {
		return 0.;
	}
	Real Q_m_sq = kin.Q_sq + 2.*sq(kin.m);
	Real diff_m = sqrt1p_1m((4.*sq(kin.m))/kin.Q_sq);
	Real sum_m = 2. + diff_m;
	Real lambda_m = kin.Q_sq*(kin.Q_sq + 4.*sq(kin.m));
	Real lambda_m_sqrt = std::sqrt(lambda_m);
	Real L_m = 1./lambda_m_sqrt*std::log(sum_m/diff_m);
	// Equation [1.38].
	Real delta = 2.*(Q_m_sq*L_m - 1.)*std::log(
		(kin.mx_sq - sq(kin.M_th))/(2.*k0_cut*kin.mx));
	return delta;
}

Real xs::delta_vac_lep(Kinematics kin) {
	// Equation [1.50].
	Real ms[3] = { MASS_E, MASS_MU, MASS_TAU };
	Real delta = 0.;
	for (unsigned idx = 0; idx < 3; ++idx) {
		Real m = ms[idx];
		Real lambda_sqrt = std::sqrt(kin.Q_sq*(kin.Q_sq + 4.*sq(m)));
		Real diff_m = sqrt1p_1m((4.*sq(m))/kin.Q_sq);
		Real sum_m = 2. + diff_m;
		Real L_m = 1/lambda_sqrt*std::log(sum_m/diff_m);
		delta += 2./3.*(kin.Q_sq
			+ 2.*sq(m))*L_m
			- 10./9.
			+ (8.*sq(m))/(3.*kin.Q_sq)*(1. - 2.*sq(m)*L_m);
	}
	return delta;
}

Real xs::delta_vac_had(Kinematics kin) {
	if (kin.Q_sq < 1.) {
		return -(2.*PI)/ALPHA*(-1.345e-9 - 2.302e-3*std::log(1. + 4.091*kin.Q_sq));
	} else if (kin.Q_sq < 64.) {
		return -(2.*PI)/ALPHA*(-1.512e-3 - 2.822e-3*std::log(1. + 1.218*kin.Q_sq));
	} else {
		return -(2.*PI)/ALPHA*(-1.1344-3 - 3.0680-3*std::log(1. + 0.99992*kin.Q_sq));
	}
}

// Born base functions.
Born::Born(Kinematics kin) :
	// Equation [1.15]. The `Q^4` factor has been absorbed into `C_1`.
	coeff((sq(ALPHA)*kin.S*sq(kin.S_x))/(8.*kin.M*kin.ph_l*kin.lambda_S)) { }

Real xs::born_uu_base(Born b, LepBornUU lep, HadUU had) {
	return b.coeff*(
		lep.theta_1*had.H_10
		+ lep.theta_2*had.H_20
		+ lep.theta_3*had.H_30
		+ lep.theta_4*had.H_40);
}
Real xs::born_ul_base(Born b, LepBornUP lep, HadUL had) {
	return b.coeff*(lep.theta_6*had.H_63 + lep.theta_8*had.H_83);
}
Real xs::born_ut1_base(Born b, LepBornUP lep, HadUT1 had) {
	return b.coeff*(lep.theta_6*had.H_61 + lep.theta_8*had.H_81);
}
Real xs::born_ut2_base(Born b, LepBornUU lep, HadUT2 had) {
	return b.coeff*(
		lep.theta_1*had.H_12
		+ lep.theta_2*had.H_22
		+ lep.theta_3*had.H_32
		+ lep.theta_4*had.H_42);
}
Real xs::born_lu_base(Born b, LepBornLU lep, HadLU had) {
	return b.coeff*lep.theta_5*had.H_50;
}
Real xs::born_ll_base(Born b, LepBornLP lep, HadLL had) {
	return b.coeff*(lep.theta_7*had.H_73 + lep.theta_9*had.H_93);
}
Real xs::born_lt1_base(Born b, LepBornLP lep, HadLT1 had) {
	return b.coeff*(lep.theta_7*had.H_71 + lep.theta_9*had.H_91);
}
Real xs::born_lt2_base(Born b, LepBornLU lep, HadLT2 had) {
	return b.coeff*lep.theta_5*had.H_52;
}

// AMM base functions.
Amm::Amm(Kinematics kin) {
	// Equation [1.53]. The `Q^4` factor has been absorbed into `C_1`.
	Real lambda_m = kin.Q_sq*(kin.Q_sq + 4.*sq(kin.m));
	Real lambda_m_sqrt = std::sqrt(lambda_m);
	Real diff_m = sqrt1p_1m((4.*sq(kin.m))/kin.Q_sq);
	Real sum_m = 2. + diff_m;
	Real L_m = 1./lambda_m_sqrt*std::log(sum_m/diff_m);
	coeff = L_m*kin.Q_sq*(std::pow(ALPHA, 3)*sq(kin.m)*kin.S*sq(kin.S_x))
		/(16.*PI*kin.M*kin.ph_l*kin.lambda_S);
}

Real xs::amm_uu_base(Amm b, LepAmmUU lep, HadUU had) {
	return b.coeff*(
		lep.theta_1*had.H_10
		+ lep.theta_2*had.H_20
		+ lep.theta_3*had.H_30
		+ lep.theta_4*had.H_40);
}
Real xs::amm_ul_base(Amm b, LepAmmUP lep, HadUL had) {
	return b.coeff*(lep.theta_6*had.H_63 + lep.theta_8*had.H_83);
}
Real xs::amm_ut1_base(Amm b, LepAmmUP lep, HadUT1 had) {
	return b.coeff*(lep.theta_6*had.H_61 + lep.theta_8*had.H_81);
}
Real xs::amm_ut2_base(Amm b, LepAmmUU lep, HadUT2 had) {
	return b.coeff*(
		lep.theta_1*had.H_12
		+ lep.theta_2*had.H_22
		+ lep.theta_3*had.H_32
		+ lep.theta_4*had.H_42);
}
Real xs::amm_lu_base(Amm b, LepAmmLU lep, HadLU had) {
	return b.coeff*lep.theta_5*had.H_50;
}
Real xs::amm_ll_base(Amm b, LepAmmLP lep, HadLL had) {
	return b.coeff*(lep.theta_7*had.H_73 + lep.theta_9*had.H_93);
}
Real xs::amm_lt1_base(Amm b, LepAmmLP lep, HadLT1 had) {
	return b.coeff*(lep.theta_7*had.H_71 + lep.theta_9*had.H_91);
}
Real xs::amm_lt2_base(Amm b, LepAmmLU lep, HadLT2 had) {
	return b.coeff*lep.theta_5*had.H_52;
}

// Non-radiative infrared-divergence-free base functions.
NRadIR::NRadIR(Kinematics kin, Real k0_cut) {
	Born born(kin);
	Amm amm(kin);
	Real born_factor = 1. + ALPHA/PI*(
		delta_vert_rad_ir(kin, k0_cut)
		+ delta_vac_lep(kin)
		+ delta_vac_had(kin));
	coeff_born = born_factor*born.coeff;
	coeff_amm = amm.coeff;
}
Real xs::nrad_ir_uu_base(NRadIR b, LepBornUU lep_born, LepAmmUU lep_amm, HadUU had) {
	return
		(b.coeff_born*lep_born.theta_1 + b.coeff_amm*lep_amm.theta_1)*had.H_10
		+ (b.coeff_born*lep_born.theta_2 + b.coeff_amm*lep_amm.theta_2)*had.H_20
		+ (b.coeff_born*lep_born.theta_3 + b.coeff_amm*lep_amm.theta_3)*had.H_30
		+ (b.coeff_born*lep_born.theta_4 + b.coeff_amm*lep_amm.theta_4)*had.H_40;
}
Real xs::nrad_ir_ul_base(NRadIR b, LepBornUP lep_born, LepAmmUP lep_amm, HadUL had) {
	return
		(b.coeff_born*lep_born.theta_6 + b.coeff_amm*lep_amm.theta_6)*had.H_63
		+ (b.coeff_born*lep_born.theta_8 + b.coeff_amm*lep_amm.theta_8)*had.H_83;
}
Real xs::nrad_ir_ut1_base(NRadIR b, LepBornUP lep_born, LepAmmUP lep_amm, HadUT1 had) {
	return
		(b.coeff_born*lep_born.theta_6 + b.coeff_amm*lep_amm.theta_6)*had.H_61
		+ (b.coeff_born*lep_born.theta_8 + b.coeff_amm*lep_amm.theta_8)*had.H_81;
}
Real xs::nrad_ir_ut2_base(NRadIR b, LepBornUU lep_born, LepAmmUU lep_amm, HadUT2 had) {
	return
		(b.coeff_born*lep_born.theta_1 + b.coeff_amm*lep_amm.theta_1)*had.H_12
		+ (b.coeff_born*lep_born.theta_2 + b.coeff_amm*lep_amm.theta_2)*had.H_22
		+ (b.coeff_born*lep_born.theta_3 + b.coeff_amm*lep_amm.theta_3)*had.H_32
		+ (b.coeff_born*lep_born.theta_4 + b.coeff_amm*lep_amm.theta_4)*had.H_42;
}
Real xs::nrad_ir_lu_base(NRadIR b, LepBornLU lep_born, LepAmmLU lep_amm, HadLU had) {
	return (b.coeff_born*lep_born.theta_5 + b.coeff_amm*lep_amm.theta_5)*had.H_50;
}
Real xs::nrad_ir_ll_base(NRadIR b, LepBornLP lep_born, LepAmmLP lep_amm, HadLL had) {
	return
		(b.coeff_born*lep_born.theta_7 + b.coeff_amm*lep_amm.theta_7)*had.H_73
		+ (b.coeff_born*lep_born.theta_9 + b.coeff_amm*lep_amm.theta_9)*had.H_93;
}
Real xs::nrad_ir_lt1_base(NRadIR b, LepBornLP lep_born, LepAmmLP lep_amm, HadLT1 had) {
	return
		(b.coeff_born*lep_born.theta_7 + b.coeff_amm*lep_amm.theta_7)*had.H_71
		+ (b.coeff_born*lep_born.theta_9 + b.coeff_amm*lep_amm.theta_9)*had.H_91;
}
Real xs::nrad_ir_lt2_base(NRadIR b, LepBornLU lep_born, LepAmmLU lep_amm, HadLT2 had) {
	return (b.coeff_born*lep_born.theta_5 + b.coeff_amm*lep_amm.theta_5)*had.H_52;
}

// Radiative base functions.
Rad::Rad(KinematicsRad kin) {
	// Equation [1.43].
	coeff = -(std::pow(ALPHA, 3)*kin.S*sq(kin.S_x))
		/(64.*sq(PI)*kin.M*kin.ph_l*kin.lambda_S*kin.lambda_Y_sqrt);
	R = kin.R;
}

Real xs::rad_uu_base(Rad b, LepRadUU lep, HadUU shift_had) {
	return b.coeff*(
		1./b.R*(
			lep.theta_011*shift_had.H_10
			+ lep.theta_021*shift_had.H_20
			+ lep.theta_031*shift_had.H_30
			+ lep.theta_041*shift_had.H_40)
		+ (
			lep.theta_012*shift_had.H_10
			+ lep.theta_022*shift_had.H_20
			+ lep.theta_032*shift_had.H_30
			+ lep.theta_042*shift_had.H_40)
		+ b.R*(
			lep.theta_013*shift_had.H_10
			+ lep.theta_023*shift_had.H_20
			+ lep.theta_033*shift_had.H_30
			+ lep.theta_043*shift_had.H_40));
}
Vec3 xs::rad_up_base(Rad b, LepRadUX lep, HadUP shift_had, Transform3 shift_rot) {
	return b.coeff*dot(shift_rot, Vec3(
		// UT1.
		1./b.R*(lep.theta_061*shift_had.H_61 + lep.theta_081*shift_had.H_81)
		+ (lep.theta_062*shift_had.H_61 + lep.theta_082*shift_had.H_81)
		+ b.R*(lep.theta_063*shift_had.H_61 + lep.theta_083*shift_had.H_81)
		+ b.R*b.R*(lep.theta_064*shift_had.H_61 + lep.theta_084*shift_had.H_81),
		// UT2.
		1./b.R*(
			lep.theta_011*shift_had.H_12
			+ lep.theta_021*shift_had.H_22
			+ lep.theta_031*shift_had.H_32
			+ lep.theta_041*shift_had.H_42)
		+ (
			lep.theta_012*shift_had.H_12
			+ lep.theta_022*shift_had.H_22
			+ lep.theta_032*shift_had.H_32
			+ lep.theta_042*shift_had.H_42)
		+ b.R*(
			lep.theta_013*shift_had.H_12
			+ lep.theta_023*shift_had.H_22
			+ lep.theta_033*shift_had.H_32
			+ lep.theta_043*shift_had.H_42),
		// UL.
		1./b.R*(lep.theta_061*shift_had.H_63 + lep.theta_081*shift_had.H_83)
		+ (lep.theta_062*shift_had.H_63 + lep.theta_082*shift_had.H_83)
		+ b.R*(lep.theta_063*shift_had.H_63 + lep.theta_083*shift_had.H_83)
		+ b.R*b.R*(lep.theta_064*shift_had.H_63 + lep.theta_084*shift_had.H_83)));
}
Real xs::rad_lu_base(Rad b, LepRadLU lep, HadLU shift_had) {
	return b.coeff*(
		1./b.R*(lep.theta_051 + lep.theta_151)*shift_had.H_50
		+ (lep.theta_052 + lep.theta_152)*shift_had.H_50
		+ b.R*(lep.theta_053 + lep.theta_153)*shift_had.H_50);
}
Vec3 xs::rad_lp_base(Rad b, LepRadLX lep, HadLP shift_had, Transform3 shift_rot) {
	return b.coeff*dot(shift_rot, Vec3(
		// LT1.
		1./b.R*(
			(lep.theta_071 + lep.theta_171)*shift_had.H_71
			+ (lep.theta_091 + lep.theta_191)*shift_had.H_91)
		+ (
			(lep.theta_072 + lep.theta_172)*shift_had.H_71
			+ (lep.theta_092 + lep.theta_192)*shift_had.H_91)
		+ b.R*(
			(lep.theta_073 + lep.theta_173)*shift_had.H_71
			+ (lep.theta_093 + lep.theta_193)*shift_had.H_91)
		+ b.R*b.R*(
			(lep.theta_074 + lep.theta_174)*shift_had.H_71
			+ (lep.theta_094 + lep.theta_194)*shift_had.H_91),
		// LT2.
		1./b.R*(lep.theta_051 + lep.theta_151)*shift_had.H_52
		+ (lep.theta_052 + lep.theta_152)*shift_had.H_52
		+ b.R*(lep.theta_053 + lep.theta_153)*shift_had.H_52,
		// LL.
		1./b.R*(
			(lep.theta_071 + lep.theta_171)*shift_had.H_73
			+ (lep.theta_091 + lep.theta_191)*shift_had.H_93)
		+ (
			(lep.theta_072 + lep.theta_172)*shift_had.H_73
			+ (lep.theta_092 + lep.theta_192)*shift_had.H_93)
		+ b.R*(
			(lep.theta_073 + lep.theta_173)*shift_had.H_73
			+ (lep.theta_093 + lep.theta_193)*shift_had.H_93)
		+ b.R*b.R*(
			(lep.theta_074 + lep.theta_174)*shift_had.H_73
			+ (lep.theta_094 + lep.theta_194)*shift_had.H_93)));
}

// "Hard" infrared-divergent-free radiative base functions.
Real xs::rad_f_hard_uu_base(Rad b, LepRadUU lep, HadUU shift_had) {
	return b.coeff*(
		(
			lep.theta_012*shift_had.H_10
			+ lep.theta_022*shift_had.H_20
			+ lep.theta_032*shift_had.H_30
			+ lep.theta_042*shift_had.H_40)
		+ b.R*(
			lep.theta_013*shift_had.H_10
			+ lep.theta_023*shift_had.H_20
			+ lep.theta_033*shift_had.H_30
			+ lep.theta_043*shift_had.H_40));
}
Vec3 xs::rad_f_hard_up_base(Rad b, LepRadUX lep, HadUP shift_had, Transform3 shift_rot) {
	return b.coeff*dot(shift_rot, Vec3(
		// UT1.
		(lep.theta_062*shift_had.H_61 + lep.theta_082*shift_had.H_81)
		+ b.R*(lep.theta_063*shift_had.H_61 + lep.theta_083*shift_had.H_81)
		+ b.R*b.R*(lep.theta_064*shift_had.H_61 + lep.theta_084*shift_had.H_81),
		// UT2.
		(
			lep.theta_012*shift_had.H_12
			+ lep.theta_022*shift_had.H_22
			+ lep.theta_032*shift_had.H_32
			+ lep.theta_042*shift_had.H_42)
		+ b.R*(
			lep.theta_013*shift_had.H_12
			+ lep.theta_023*shift_had.H_22
			+ lep.theta_033*shift_had.H_32
			+ lep.theta_043*shift_had.H_42),
		// UL.
		(lep.theta_062*shift_had.H_63 + lep.theta_082*shift_had.H_83)
		+ b.R*(lep.theta_063*shift_had.H_63 + lep.theta_083*shift_had.H_83)
		+ b.R*b.R*(lep.theta_064*shift_had.H_63 + lep.theta_084*shift_had.H_83)));
}
Real xs::rad_f_hard_lu_base(Rad b, LepRadLU lep, HadLU shift_had) {
	return b.coeff*(
		(lep.theta_052 + lep.theta_152)*shift_had.H_50
		+ b.R*(lep.theta_053 + lep.theta_153)*shift_had.H_50);
}
Vec3 xs::rad_f_hard_lp_base(Rad b, LepRadLX lep, HadLP shift_had, Transform3 shift_rot) {
	return b.coeff*dot(shift_rot, Vec3(
		// LT1.
		(
			(lep.theta_072 + lep.theta_172)*shift_had.H_71
			+ (lep.theta_092 + lep.theta_192)*shift_had.H_91)
		+ b.R*(
			(lep.theta_073 + lep.theta_173)*shift_had.H_71
			+ (lep.theta_093 + lep.theta_193)*shift_had.H_91)
		+ b.R*b.R*(
			(lep.theta_074 + lep.theta_174)*shift_had.H_71
			+ (lep.theta_094 + lep.theta_194)*shift_had.H_91),
		// LT2.
		(lep.theta_052 + lep.theta_152)*shift_had.H_52
		+ b.R*(lep.theta_053 + lep.theta_153)*shift_had.H_52,
		// LL.
		(
			(lep.theta_072 + lep.theta_172)*shift_had.H_73
			+ (lep.theta_092 + lep.theta_192)*shift_had.H_93)
		+ b.R*(
			(lep.theta_073 + lep.theta_173)*shift_had.H_73
			+ (lep.theta_093 + lep.theta_193)*shift_had.H_93)
		+ b.R*b.R*(
			(lep.theta_074 + lep.theta_174)*shift_had.H_73
			+ (lep.theta_094 + lep.theta_194)*shift_had.H_93)));
}

// "Soft" infrared-divergent-free radiative base functions.
Real xs::rad_f_soft_uu_base(Rad b, LepRadUU lep, HadUU had, HadUU shift_had) {
	return 1./b.R*(
		rad_f_soft_uu_base_R(b, lep, shift_had)
		- rad_f_soft_uu_base_R0(b, lep, had));
}
Vec3 xs::rad_f_soft_up_base(Rad b, lep::LepRadUX lep, had::HadUP had, had::HadUP shift_had, math::Transform3 shift_rot) {
	Vec3 up_R0(
		rad_f_soft_ut1_base_R0(b, lep, had),
		rad_f_soft_ut2_base_R0(b, lep, had),
		rad_f_soft_ul_base_R0(b, lep, had));
	return 1./b.R*(rad_f_soft_up_base_R(b, lep, shift_had, shift_rot) - up_R0);
}
Real xs::rad_f_soft_lu_base(Rad b, lep::LepRadLU lep, had::HadLU had, had::HadLU shift_had) {
	return 1./b.R*(
		rad_f_soft_lu_base_R(b, lep, shift_had)
		- rad_f_soft_lu_base_R0(b, lep, had));
}
Vec3 xs::rad_f_soft_lp_base(Rad b, lep::LepRadLX lep, had::HadLP had, had::HadLP shift_had, math::Transform3 shift_rot) {
	Vec3 lp_R0(
		rad_f_soft_lt1_base_R0(b, lep, had),
		rad_f_soft_lt2_base_R0(b, lep, had),
		rad_f_soft_ll_base_R0(b, lep, had));
	return 1./b.R*(rad_f_soft_lp_base_R(b, lep, shift_had, shift_rot) - lp_R0);
}

Real xs::rad_f_soft_uu_base_R(Rad b, LepRadUU lep, HadUU shift_had) {
	return b.coeff*(
			lep.theta_011*shift_had.H_10
			+ lep.theta_021*shift_had.H_20
			+ lep.theta_031*shift_had.H_30
			+ lep.theta_041*shift_had.H_40);
}
Vec3 xs::rad_f_soft_up_base_R(Rad b, LepRadUX lep, HadUP shift_had, Transform3 shift_rot) {
	return b.coeff*dot(shift_rot, Vec3(
		// UT1.
		lep.theta_061*shift_had.H_61 + lep.theta_081*shift_had.H_81,
		// UT2.
		lep.theta_011*shift_had.H_12
		+ lep.theta_021*shift_had.H_22
		+ lep.theta_031*shift_had.H_32
		+ lep.theta_041*shift_had.H_42,
		// UL.
		lep.theta_061*shift_had.H_63 + lep.theta_081*shift_had.H_83));
}
Real xs::rad_f_soft_lu_base_R(Rad b, LepRadLU lep, HadLU shift_had) {
	return b.coeff*(lep.theta_051 + lep.theta_151)*shift_had.H_50;
}
Vec3 xs::rad_f_soft_lp_base_R(Rad b, LepRadLX lep, HadLP shift_had, Transform3 shift_rot) {
	return b.coeff*dot(shift_rot, Vec3(
		// LT1.
		(lep.theta_071 + lep.theta_171)*shift_had.H_71
		+ (lep.theta_091 + lep.theta_191)*shift_had.H_91,
		// LT2.
		(lep.theta_051 + lep.theta_151)*shift_had.H_52,
		// LL.
		(lep.theta_071 + lep.theta_171)*shift_had.H_73
		+ (lep.theta_091 + lep.theta_191)*shift_had.H_93));
}

// "Soft" infrared-divergent-free radiative base functions with hadronic part at
// `R=0`.
Real xs::rad_f_soft_uu_base_R0(Rad b, LepRadUU lep, HadUU had) {
	return rad_f_soft_uu_base_R(b, lep, had);
}
Real xs::rad_f_soft_ul_base_R0(Rad b, LepRadUP lep, HadUL had) {
	return b.coeff*(lep.theta_061*had.H_63 + lep.theta_081*had.H_83);
}
Real xs::rad_f_soft_ut1_base_R0(Rad b, LepRadUP lep, HadUT1 had) {
	return b.coeff*(lep.theta_061*had.H_61 + lep.theta_081*had.H_81);
}
Real xs::rad_f_soft_ut2_base_R0(Rad b, LepRadUU lep, HadUT2 had) {
	return b.coeff*(
		lep.theta_011*had.H_12
		+ lep.theta_021*had.H_22
		+ lep.theta_031*had.H_32
		+ lep.theta_041*had.H_42);
}
Real xs::rad_f_soft_lu_base_R0(Rad b, LepRadLU lep, HadLU had) {
	return rad_f_soft_lu_base_R(b, lep, had);
}
Real xs::rad_f_soft_ll_base_R0(Rad b, LepRadLP lep, HadLL had) {
	return b.coeff*(
		(lep.theta_071 + lep.theta_171)*had.H_73
		+ (lep.theta_091 + lep.theta_191)*had.H_93);
}
Real xs::rad_f_soft_lt1_base_R0(Rad b, LepRadLP lep, HadLT1 had) {
	return b.coeff*(
		(lep.theta_071 + lep.theta_171)*had.H_71
		+ (lep.theta_091 + lep.theta_191)*had.H_91);
}
Real xs::rad_f_soft_lt2_base_R0(Rad b, LepRadLU lep, HadLT2 had) {
	return b.coeff*(lep.theta_051 + lep.theta_151)*had.H_52;
}

