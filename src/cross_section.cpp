#include "sidis/cross_section.hpp"

#include <cmath>
#include <limits>
#include <string>

#include <cubature.hpp>

#include "sidis/constant.hpp"
#include "sidis/cut.hpp"
#include "sidis/frame.hpp"
#include "sidis/hadronic_coeff.hpp"
#include "sidis/kinematics.hpp"
#include "sidis/leptonic_coeff.hpp"
#include "sidis/structure_function.hpp"
#include "sidis/extra/bounds.hpp"
#include "sidis/extra/math.hpp"
#include "sidis/extra/integrate.hpp"
#include "sidis/extra/vector.hpp"

using namespace sidis;
using namespace sidis::constant;
using namespace sidis::cut;
using namespace sidis::had;
using namespace sidis::integ;
using namespace sidis::kin;
using namespace sidis::lep;
using namespace sidis::math;
using namespace sidis::sf;
using namespace sidis::xs;

namespace {

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
		lambda_S_prime*sq(L_S_prime)/4. - lambda_X_prime*sq(L_X_prime)/4.
		+ dilog(1. - 1./(sum_S_prime*S_prime)*rho)
		+ dilog(1. - (sum_S_prime*S_prime)/(4.*sq(kin.m)*kin.mx_sq)*rho)
		- dilog(1. - (sum_X_prime*X_prime)/(kin.mx_sq*sq(sum_m)*kin.Q_sq)*rho)
		- dilog(1. - (4.*sq(kin.m))/(sum_X_prime*X_prime*sq(sum_m)*kin.Q_sq)*rho));

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

Real xs::born(Real lambda_e, Vec3 eta, Kinematics kin, SfSet const& model) {
	Born b(kin);
	LepBornXX lep(kin);
	HadXX had(kin, model);
	return born_xx_base(lambda_e, eta, b, lep, had);
}

Real xs::amm(Real lambda_e, Vec3 eta, Kinematics kin, SfSet const& model) {
	Amm b(kin);
	LepAmmXX lep(kin);
	HadXX had(kin, model);
	return amm_xx_base(lambda_e, eta, b, lep, had);
}

Real xs::nrad_ir(Real lambda_e, Vec3 eta, Kinematics kin, SfSet const& model, Real k0_cut) {
	NRadIR b(kin, k0_cut);
	LepBornXX lep_born(kin);
	LepAmmXX lep_amm(kin);
	HadXX had(kin, model);
	return nrad_ir_xx_base(lambda_e, eta, b, lep_born, lep_amm, had);
}

Real xs::rad(Real lambda_e, Vec3 eta, KinematicsRad kin, SfSet const& model) {
	Rad b(kin);
	LepRadXX lep(kin);
	HadRadXX had(kin, model);
	return rad_xx_base(lambda_e, eta, b, lep, had);
}

Real xs::rad_f(Real lambda_e, Vec3 eta, KinematicsRad kin, SfSet const& model) {
	Rad b(kin);
	LepRadXX lep(kin);
	HadRadFXX had(kin, model);
	return rad_f_xx_base(lambda_e, eta, b, lep, had);
}

Real xs::nrad(Real lambda_e, Vec3 eta, Kinematics kin, SfSet const& model, Real k0_cut) {
	// The soft part of the radiative cross-section (below `k0_cut`) is bundled
	// into the return value here.
	Real xs_nrad_ir = nrad_ir(lambda_e, eta, kin, model, k0_cut);
	Real xs_rad_f = rad_f_integ(lambda_e, eta, kin, model, k0_cut);
	return xs_nrad_ir + xs_rad_f;
}

Real xs::rad_f_integ(Real lambda_e, Vec3 eta, Kinematics kin, SfSet const& model, Real k0_cut) {
	HadXX had_0(kin, model);
	CutRad cut;
	cut.k_0_bar = Bounds(0., k0_cut);
	cubature::EstErr<Real> xs_integ = cubature::cubature<3>(
		[&](cubature::Point<3, Real> x) {
			Real point[3] = { x[0], x[1], x[2] };
			KinematicsRad kin_rad;
			Real jacobian;
			while (!take(cut, kin, point, &kin_rad, &jacobian)) { }

			Rad b(kin_rad);
			LepRadXX lep(kin_rad);
			HadRadFXX had(kin_rad, model, had_0);
			Real xs = rad_f_xx_base(lambda_e, eta, b, lep, had);
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
		1000000, 0., 1e-6);
	return xs_integ.val;
}

Real xs::rad_integ(Real lambda_e, Vec3 eta, Kinematics kin, SfSet const& model, Real k0_cut) {
	CutRad cut;
	cut.k_0_bar = Bounds(k0_cut, INF);
	cubature::EstErr<Real> xs_integ = cubature::cubature<3>(
		[&](cubature::Point<3, Real> x) {
			Real point[3] = { x[0], x[1], x[2] };
			KinematicsRad kin_rad;
			Real jacobian;
			while (!take(cut, kin, point, &kin_rad, &jacobian)) { }

			Rad b(kin_rad);
			LepRadXX lep(kin_rad);
			HadRadXX had(kin_rad, model);
			Real xs = rad_xx_base(lambda_e, eta, b, lep, had);
			if (std::isnan(xs)) {
				// TODO: Handle this case more correctly.
				return 0.;
			} else {
				return jacobian * xs;
			}
		},
		cubature::Point<3, Real>{ 0., 0., 0. },
		cubature::Point<3, Real>{ 1., 1., 1. },
		1000000, 0., 1e-6);
	return xs_integ.val;
}

// Radiative corrections to Born cross-section.
Real xs::delta_vert_rad_ir(Kinematics kin, Real k0_cut) {
	// Paragraph following equation [1.C17].
	Real k0_max = (kin.mx_sq - sq(kin.Mth))/(2.*kin.mx);
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
		(kin.mx_sq - sq(kin.Mth))/(kin.m*kin.mx));
	return delta_0 + delta_shift;
}
Real xs::delta_rad_ir_hard(Kinematics kin, Real k0_cut) {
	Real k0_max = (kin.mx_sq - sq(kin.Mth))/(2.*kin.mx);
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
		(kin.mx_sq - sq(kin.Mth))/(2.*k0_cut*kin.mx));
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

Real xs::born_xx_base(Real lambda_e, Vec3 eta, Born b, LepBornXX lep, HadXX had) {
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
Real xs::born_ut1_base(Born b, LepBornUP lep, HadUT had) {
	return b.coeff*(lep.theta_6*had.H_61 + lep.theta_8*had.H_81);
}
Real xs::born_ut2_base(Born b, LepBornUU lep, HadUT had) {
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
Real xs::born_lt1_base(Born b, LepBornLP lep, HadLT had) {
	return b.coeff*(lep.theta_7*had.H_71 + lep.theta_9*had.H_91);
}
Real xs::born_lt2_base(Born b, LepBornLU lep, HadLT had) {
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

Real xs::amm_xx_base(Real lambda_e, Vec3 eta, Amm b, LepAmmXX lep, HadXX had) {
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
Real xs::amm_ut1_base(Amm b, LepAmmUP lep, HadUT had) {
	return b.coeff*(lep.theta_6*had.H_61 + lep.theta_8*had.H_81);
}
Real xs::amm_ut2_base(Amm b, LepAmmUU lep, HadUT had) {
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
Real xs::amm_lt1_base(Amm b, LepAmmLP lep, HadLT had) {
	return b.coeff*(lep.theta_7*had.H_71 + lep.theta_9*had.H_91);
}
Real xs::amm_lt2_base(Amm b, LepAmmLU lep, HadLT had) {
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

Real xs::nrad_ir_xx_base(Real lambda_e, math::Vec3 eta, NRadIR b, LepBornXX lep_born, LepAmmXX lep_amm, HadXX had) {
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
Real xs::nrad_ir_ut1_base(NRadIR b, LepBornUP lep_born, LepAmmUP lep_amm, HadUT had) {
	return
		(b.coeff_born*lep_born.theta_6 + b.coeff_amm*lep_amm.theta_6)*had.H_61
		+ (b.coeff_born*lep_born.theta_8 + b.coeff_amm*lep_amm.theta_8)*had.H_81;
}
Real xs::nrad_ir_ut2_base(NRadIR b, LepBornUU lep_born, LepAmmUU lep_amm, HadUT had) {
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
Real xs::nrad_ir_lt1_base(NRadIR b, LepBornLP lep_born, LepAmmLP lep_amm, HadLT had) {
	return
		(b.coeff_born*lep_born.theta_7 + b.coeff_amm*lep_amm.theta_7)*had.H_71
		+ (b.coeff_born*lep_born.theta_9 + b.coeff_amm*lep_amm.theta_9)*had.H_91;
}
Real xs::nrad_ir_lt2_base(NRadIR b, LepBornLU lep_born, LepAmmLU lep_amm, HadLT had) {
	return (b.coeff_born*lep_born.theta_5 + b.coeff_amm*lep_amm.theta_5)*had.H_52;
}

// Radiative base functions.
Rad::Rad(KinematicsRad kin) {
	// Equation [1.43].
	coeff = -(std::pow(ALPHA, 3)*kin.S*sq(kin.S_x))
		/(64.*sq(PI)*kin.M*kin.ph_l*kin.lambda_S*kin.lambda_Y_sqrt);
	R = kin.R;
}

Real xs::rad_xx_base(Real lambda_e, Vec3 eta, Rad b, LepRadXX lep, HadRadXX had) {
	Real uu = rad_uu_base(b, lep, had);
	Vec3 up = rad_up_base(b, lep, had);
	Real lu = rad_lu_base(b, lep, had);
	Vec3 lp = rad_lp_base(b, lep, had);
	return uu + dot(eta, up) + lambda_e * (lu + dot(eta, lp));
}
Real xs::rad_uu_base(Rad b, LepRadUU lep, HadRadUU had) {
	return b.coeff*(
		1./b.R*(
			lep.theta_011*had.H_10
			+ lep.theta_021*had.H_20
			+ lep.theta_031*had.H_30
			+ lep.theta_041*had.H_40)
		+ (
			lep.theta_012*had.H_10
			+ lep.theta_022*had.H_20
			+ lep.theta_032*had.H_30
			+ lep.theta_042*had.H_40)
		+ b.R*(
			lep.theta_013*had.H_10
			+ lep.theta_023*had.H_20
			+ lep.theta_033*had.H_30
			+ lep.theta_043*had.H_40));
}
Vec3 xs::rad_up_base(Rad b, LepRadUX lep, HadRadUP had) {
	return b.coeff*(
		1./b.R*(
			lep.theta_011*had.H_1
			+ lep.theta_021*had.H_2
			+ lep.theta_031*had.H_3
			+ lep.theta_041*had.H_4
			+ lep.theta_061*had.H_6
			+ lep.theta_081*had.H_8)
		+ (
			lep.theta_012*had.H_1
			+ lep.theta_022*had.H_2
			+ lep.theta_032*had.H_3
			+ lep.theta_042*had.H_4
			+ lep.theta_062*had.H_6
			+ lep.theta_082*had.H_8)
		+ b.R*(
			lep.theta_013*had.H_1
			+ lep.theta_023*had.H_2
			+ lep.theta_033*had.H_3
			+ lep.theta_043*had.H_4
			+ lep.theta_063*had.H_6
			+ lep.theta_083*had.H_8)
		+ b.R*b.R*(
			lep.theta_064*had.H_6
			+ lep.theta_084*had.H_8));
}
Real xs::rad_lu_base(Rad b, LepRadLU lep, HadRadLU had) {
	return b.coeff*(
		1./b.R*(lep.theta_051 + lep.theta_151)*had.H_50
		+ (lep.theta_052 + lep.theta_152)*had.H_50
		+ b.R*(lep.theta_053 + lep.theta_153)*had.H_50);
}
Vec3 xs::rad_lp_base(Rad b, LepRadLX lep, HadRadLP had) {
	return b.coeff*(
		1./b.R*(
			(lep.theta_051 + lep.theta_151)*had.H_5
			+ (lep.theta_071 + lep.theta_171)*had.H_7
			+ (lep.theta_091 + lep.theta_191)*had.H_9)
		+ (
		 	(lep.theta_052 + lep.theta_152)*had.H_5
			+ (lep.theta_072 + lep.theta_172)*had.H_7
			+ (lep.theta_092 + lep.theta_192)*had.H_9)
		+ b.R*(
		 	(lep.theta_053 + lep.theta_153)*had.H_5
			+ (lep.theta_073 + lep.theta_173)*had.H_7
			+ (lep.theta_093 + lep.theta_193)*had.H_9)
		+ b.R*b.R*(
			(lep.theta_074 + lep.theta_174)*had.H_7
			+ (lep.theta_094 + lep.theta_194)*had.H_9));
}

Real xs::rad_f_xx_base(Real lambda_e, Vec3 eta, Rad b, LepRadXX lep, HadRadFXX had) {
	Real uu = rad_f_uu_base(b, lep, had);
	Vec3 up = rad_f_up_base(b, lep, had);
	Real lu = rad_f_lu_base(b, lep, had);
	Vec3 lp = rad_f_lp_base(b, lep, had);
	return uu + dot(eta, up) + lambda_e * (lu + dot(eta, lp));
}
Real xs::rad_f_uu_base(Rad b, LepRadUU lep, HadRadFUU had) {
	return b.coeff*(
		(
			lep.theta_011*had.H_10_diff
			+ lep.theta_021*had.H_20_diff
			+ lep.theta_031*had.H_30_diff
			+ lep.theta_041*had.H_40_diff)
		+ (
			lep.theta_012*had.H_10
			+ lep.theta_022*had.H_20
			+ lep.theta_032*had.H_30
			+ lep.theta_042*had.H_40)
		+ b.R*(
			lep.theta_013*had.H_10
			+ lep.theta_023*had.H_20
			+ lep.theta_033*had.H_30
			+ lep.theta_043*had.H_40));
}
Vec3 xs::rad_f_up_base(Rad b, LepRadUX lep, HadRadFUP had) {
	return b.coeff*(
		(
			lep.theta_011*had.H_1_diff
			+ lep.theta_021*had.H_2_diff
			+ lep.theta_031*had.H_3_diff
			+ lep.theta_041*had.H_4_diff
			+ lep.theta_061*had.H_6_diff
			+ lep.theta_081*had.H_8_diff)
		+ (
			lep.theta_012*had.H_1
			+ lep.theta_022*had.H_2
			+ lep.theta_032*had.H_3
			+ lep.theta_042*had.H_4
			+ lep.theta_062*had.H_6
			+ lep.theta_082*had.H_8)
		+ b.R*(
			lep.theta_013*had.H_1
			+ lep.theta_023*had.H_2
			+ lep.theta_033*had.H_3
			+ lep.theta_043*had.H_4
			+ lep.theta_063*had.H_6
			+ lep.theta_083*had.H_8)
		+ b.R*b.R*(
			lep.theta_064*had.H_6
			+ lep.theta_084*had.H_8));
}
Real xs::rad_f_lu_base(Rad b, LepRadLU lep, HadRadFLU had) {
	return b.coeff*(
		(lep.theta_051 + lep.theta_151)*had.H_50_diff
		+ (lep.theta_052 + lep.theta_152)*had.H_50
		+ b.R*(lep.theta_053 + lep.theta_153)*had.H_50);
}
Vec3 xs::rad_f_lp_base(Rad b, LepRadLX lep, HadRadFLP had) {
	return b.coeff*(
		(
			(lep.theta_051 + lep.theta_151)*had.H_5_diff
			+ (lep.theta_071 + lep.theta_171)*had.H_7_diff
			+ (lep.theta_091 + lep.theta_191)*had.H_9_diff)
		+ (
		 	(lep.theta_052 + lep.theta_152)*had.H_5
			+ (lep.theta_072 + lep.theta_172)*had.H_7
			+ (lep.theta_092 + lep.theta_192)*had.H_9)
		+ b.R*(
		 	(lep.theta_053 + lep.theta_153)*had.H_5
			+ (lep.theta_073 + lep.theta_173)*had.H_7
			+ (lep.theta_093 + lep.theta_193)*had.H_9)
		+ b.R*b.R*(
			(lep.theta_074 + lep.theta_174)*had.H_7
			+ (lep.theta_094 + lep.theta_194)*had.H_9));
}

