/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2020-2024 Danny van Dyk
 * Copyright (c) 2024 Matthew J. Kirk
 *
 * This file is part of the EOS project. EOS is free software;
 * you can redistribute it and/or modify it under the terms of the GNU General
 * Public License version 2, as published by the Free Software Foundation.
 *
 * EOS is distributed in the hope that it will be useful, but WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
 * details.
 *
 * You should have received a copy of the GNU General Public License along with
 * this program; if not, write to the Free Software Foundation, Inc., 59 Temple
 * Place, Suite 330, Boston, MA  02111-1307  USA
 */

#include <eos/form-factors/parametric-kkrvd2024.hh>
#include <eos/maths/power-of.hh>
#include <eos/utils/parameters.hh>
#include <eos/utils/options.hh>
#include <eos/utils/options-impl.hh>
#include <eos/utils/qualified-name.hh>
#include <eos/utils/stringify.hh>

#include <numeric>

namespace eos
{
    KKRvD2024FormFactors<PiToPi>::KKRvD2024FormFactors(const Parameters & p, const Options & /*o*/) :
        _b_fp_I1{{
            UsedParameter(p[_par_name("+", "1", "0")], *this),
            UsedParameter(p[_par_name("+", "1", "1")], *this),
            UsedParameter(p[_par_name("+", "1", "2")], *this),
            UsedParameter(p[_par_name("+", "1", "3")], *this),
            UsedParameter(p[_par_name("+", "1", "4")], *this),
            UsedParameter(p[_par_name("+", "1", "5")], *this),
            UsedParameter(p[_par_name("+", "1", "6")], *this),
            UsedParameter(p[_par_name("+", "1", "7")], *this),
            UsedParameter(p[_par_name("+", "1", "8")], *this),
            UsedParameter(p[_par_name("+", "1", "9")], *this)
        }},
        _re_cr(p["pi->pi::Re{c}_(+,1)@KKRvD2024"], *this),
        _im_cr(p["pi->pi::Im{c}_(+,1)@KKRvD2024"], *this),
        _Mr(p["pi->pi::M_(+,1)@KKRvD2024"], *this),
        _Gr(p["pi->pi::Gamma_(+,1)@KKRvD2024"], *this),
        _m_pi(p["mass::pi^+"], *this),
        _t_0(p["pi->pi::t_0@KKRvD2024"], *this)
    {
    }

    KKRvD2024FormFactors<PiToPi>::~KKRvD2024FormFactors() = default;

    FormFactors<PToP> *
    KKRvD2024FormFactors<PiToPi>::make(const Parameters & p, const Options & o)
    {
        return new KKRvD2024FormFactors<PiToPi>(p, o);
    }

    double
    KKRvD2024FormFactors<PiToPi>::z(const double & q2) const
    {
        return _z(q2, this->_t_0());
    }

    double
    KKRvD2024FormFactors<PiToPi>::phi_p(const double & z, const double & chi) const
    {
        const double t_p     = this->_t_p();
        const double t_0     = this->_t_0();
        const double tfactor = 1.0 - t_0 / t_p;
        const double Q2      = 2.0;

        // cf. [BL:1998A], eq. (4.6), p. 9.
        // note that the weight function factor ``(1.0 + z)^2 * (1.0 - z)^(+1/2)`` has been cancelled against the factor
        // in the redefined series expansion, i.e., ``(1.0 + z)^2 * (1.0 - z)^(+1/2) \sum_n b_n z^n = \sum_n a_n z^n``.
        return 1.0
           / sqrt(12.0 * M_PI * t_p * chi)
           * pow(tfactor, 5.0 / 4.0) * pow(sqrt(tfactor) * (1.0 + z) + (1.0 - z), -0.5)
           * pow(sqrt(1.0 + Q2 / t_p) * (1.0 - z) + sqrt(tfactor) * (1.0 + z), -3.0);
    }

    double
    KKRvD2024FormFactors<PiToPi>::series_m(const double & z, const std::array<double, 10u> & c) const
    {
        std::array<double, 10> zvalues;
        double current_zv = 1.0;
        for (double & zv: zvalues)
        {
            zv = current_zv;
            current_zv *= z;
        }

        return std::inner_product(c.cbegin(), c.cend(), zvalues.cbegin(), 0.0);
    }

    double
    KKRvD2024FormFactors<PiToPi>::f_p(const double & q2) const
    {
        // prepare expansion coefficients
        std::array<double, 10> b;
        std::copy(_b_fp_I1.cbegin(), _b_fp_I1.cend(), b.begin());

        const auto z           = this->z(q2);
        const auto chi         = 3.52e-3; // GeV^-2, cf. [BL:1998A], p. 13
        const auto phi         = this->phi_p(z, chi);
        // the weight function has been absorbed into the outer function to cancel a superficial divergence
        // as z -> +/- 1.0
        const auto series      = this->series_m(z, b);

        // Super-threshold pole location
        const auto zr =  this->_zr(this->_Mr(), this->_Gr());
        // Artificial FF zero interpolation value
        const auto cr = complex<double>(this->_re_cr(), this->_im_cr());

        // Inverse Blaschke factors
        const auto B1 = this->_inverseblaschke(z, zr);

        return (series * power_of<2>(std::abs(B1)) + 2 * (cr * B1).real()) /  phi; // the weight factor has been absorbed into 1 / phi
    }

    double
    KKRvD2024FormFactors<PiToPi>::f_plus_T(const double & /*q2*/) const
    {
        throw InternalError("Not implemented!");
        return 1.0;
    }

    double
    KKRvD2024FormFactors<PiToPi>::f_t(const double & /*q2*/) const
    {
        throw InternalError("Not implemented!");
        return 1.0;
    }

    double
    KKRvD2024FormFactors<PiToPi>::f_0(const double & /*q2*/) const
    {
        return 0.0; // vanishes in our approximation
    }

    template class KKRvD2024FormFactors<PiToPi>;

    /* Vacuum -> pi pi */
    KKRvD2024FormFactors<VacuumToPiPi>::KKRvD2024FormFactors(const Parameters & p, const Options & /*o*/) :
        _b_fp_I1{{
            UsedParameter(p[_par_name("+", "1", "0")], *this),
            UsedParameter(p[_par_name("+", "1", "1")], *this),
            UsedParameter(p[_par_name("+", "1", "2")], *this),
            UsedParameter(p[_par_name("+", "1", "3")], *this),
            UsedParameter(p[_par_name("+", "1", "4")], *this),
            UsedParameter(p[_par_name("+", "1", "5")], *this),
            UsedParameter(p[_par_name("+", "1", "6")], *this),
            UsedParameter(p[_par_name("+", "1", "7")], *this),
            UsedParameter(p[_par_name("+", "1", "8")], *this),
            UsedParameter(p[_par_name("+", "1", "9")], *this)
        }},
        _re_cr(p["pi->pi::Re{c}_(+,1)@KKRvD2024"], *this),
        _im_cr(p["pi->pi::Im{c}_(+,1)@KKRvD2024"], *this),
        _Mr(p["pi->pi::M_(+,1)@KKRvD2024"], *this),
        _Gr(p["pi->pi::Gamma_(+,1)@KKRvD2024"], *this),
        _m_pi(p["mass::pi^+"], *this),
        _t_0(p["pi->pi::t_0@KKRvD2024"], *this)
    {
    }

    KKRvD2024FormFactors<VacuumToPiPi>::~KKRvD2024FormFactors() = default;

    FormFactors<VacuumToPP> *
    KKRvD2024FormFactors<VacuumToPiPi>::make(const Parameters & p, const Options & o)
    {
        return new KKRvD2024FormFactors<VacuumToPiPi>(p, o);
    }

    complex<double>
    KKRvD2024FormFactors<VacuumToPiPi>::z(const double & q2) const
    {
        return _z(q2, this->_t_0());
    }

    complex<double>
    KKRvD2024FormFactors<VacuumToPiPi>::phi_p(const complex<double> & z, const double & chi) const
    {
        const double t_p     = this->_t_p();
        const double t_0     = this->_t_0();
        const double tfactor = 1.0 - t_0 / t_p;
        const double Q2      = 2.0;

        // cf. [BL:1998A], eq. (4.6), p. 9.
        // note that the asymptotic factor ``(1.0 + z)^2 * (1.0 - z)^(+1/2)`` has been cancelled against the factor
        // in the redefined series expansion, i.e., ``(1.0 + z)^2 * (1.0 - z)^(-1/2) \sum_n b_n p_n(z) = \sum_n a_n z^n``.
        return 1.0
            / sqrt(12.0 * M_PI * t_p * chi)
            * pow(tfactor, 5.0 / 4.0) * pow(sqrt(tfactor) * (1.0 + z) + (1.0 - z), -0.5)
            * pow(sqrt(1.0 + Q2 / t_p) * (1.0 - z) + sqrt(tfactor) * (1.0 + z), -3.0);
    }

    complex<double>
    KKRvD2024FormFactors<VacuumToPiPi>::series_m(const complex<double> & z, const std::array<double, 10u> & c) const
    {
        std::array<complex<double>, 10> zvalues;
        complex<double> current_zv = 1.0;
        for (complex<double> & zv: zvalues)
        {
            zv = current_zv;
            current_zv *= z;
        }

        return std::inner_product(c.cbegin(), c.cend(), zvalues.cbegin(), complex<double>(0.0, 0.0));
    }

    complex<double>
    KKRvD2024FormFactors<VacuumToPiPi>::f_p(const double & q2) const
    {
        // prepare expansion coefficients
        std::array<double, 10> b;
        std::copy(_b_fp_I1.cbegin(), _b_fp_I1.cend(), b.begin());

        const auto z           = this->z(q2);
        const auto chi         = 3.52e-3; // GeV^-2, cf. [BL:1998A], p. 13
        const auto phi         = this->phi_p(z, chi);
        const auto series      = this->series_m(z, b);

        // Super-threshold pole location
        const auto zr =  this->_zr(this->_Mr(), this->_Gr());
        // Artificial FF zero interpolation value
        const auto cr = complex<double>(this->_re_cr(), this->_im_cr());

        // Inverse Blaschke factors
        const auto B1 = this->_inverseblaschke(z, zr);
        const auto B2 = this->_inverseblaschke(z, std::conj(zr));

        return (series * B1 * B2 + cr * B1 + std::conj(cr) * B2) /  phi; // the weight factor has been absorbed into 1 / phi
    }

    complex<double>
    KKRvD2024FormFactors<VacuumToPiPi>::f_t(const double & /*q2*/) const
    {
        throw InternalError("Not implemented!");
        return 0.0;
    }

    complex<double>
    KKRvD2024FormFactors<VacuumToPiPi>::f_0(const double & /*q2*/) const
    {
        return 0.0; // vanishes in our approximation
    }

    template class KKRvD2024FormFactors<VacuumToPiPi>;
}
