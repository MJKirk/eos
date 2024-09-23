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
            UsedParameter(p[_par_name("+", "1", "2")], *this),
            UsedParameter(p[_par_name("+", "1", "3")], *this),
            UsedParameter(p[_par_name("+", "1", "4")], *this),
            UsedParameter(p[_par_name("+", "1", "5")], *this),
            UsedParameter(p[_par_name("+", "1", "6")], *this),
            UsedParameter(p[_par_name("+", "1", "7")], *this),
            UsedParameter(p[_par_name("+", "1", "8")], *this),
            UsedParameter(p[_par_name("+", "1", "9")], *this)
        }},
        _M_fp_I1(p["pi->pi::M_(+,1)@KKRvD2024"], *this),
        _G_fp_I1(p["pi->pi::Gamma_(+,1)@KKRvD2024"], *this),
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
        const double Q2      = 1.0;

        // cf. [BL:1998A], eq. (4.6), p. 9.
        // note that the weight function factor ``(1.0 + z)^2 * (1.0 - z)^(+1/2)`` has been cancelled against the factor
        // in the redefined series expansion, i.e., ``(1.0 + z)^2 * (1.0 - z)^(+1/2) \sum_n b_n z^n = \sum_n a_n z^n``.
        return 1.0
           / sqrt(12.0 * M_PI * t_p * chi)
           * pow(tfactor, 5.0 / 4.0) * pow(sqrt(tfactor) * (1.0 + z) + (1.0 - z), -0.5)
           * pow(sqrt(1.0 + Q2 / t_p) * (1.0 - z) + sqrt(tfactor) * (1.0 + z), -3.0);
    }

    double
    KKRvD2024FormFactors<PiToPi>::phiprime_p_m1(const double & chi) const
    {
        const double t_p     = this->_t_p();
        const double t_0     = this->_t_0();
        const double tfactor = 1.0 - t_0 / t_p;
        const double Q2      = 1.0;

        return  -1 / (64 * std::sqrt(6 * M_PI * t_p * chi) * power_of<2>(Q2 + t_p))
                   * pow(tfactor, 5.0/4.0) * power_of<2>(t_p) * (-7 * std::sqrt(1.0 + Q2 / t_p) + sqrt(tfactor) * (6.0 + std::sqrt(1.0 + Q2 / t_p)));
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
    KKRvD2024FormFactors<PiToPi>::_b0_fp_I1(const double & chi, const complex<double> & zr) const
    {
        const double z0          = this->z(0.0);
        const double phi_z0      = this->phi_p(z0, chi);
        const double phi_m1      = this->phi_p(-1.0, chi);
        const double phiprime_m1 = this->phiprime_p_m1(chi);

        const double x_z0 = 1 / (phi_z0 * std::norm(z0 - zr));
        const double x_m1 = 1 / (phi_m1 * std::norm(1.0 + zr));
        const double xprime_m1 = (2.0 * (1 + std::real(zr)) * phi_m1 - std::norm(1.0 + zr) * phiprime_m1 ) / power_of<2>(std::norm(1.0 + zr) * phi_m1);

        std::array<double, 10> b;
        b[0] = 0.0;
        b[1] = 0.0;
        std::copy(_b_fp_I1.cbegin(), _b_fp_I1.cend(), b.begin()+2);

        double sum = 0.0;
        for (auto i = 0u; i < b.size(); i++)
        {
            sum += b[i] * ( (pow(-1, i) * i + pow(z0, i-1)) * x_m1 - (pow(-1, i) + pow(z0, i-1)) * xprime_m1 );
        };

        const double b0 = std::real((x_m1 - xprime_m1 - x_z0 * z0 * sum) / (x_z0 * (x_m1 - (1.0 + z0) * xprime_m1)));
        return b0;
    }

    double
    KKRvD2024FormFactors<PiToPi>::_b1_fp_I1(const double & chi, const complex<double> & zr) const
    {
        const double z0          = this->z(0.0);
        const double phi_z0      = this->phi_p(z0, chi);
        const double phi_m1      = this->phi_p(-1.0, chi);
        const double phiprime_m1 = this->phiprime_p_m1(chi);

        const double x_z0 = 1 / (phi_z0 * std::norm(z0 - zr));
        const double x_m1 = 1 / (phi_m1 * std::norm(1.0 + zr));
        const double xprime_m1 = (2.0 * (1 + std::real(zr)) * phi_m1 - std::norm(1.0 + zr) * phiprime_m1 ) / power_of<2>(std::norm(1.0 + zr) * phi_m1);

        std::array<double, 10> b;
        b[0] = 0.0;
        b[1] = 0.0;
        std::copy(_b_fp_I1.cbegin(), _b_fp_I1.cend(), b.begin()+2);

        double sum = 0.0;
        for (auto i = 0u; i < b.size(); i++)
        {
            sum += b[i] * (pow(-1, i) * i * x_m1 + (pow(-1, i+1) + pow(z0, i)) * xprime_m1);
        };

        const double b1 = (sum * x_z0 -  xprime_m1) / ( x_z0 * (x_m1 - (1.0 + z0) * xprime_m1));
        return b1;
    }

    double
    KKRvD2024FormFactors<PiToPi>::f_p(const double & q2) const
    {
        const double z   = this->z(q2);
        const double chi = 0.0258815; // GeV^-2, from Meril at Q^2 = 1 GeV^2
        const double phi = this->phi_p(z, chi);
        // the weight function has been absorbed into the outer function to cancel a superficial divergence
        // as z -> +/- 1.0

        // Super-threshold pole location
        const auto zr =  this->_zr(this->_M_fp_I1(), this->_G_fp_I1());

        // prepare expansion coefficients
        std::array<double, 10> b;
        std::copy(_b_fp_I1.cbegin(), _b_fp_I1.cend(), b.begin()+2);
        // Fix b[0] and b[1] to enforce F(q2=0) = 1 and F'(q2=t+) = 0
        b[0] = _b0_fp_I1(chi, zr);
        b[1] = _b1_fp_I1(chi, zr);

        const auto series      = this->series_m(z, b);

        return (series) / (phi * std::norm(z - zr)); // the weight factor has been absorbed into 1 / phi
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
            UsedParameter(p[_par_name("+", "1", "2")], *this),
            UsedParameter(p[_par_name("+", "1", "3")], *this),
            UsedParameter(p[_par_name("+", "1", "4")], *this),
            UsedParameter(p[_par_name("+", "1", "5")], *this),
            UsedParameter(p[_par_name("+", "1", "6")], *this),
            UsedParameter(p[_par_name("+", "1", "7")], *this),
            UsedParameter(p[_par_name("+", "1", "8")], *this),
            UsedParameter(p[_par_name("+", "1", "9")], *this)
        }},
        _M_fp_I1(p["pi->pi::M_(+,1)@KKRvD2024"], *this),
        _G_fp_I1(p["pi->pi::Gamma_(+,1)@KKRvD2024"], *this),
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
    KKRvD2024FormFactors<VacuumToPiPi>::z(const complex<double> & q2) const
    {
        return _z(q2, this->_t_0());
    }

    complex<double>
    KKRvD2024FormFactors<VacuumToPiPi>::phi_p(const complex<double> & z, const double & chi) const
    {
        const double t_p     = this->_t_p();
        const double t_0     = this->_t_0();
        const double tfactor = 1.0 - t_0 / t_p;
        const double Q2      = 1.0;

        // cf. [BL:1998A], eq. (4.6), p. 9.
        // note that the asymptotic factor ``(1.0 + z)^2 * (1.0 - z)^(+1/2)`` has been cancelled against the factor
        // in the redefined series expansion, i.e., ``(1.0 + z)^2 * (1.0 - z)^(-1/2) \sum_n b_n p_n(z) = \sum_n a_n z^n``.
        return 1.0
            / sqrt(12.0 * M_PI * t_p * chi)
            * pow(tfactor, 5.0 / 4.0) * pow(sqrt(tfactor) * (1.0 + z) + (1.0 - z), -0.5)
            * pow(sqrt(1.0 + Q2 / t_p) * (1.0 - z) + sqrt(tfactor) * (1.0 + z), -3.0);
    }

    complex<double>
    KKRvD2024FormFactors<VacuumToPiPi>::phiprime_p_m1(const double & chi) const
    {
        const double t_p     = this->_t_p();
        const double t_0     = this->_t_0();
        const double tfactor = 1.0 - t_0 / t_p;
        const double Q2      = 1.0;

        return  -1 / (64 * std::sqrt(6 * M_PI * t_p * chi) * power_of<2>(Q2 + t_p))
                   * pow(tfactor, 5.0/4.0) * power_of<2>(t_p) * (-7 * std::sqrt(1.0 + Q2 / t_p) + sqrt(tfactor) * (6.0 + std::sqrt(1.0 + Q2 / t_p)));
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

    double
    KKRvD2024FormFactors<VacuumToPiPi>::_b0_fp_I1(const double & chi, const complex<double> & zr) const
    {
        const complex<double> z0          = this->z(0.0);
        const complex<double> phi_z0      = this->phi_p(z0, chi);
        const complex<double> phi_m1      = this->phi_p(-1.0, chi);
        const complex<double> phiprime_m1 = this->phiprime_p_m1(chi);

        const complex<double> x_z0 = 1.0 / (phi_z0 * std::norm(z0 - zr));
        const complex<double> x_m1 = 1.0 / (phi_m1 * std::norm(1.0 + zr));
        const complex<double> xprime_m1 = (2.0 * (1 + std::real(zr)) * phi_m1 - std::norm(1.0 + zr) * phiprime_m1 ) / power_of<2>(std::norm(1.0 + zr) * phi_m1);

        std::array<double, 10> b;
        b[0] = 0.0;
        b[1] = 0.0;
        std::copy(_b_fp_I1.cbegin(), _b_fp_I1.cend(), b.begin()+2);

        complex<double> sum = 0.0;
        for (auto i = 0u; i < b.size(); i++)
        {
            sum += b[i] * ( (pow(-1, i) * i + pow(z0, i-1)) * x_m1 - (pow(-1, i) + pow(z0, i-1)) * xprime_m1 );
        };

        const double b0 = std::real((x_m1 - xprime_m1 - x_z0 * z0 * sum) / (x_z0 * (x_m1 - (1.0 + z0) * xprime_m1)));
        return b0;
    }

    double
    KKRvD2024FormFactors<VacuumToPiPi>::_b1_fp_I1(const double & chi, const complex<double> & zr) const
    {
        const complex<double> z0          = this->z(0.0);
        const complex<double> phi_z0      = this->phi_p(z0, chi);
        const complex<double> phi_m1      = this->phi_p(-1.0, chi);
        const complex<double> phiprime_m1 = this->phiprime_p_m1(chi);

        const complex<double> x_z0 = 1.0 / (phi_z0 * std::norm(z0 - zr));
        const complex<double> x_m1 = 1.0 / (phi_m1 * std::norm(1.0 + zr));
        const complex<double> xprime_m1 = (2.0 * (1 + std::real(zr)) * phi_m1 - std::norm(1.0 + zr) * phiprime_m1 ) / power_of<2>(std::norm(1.0 + zr) * phi_m1);

        std::array<double, 10> b;
        b[0] = 0.0;
        b[1] = 0.0;
        std::copy(_b_fp_I1.cbegin(), _b_fp_I1.cend(), b.begin()+2);

        complex<double> sum = 0.0;
        for (auto i = 0u; i < b.size(); i++)
        {
            sum += b[i] * (pow(-1, i) * i * x_m1 + (pow(-1, i+1) + pow(z0, i)) * xprime_m1);
        };

        const double b1 = std::real((sum * x_z0 -  xprime_m1) / ( x_z0 * (x_m1 - (1.0 + z0) * xprime_m1)));
        return b1;
    }

    complex<double>
    KKRvD2024FormFactors<VacuumToPiPi>::f_p(const double & q2) const
    {
        static const double eps = 1.0e-12;
        return f_p(complex<double>(q2, eps));
    }

    complex<double>
    KKRvD2024FormFactors<VacuumToPiPi>::f_p(const complex<double> & q2) const
    {
        const auto z           = this->z(q2);
        const auto chi         = 0.0258815; // GeV^-2, from Meril at Q^2 = 1 GeV^2
        const auto phi         = this->phi_p(z, chi);

        // Super-threshold pole location
        const auto zr =  this->_zr(this->_M_fp_I1(), this->_G_fp_I1());

        // prepare expansion coefficients
        std::array<double, 10> b;
        std::copy(_b_fp_I1.cbegin(), _b_fp_I1.cend(), b.begin()+2);
        // Fix b[0] and b[1] to enforce F(q2=0) = 1 and F'(q2=t+) = 0
        b[0] = _b0_fp_I1(chi, zr);
        b[1] = _b1_fp_I1(chi, zr);

        const auto series      = this->series_m(z, b);

        return series / (z - zr) / (z - std::conj(zr)) /  phi; // the weight factor has been absorbed into 1 / phi
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

    double KKRvD2024FormFactors<VacuumToPiPi>::get_b0() const
    {
        const auto chi = 0.0258815; // GeV^-2, from Meril at Q^2 = 1 GeV^2
        const auto zr  =  this->_zr(this->_M_fp_I1(), this->_G_fp_I1());
        return _b0_fp_I1(chi, zr);
    }

    complex<double> KKRvD2024FormFactors<VacuumToPiPi>::get_residue() const
    {
        // Super-threshold pole location
        const auto zr =  this->_zr(this->_M_fp_I1(), this->_G_fp_I1());

        const auto chi         = 0.0258815; // GeV^-2, from Meril at Q^2 = 1 GeV^2
        const auto phi         = this->phi_p(zr, chi);

        // prepare expansion coefficients
        std::array<double, 10> b;
        std::copy(_b_fp_I1.cbegin(), _b_fp_I1.cend(), b.begin()+1);
        // Fix b[0] to enforce F(q2=0) = 1
        b[0] = _b0_fp_I1(chi, zr);

        const auto series      = this->series_m(zr, b);

        return series / (zr - std::conj(zr)) /  phi; // the weight factor has been absorbed into 1 / phi
    }

    double KKRvD2024FormFactors<VacuumToPiPi>::get_residue_real() const
    {
        return std::real(get_residue());
    }

    double KKRvD2024FormFactors<VacuumToPiPi>::get_residue_imag() const
    {
        return std::imag(get_residue());
    }

    const std::set<ReferenceName>
    KKRvD2024FormFactors<VacuumToPiPi>::references
    {
    };

    const std::vector<OptionSpecification>
    KKRvD2024FormFactors<VacuumToPiPi>::option_specifications
    {
    };

    std::vector<OptionSpecification>::const_iterator
    KKRvD2024FormFactors<VacuumToPiPi>::begin_options()
    {
        return option_specifications.cbegin();
    }

    std::vector<OptionSpecification>::const_iterator
    KKRvD2024FormFactors<VacuumToPiPi>::end_options()
    {
        return option_specifications.cend();
    }

    template class KKRvD2024FormFactors<VacuumToPiPi>;
}
