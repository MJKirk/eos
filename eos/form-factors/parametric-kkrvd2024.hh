/* vim: set sw=4 sts=4 et tw=120 foldmethod=syntax : */

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

#ifndef EOS_GUARD_EOS_FORM_FACTORS_PARAMETRIC_KKRvD2024_HH
#define EOS_GUARD_EOS_FORM_FACTORS_PARAMETRIC_KKRvD2024_HH 1

#include <eos/form-factors/mesonic.hh>
#include <eos/form-factors/mesonic-processes.hh>
#include <eos/maths/complex.hh>
#include <eos/maths/power-of.hh>
#include <eos/utils/parameters.hh>
#include <eos/utils/reference-name.hh>

#include <array>

namespace eos
{
    template <typename Process_> class KKRvD2024FormFactors;

    template <> class KKRvD2024FormFactors<PiToPi> :
        public FormFactors<PToP>
    {
        private:
            // parameters for form factor f_+ (I=1 projection)
            std::array<UsedParameter, 10u> _b_fp_I1;
            UsedParameter _re_c_fp_I1;
            UsedParameter _im_c_fp_I1;
            UsedParameter _M_fp_I1;
            UsedParameter _G_fp_I1;

            // hadron masses
            UsedParameter _m_pi;

            // parameter for zero point of z
            UsedParameter _t_0;

            inline std::string _par_name(const std::string & ff, const std::string & isospin, const std::string & index) const
            {
                return "pi->pi::b_(" + ff + "," + isospin + ")^" + index + "@KKRvD2024";
            }

            inline double _t_p() const
            {
                return 4.0 * _m_pi() * _m_pi();
            }

            inline double _z(const double & q2, const double & t_0) const
            {
                const auto t_p = _t_p();
                return (sqrt(t_p - q2) - sqrt(t_p - t_0)) / (sqrt(t_p - q2) + sqrt(t_p - t_0));
            }

            inline complex<double> _z(const complex<double> & q2, const double & t_0) const
            {
                const auto t_p = _t_p();
                return (sqrt(t_p - q2) - sqrt(t_p - t_0)) / (sqrt(t_p - q2) + sqrt(t_p - t_0));
            }

            inline complex<double> _zr(const double & M, const double & Gamma) const
            {
                return 1.0 / _z(power_of<2>(complex<double>(M, -Gamma/2)), _t_0());
            }

            inline complex<double> _inverseblaschke(const complex<double> & z, const complex<double> & x) const
            {
                return (x / abs(x)) * (1.0 - z * std::conj(x)) / (z - x);
            }

        public:
            KKRvD2024FormFactors(const Parameters & p, const Options & o);
            ~KKRvD2024FormFactors();

            static FormFactors<PToP> * make(const Parameters & p, const Options & o);

            /* auxiliary functions */
            double z(const double & q2) const;
            double phi_p(const double & z, const double & chi) const;
            double series_m(const double & z, const std::array<double, 10u> & c) const;

            virtual double f_p(const double & q2) const override;
            virtual double f_plus_T(const double & q2) const override;
            virtual double f_t(const double & q2) const override;
            virtual double f_0(const double & q2) const override;
    };

    extern template class KKRvD2024FormFactors<PiToPi>;

    template <> class KKRvD2024FormFactors<VacuumToPiPi> :
        public FormFactors<VacuumToPP>
    {
        private:
            // parameters for form factor f_+ (I=1 projection)
            std::array<UsedParameter, 10u> _b_fp_I1;
            UsedParameter _re_c_fp_I1;
            UsedParameter _im_c_fp_I1;
            UsedParameter _M_fp_I1;
            UsedParameter _G_fp_I1;

            // hadron masses
            UsedParameter _m_pi;

            // parameter for zero point of z
            UsedParameter _t_0;

            inline std::string _par_name(const std::string & ff, const std::string & isospin, const std::string & index) const
            {
                return "pi->pi::b_(" + ff + "," + isospin + ")^" + index + "@KKRvD2024";
            }

            inline double _t_p() const
            {
                return 4.0 * _m_pi() * _m_pi();
            }

            inline complex<double> _z(const double & q2, const double & t_0) const
            {
                const complex<double> t_p = complex<double>{ _t_p(), 0.0};

                return (sqrt(t_p - q2) - sqrt(t_p - t_0)) / (sqrt(t_p - q2) + sqrt(t_p - t_0));
            }

            inline complex<double> _z(const complex<double> & q2, const double & t_0) const
            {
                const auto t_p = _t_p();
                return (sqrt(t_p - q2) - sqrt(t_p - t_0)) / (sqrt(t_p - q2) + sqrt(t_p - t_0));
            }

            inline complex<double> _zr(const double & M, const double & Gamma) const
            {
                return 1.0 / _z(power_of<2>(complex<double>(M, -Gamma/2)), _t_0());
            }

            inline complex<double> _inverseblaschke(const complex<double> & z, const complex<double> & x) const
            {
                return (x / abs(x)) * (1.0 - z * std::conj(x)) / (z - x);
            }

        public:
            KKRvD2024FormFactors(const Parameters & p, const Options & o);
            ~KKRvD2024FormFactors();

            static FormFactors<VacuumToPP> * make(const Parameters & p, const Options & o);

            /* auxiliary functions */
            complex<double> z(const double & q2) const;
            complex<double> phi_p(const complex<double> & z, const double & chi) const;
            complex<double> series_m(const complex<double> & z, const std::array<double, 10u> & c) const;

            virtual complex<double> f_p(const double & q2) const;
            virtual complex<double> f_t(const double & q2) const;
            virtual complex<double> f_0(const double & q2) const;
    };

    extern template class KKRvD2024FormFactors<VacuumToPiPi>;
}

#endif
