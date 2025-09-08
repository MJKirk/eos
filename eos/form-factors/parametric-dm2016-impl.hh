/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2014-2023 Danny van Dyk
 * Copyright (c) 2022 Méril Reboud
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

#ifndef EOS_GUARD_EOS_FORM_FACTORS_PARAMETRIC_DM2016_IMPL_HH
#define EOS_GUARD_EOS_FORM_FACTORS_PARAMETRIC_DM2016_IMPL_HH 1

#include <eos/form-factors/parametric-dm2016.hh>

namespace eos
{
    template <typename Process_>
    class DM2016FormFactorTraits :
        public virtual ParameterUser
    {
        public:
            // The following parameters are part of the parameterization and must match the
            // the ones used for the extraction of the coefficients of the z-expension
            const double m_1, m_2; // m_1 is the mass of the heavier particle, m_2 the mass of the lighter particle
            const double m_R_0m, m_R_0p, m_R_1m, m_R_1p;
            const double tp; // pair production thresholds
            const double tm; // kinematic endpoint
            const double t0; // z(t_0) = 0, t_m is the endpoint of the semileptonic process

            static const std::map<std::tuple<QuarkFlavor, QuarkFlavor>, double> resonance_0m_masses;
            static const std::map<std::tuple<QuarkFlavor, QuarkFlavor>, double> resonance_0p_masses;
            static const std::map<std::tuple<QuarkFlavor, QuarkFlavor>, double> resonance_1m_masses;
            static const std::map<std::tuple<QuarkFlavor, QuarkFlavor>, double> resonance_1p_masses;
            static const std::map<std::tuple<QuarkFlavor, QuarkFlavor>, double> threshold_tp_values;
            static const std::map<std::tuple<QuarkFlavor, QuarkFlavor>, double> config_t0_values;

            DM2016FormFactorTraits(const Parameters & /* p */) :
                m_1(Process_::m1),
                m_2(Process_::m2),
                m_R_0m(resonance_0m_masses.at(Process_::partonic_transition)),
                m_R_0p(resonance_0p_masses.at(Process_::partonic_transition)),
                m_R_1m(resonance_1m_masses.at(Process_::partonic_transition)),
                m_R_1p(resonance_1p_masses.at(Process_::partonic_transition)),
                tp(threshold_tp_values.at(Process_::partonic_transition)),
                tm((m_1 - m_2) * (m_1 - m_2)),
                t0(tm)
            {
            }

            complex<double> calc_z(const complex<double> & t, const complex<double> & tp, const complex<double> & t0) const
            {
                return (std::sqrt(tp - t) - std::sqrt(tp - t0)) / (std::sqrt(tp - t) + std::sqrt(tp - t0));
            }

            double calc_z(const double & t, const double & tp, const double & t0) const
            {
                if (t > tp)
                    throw InternalError("The real conformal mapping is used above threshold: " + stringify(t) + " > " + stringify(tp));

                return real(calc_z(complex<double>(t, 0.0), complex<double>(tp, 0.0), complex<double>(t0, 0.0)));
            }

            double calc_z(const double & t) const
            {
                return calc_z(t, this->tp, this->tm);
            }
    };

    template <typename Process_>
    const std::map<std::tuple<QuarkFlavor, QuarkFlavor>, double>
    DM2016FormFactorTraits<Process_>::resonance_0m_masses
    {
        { std::make_tuple(QuarkFlavor::bottom, QuarkFlavor::strange), 5.367 },
    };

    template <typename Process_>
    const std::map<std::tuple<QuarkFlavor, QuarkFlavor>, double>
    DM2016FormFactorTraits<Process_>::resonance_0p_masses
    {
        { std::make_tuple(QuarkFlavor::bottom, QuarkFlavor::strange), 5.711 },
    };

    template <typename Process_>
    const std::map<std::tuple<QuarkFlavor, QuarkFlavor>, double>
    DM2016FormFactorTraits<Process_>::resonance_1m_masses
    {
        { std::make_tuple(QuarkFlavor::bottom, QuarkFlavor::strange), 5.416 },
    };

    template <typename Process_>
    const std::map<std::tuple<QuarkFlavor, QuarkFlavor>, double>
    DM2016FormFactorTraits<Process_>::resonance_1p_masses
    {
        { std::make_tuple(QuarkFlavor::bottom, QuarkFlavor::strange), 5.750 },
    };

    template <typename Process_>
    const std::map<std::tuple<QuarkFlavor, QuarkFlavor>, double>
    DM2016FormFactorTraits<Process_>::threshold_tp_values
    {
        { std::make_tuple(QuarkFlavor::bottom, QuarkFlavor::strange), (5.279 + 0.494) * (5.279 + 0.494) },
    };

    template <typename Process_>
    DM2016FormFactors<Process_>::DM2016FormFactors(const Parameters & p, const Options &) :
        // time, V
        _alpha_0_time_v(p[stringify(Process_::label) + "::a_0_time^V@DM2016"], *this),
        _alpha_1_time_v(p[stringify(Process_::label) + "::a_1_time^V@DM2016"], *this),
        _alpha_2_time_v(p[stringify(Process_::label) + "::a_2_time^V@DM2016"], *this),
        // time, A
        _alpha_0_time_a(p[stringify(Process_::label) + "::a_0_time^A@DM2016"], *this),
        _alpha_1_time_a(p[stringify(Process_::label) + "::a_1_time^A@DM2016"], *this),
        _alpha_2_time_a(p[stringify(Process_::label) + "::a_2_time^A@DM2016"], *this),

        // long, V
        _alpha_0_long_v(p[stringify(Process_::label) + "::a_0_long^V@DM2016"], *this),
        _alpha_1_long_v(p[stringify(Process_::label) + "::a_1_long^V@DM2016"], *this),
        _alpha_2_long_v(p[stringify(Process_::label) + "::a_2_long^V@DM2016"], *this),
        // long, A
        _alpha_0_long_a(p[stringify(Process_::label) + "::a_0_long^A@DM2016"], *this),
        _alpha_1_long_a(p[stringify(Process_::label) + "::a_1_long^A@DM2016"], *this),
        _alpha_2_long_a(p[stringify(Process_::label) + "::a_2_long^A@DM2016"], *this),
        // perp, V
        _alpha_0_perp_v(p[stringify(Process_::label) + "::a_0_perp^V@DM2016"], *this),
        _alpha_1_perp_v(p[stringify(Process_::label) + "::a_1_perp^V@DM2016"], *this),
        _alpha_2_perp_v(p[stringify(Process_::label) + "::a_2_perp^V@DM2016"], *this),
        // perp, A
        _alpha_1_perp_a(p[stringify(Process_::label) + "::a_1_perp^A@DM2016"], *this),
        _alpha_2_perp_a(p[stringify(Process_::label) + "::a_2_perp^A@DM2016"], *this),

        // long, T
        _alpha_0_long_t(p[stringify(Process_::label) + "::a_0_long^T@DM2016"], *this),
        _alpha_1_long_t(p[stringify(Process_::label) + "::a_1_long^T@DM2016"], *this),
        _alpha_2_long_t(p[stringify(Process_::label) + "::a_2_long^T@DM2016"], *this),
        // long, T5
        _alpha_0_long_t5(p[stringify(Process_::label) + "::a_0_long^T5@DM2016"], *this),
        _alpha_1_long_t5(p[stringify(Process_::label) + "::a_1_long^T5@DM2016"], *this),
        _alpha_2_long_t5(p[stringify(Process_::label) + "::a_2_long^T5@DM2016"], *this),
        // perp, T
        _alpha_0_perp_t(p[stringify(Process_::label) + "::a_0_perp^T@DM2016"], *this),
        _alpha_1_perp_t(p[stringify(Process_::label) + "::a_1_perp^T@DM2016"], *this),
        _alpha_2_perp_t(p[stringify(Process_::label) + "::a_2_perp^T@DM2016"], *this),
        // perp, T5
        _alpha_1_perp_t5(p[stringify(Process_::label) + "::a_1_perp^T5@DM2016"], *this),
        _alpha_2_perp_t5(p[stringify(Process_::label) + "::a_2_perp^T5@DM2016"], *this),
        // traits
        _traits(new DM2016FormFactorTraits<Process_>(p))
    {
    }

    template <typename Process_>
    FormFactors<OneHalfPlusToOneHalfPlus> *
    DM2016FormFactors<Process_>::make(const Parameters & parameters, const Options & options)
    {
        return new DM2016FormFactors(parameters, options);
    }

    // vector current
    template <typename Process_>
    double
    DM2016FormFactors<Process_>::f_time_v(const double & s) const
    {
        static const double mR2 = power_of<2>(_traits->m_R_0p);

        const double z = _traits->calc_z(s), z2 = z * z;

        return 1.0 / (1.0 - s / mR2) * (_alpha_0_time_v() + _alpha_1_time_v() * z + _alpha_2_time_v() * z2);
    }

    template <typename Process_>
    double
    DM2016FormFactors<Process_>::f_long_v(const double & s) const
    {
        static const double mR2 = power_of<2>(_traits->m_R_1m);

        const double z = _traits->calc_z(s), z2 = z * z;

        return 1.0 / (1.0 - s / mR2) * (_alpha_0_long_v() + _alpha_1_long_v() * z + _alpha_2_long_v() * z2);
    }

    template <typename Process_>
    double
    DM2016FormFactors<Process_>::f_perp_v(const double & s) const
    {
        static const double mR2 = power_of<2>(_traits->m_R_1m);

        const double z = _traits->calc_z(s), z2 = z * z;

        return 1.0 / (1.0 - s / mR2) * (_alpha_0_perp_v() + _alpha_1_perp_v() * z + _alpha_2_perp_v() * z2);
    }

    // axial vector current
    template <typename Process_>
    double
    DM2016FormFactors<Process_>::f_time_a(const double & s) const
    {
        static const double mR2 = power_of<2>(_traits->m_R_0m);

        const double z = _traits->calc_z(s), z2 = z * z;

        return 1.0 / (1.0 - s / mR2) * (_alpha_0_time_a() + _alpha_1_time_a() * z + _alpha_2_time_a() * z2);
    }

    template <typename Process_>
    double
    DM2016FormFactors<Process_>::f_long_a(const double & s) const
    {
        static const double mR2 = power_of<2>(_traits->m_R_1p);

        const double z = _traits->calc_z(s), z2 = z * z;

        return 1.0 / (1.0 - s / mR2) * (_alpha_0_long_a() + _alpha_1_long_a() * z + _alpha_2_long_a() * z2);
    }

    template <typename Process_>
    double
    DM2016FormFactors<Process_>::f_perp_a(const double & s) const
    {
        static const double mR2 = power_of<2>(_traits->m_R_1p);

        const double z = _traits->calc_z(s), z2 = z * z;

        // Using alpha_0_long_a instead of alpha_0_perp_a, in order to
        // fulfill relation eq. (7), [DM:2016A], p. 3.
        return 1.0 / (1.0 - s / mR2) * (_alpha_0_long_a() + _alpha_1_perp_a() * z + _alpha_2_perp_a() * z2);
    }

    // tensor current
    template <typename Process_>
    double
    DM2016FormFactors<Process_>::f_long_t(const double & s) const
    {
        static const double mR2 = power_of<2>(_traits->m_R_1m);

        const double z = _traits->calc_z(s), z2 = z * z;

        return 1.0 / (1.0 - s / mR2) * (_alpha_0_long_t() + _alpha_1_long_t() * z + _alpha_2_long_t() * z2);
    }

    template <typename Process_>
    double
    DM2016FormFactors<Process_>::f_perp_t(const double & s) const
    {
        static const double mR2 = power_of<2>(_traits->m_R_1m);

        const double z = _traits->calc_z(s), z2 = z * z;

        return 1.0 / (1.0 - s / mR2) * (_alpha_0_perp_t() + _alpha_1_perp_t() * z + _alpha_2_perp_t() * z2);
    }

    // axial tensor current
    template <typename Process_>
    double
    DM2016FormFactors<Process_>::f_long_t5(const double & s) const
    {
        static const double mR2 = power_of<2>(_traits->m_R_1p);

        const double z = _traits->calc_z(s), z2 = z * z;

        return 1.0 / (1.0 - s / mR2) * (_alpha_0_long_t5() + _alpha_1_long_t5() * z + _alpha_2_long_t5() * z2);
    }

    template <typename Process_>
    double
    DM2016FormFactors<Process_>::f_perp_t5(const double & s) const
    {
        static const double mR2 = power_of<2>(_traits->m_R_1p);

        const double z = _traits->calc_z(s), z2 = z * z;

        // Using alpha_0_long_t5 instead of alpha_0_perp_t5, in order to
        // fulfill relation eq. (8), [DM:2016A], p. 3.
        return 1.0 / (1.0 - s / mR2) * (_alpha_0_long_t5() + _alpha_1_perp_t5() * z + _alpha_2_perp_t5() * z2);
    }

    template <typename Process_>
    const std::set<ReferenceName> DM2016FormFactors<Process_>::references
    {
        "DM:2016A"_rn
    };

    template <typename Process_>
    const std::vector<OptionSpecification> DM2016FormFactors<Process_>::options
    {
    };

    template <typename Process_>
    std::vector<OptionSpecification>::const_iterator
    DM2016FormFactors<Process_>::begin_options()
    {
        return options.cbegin();
    }

    template <typename Process_>
    std::vector<OptionSpecification>::const_iterator
    DM2016FormFactors<Process_>::end_options()
    {
        return options.cend();
    }
}

#endif
