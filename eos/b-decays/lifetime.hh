/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2023 Danny van Dyk
 * Copyright (c) 2024 Stefan Meiser
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

#ifndef EOS_GUARD_EOS_B_DECAYS_LIFETIME_HH
#define EOS_GUARD_EOS_B_DECAYS_LIFETIME_HH 1

#include <eos/maths/complex.hh>
#include <eos/utils/options.hh>
#include <eos/utils/parameters.hh>
#include <eos/utils/private_implementation_pattern.hh>
#include <eos/utils/reference-name.hh>

namespace eos
{
    class Lifetime :
        public ParameterUser,
        public PrivateImplementationPattern<Lifetime>
    {
        public:
            Lifetime(const Parameters & parameters, const Options & options);
            ~Lifetime();

            // Observables
            // Partial decay width contributions due to dbcu operators
            // in the heavy-quark expansion.
            double decay_width_dbcu_dim3_lo() const;
            double decay_width_dbcu_dim6_lo() const;
            // Partial decay width contributions due to sbcu operators
            // in the heavy-quark expansion
            double decay_width_sbcu_dim3_lo() const;
            double decay_width_sbcu_dim6_lo() const;

            /*!
             * References used in the computation of our observables.
             */
            static const std::set<ReferenceName> references;

            /*!
             * Options used in the computation of our observables.
             */
            static std::vector<OptionSpecification>::const_iterator begin_options();
            static std::vector<OptionSpecification>::const_iterator end_options();
    };
}

#endif
