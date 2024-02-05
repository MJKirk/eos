/* vim: set sw=4 sts=4 et tw=150 foldmethod=marker : */

/*
 * Copyright (c) 2021 Danny van Dyk
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

#include <eos/observable-impl.hh>
#include <eos/meson-mixing/bq-mixing.hh>
#include <eos/utils/concrete-cacheable-observable.hh>
#include <eos/utils/concrete_observable.hh>

namespace eos
{
    // B_q-Bbar_q mixing
    // {{{
    ObservableGroup
    make_bq_mixing_group()
    {
        auto imp = new Implementation<ObservableGroup>(
            R"(Observables in $B_q$--$\bar{B}_q$ mixing)",
            R"()",
            {
                make_observable("B_q<->Bbar_q::DeltaM", R"(\Delta M_q(B_q\leftrightarrow \bar{B}_q))",
                        Unit::InversePicoSecond(),
                        &BMixing::delta_m),
            }
        );

        return ObservableGroup(imp);
    }

    ObservableGroup
    make_bs_mixing_group()
    {
        auto imp = new Implementation<ObservableGroup>(
            R"(Observables in $B_s$--$\bar{B}_s$ mixing)",
            R"()",
            {
                make_observable("B_s<->Bbar_s::DeltaGamma", R"(\Delta \Gamma_s(B_s\leftrightarrow \bar{B}_s))",
                        Unit::InversePicoSecond(),
                        &BMixing::delta_gamma,
                        std::make_tuple(),
                        Options{ { "q", "s" } }),
            }
        );

        return ObservableGroup(imp);
    }
    // }}}

    ObservableSection
    make_meson_mixing_section()
    {
        auto imp = new Implementation<ObservableSection>(
            "Observables in neutral meson mixing",
            "",
            {
                // Observables implemented for both Bs and Bd
                make_bq_mixing_group(),
                // Observables implemented for just Bs
                make_bs_mixing_group(),
            }
        );

        return ObservableSection(imp);
    }
}
