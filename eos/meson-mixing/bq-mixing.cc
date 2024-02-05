/* vim: set sw=4 sts=4 et foldmethod=syntax : */

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

#include <eos/maths/power-of.hh>
#include <eos/maths/matrix.hh>
#include <eos/meson-mixing/bq-mixing.hh>
#include <eos/models/model.hh>
#include <eos/utils/destringify.hh>
#include <eos/utils/options-impl.hh>
#include <eos/utils/private_implementation_pattern-impl.hh>

namespace eos
{
    using std::norm;

    template <>
    struct Implementation<BMixing>
    {
        std::shared_ptr<Model> model;

        UsedParameter mu_sbsb;
        UsedParameter mu_sb;

        UsedParameter hbar;

        UsedParameter g_fermi;

        QuarkFlavorOption opt_q;

        UsedParameter m_B;
        UsedParameter m_b;
        UsedParameter m_q;

        UsedParameter f_B;

        UsedParameter tau_B;

        UsedParameter R_1;
        UsedParameter R_2;
        UsedParameter R_3;
        UsedParameter R_4;
        UsedParameter R_5;

        static const std::vector<OptionSpecification> options;

        Implementation(const Parameters & p, const Options & o, ParameterUser & u) :
            model(Model::make(o.get("model", "SM"), p, o)),
            mu_sbsb(p["sbsb::mu"], u),
            mu_sb(p["sb::mu"], u),
            hbar(p["QM::hbar"], u),
            g_fermi(p["WET::G_Fermi"], u),
            opt_q(o, options, "q"),
            m_B(p["mass::B_" + opt_q.str()], u),
            m_b(p["mass::b(MSbar)"], u),
            m_q(p["mass::" + opt_q.str() + "(2GeV)"], u),
            f_B(p["decay-constant::B_" + opt_q.str()], u),
            tau_B(p["life_time::B_" + opt_q.str()], u),
            R_1(p["B_" + opt_q.str() + "<->Bbar_" + opt_q.str() + "::R^1"], u),
            R_2(p["B_" + opt_q.str() + "<->Bbar_" + opt_q.str() + "::R^2"], u),
            R_3(p["B_" + opt_q.str() + "<->Bbar_" + opt_q.str() + "::R^3"], u),
            R_4(p["B_" + opt_q.str() + "<->Bbar_" + opt_q.str() + "::R^4"], u),
            R_5(p["B_" + opt_q.str() + "<->Bbar_" + opt_q.str() + "::R^5"], u)
        {
            u.uses(*model);
        }

        complex<double> M_12() const
        {
            const auto wc = model->wet_sbsb();

            // cf. [DDHLMSW:2019A]
            // TODO: still needs to be evolved to scale mu from reference scale 4.2 GeV.
            const std::array<complex<double>, 8> contributions{{
                wc.c1()  * R_1(),
                wc.c2()  * R_2(),
                wc.c3()  * R_3(),
                wc.c4()  * R_4(),
                wc.c5()  * R_5(),
                wc.c1p() * R_1(), // primed operators share the hadronic matrix elements of their unprimed partners
                wc.c2p() * R_2(),
                wc.c3p() * R_3(),
            }};

            complex<double> result = 0.0;
            for (auto & c : contributions)
            {
                result += c;
            }

            // cf. [BBL:1995A], eq. (XVIII.17), p. 153
            return 4.0 * g_fermi / std::sqrt(2.0) * power_of<2>(model->ckm_tb() * std::conj(model->ckm_ts()))
                * f_B() * f_B() * m_B() / 2.0 * result;

        }

        // cf. [LMPR:2022A], eqs. (2.54) to (2.59)
        // provided by A. Rusov as a Mathematica expression
        std::array<std::array<complex<double>, 20u>, 20u>
        M_cc(const double & mu, const double & sqrtrho) const
        {
            const double rho = sqrtrho * sqrtrho;

            // Matrix elements of operators [qbar Gamma b] [qbar Gamma b] as defined in eqs. (2.40)-(2.42) of [LMPR:2022A]
            // Within HQET, there are only 6 independent operators at leading order in 1/m_b, and only 4 matrix elements
            // after accounting for the parity symmetry of QCD.
            // The bag parameters are known from lattice QCD, and sum rules, and we take the averages
            // from appendix D of GSSS:2022A
            const double B1 = 0.886; // B_1^GSSS
            const double B2 = 0.788; // B_2^GSSS
            const double B3 = 0.926; // B_5^GSSS
            const double B4 = 0.979; // B_4^GSSS
            const double mass_ratio_sq = power_of<2>(m_B / (m_b + m_q));

            const double me  =  f_B * f_B * m_B * m_B;
            const double me1 =  8/3 * B1 *                               me;
            const double me2 = -5/3 * B2 *                               me;
            const double me3 =   -1 * B3 * (  2 + 4/3 * mass_ratio_sq) * me;
            const double me4 =    1 * B4 * (1/3 +   2 * mass_ratio_sq) * me;

            std::array<std::array<complex<double>, 20u>, 20u> result
            {{
                {-(me1*(1 - 4*rho)) + 2*me2*(1 + 2*rho),-(me1*(1 - 4*rho)) + 2*(-me1/2. - me2)*(1 + 2*rho),-3*(me1 + 2*me2)*sqrtrho,-3*(me1 + 2*(-me1/2. - me2))*sqrtrho,3*me1*rho,3*me1*rho,(3*me1*sqrtrho)/2.,(3*me1*sqrtrho)/2.,6*(me1 + 4*me2)*sqrtrho,6*(me1 + 4*(-me1/2. - me2))*sqrtrho,12*me4*rho,-6*me3*rho,3*me3*sqrtrho,-6*me4*sqrtrho,(-4*me4*(1 - rho) - me3*(1 + 2*rho))/2.,(2*me3*(1 - rho) + 2*me4*(1 + 2*rho))/2.,(-3*(me3 + 2*me4)*sqrtrho)/2.,(-3*(-me3 - 2*me4)*sqrtrho)/2.,6*(-me3 + 2*me4)*sqrtrho,6*(-me3 + 2*me4)*sqrtrho},
                {-(me1*(1 - 4*rho)) + 2*(-me1/2. - me2)*(1 + 2*rho),3*(-(me1*(1 - 4*rho)) + 2*(-me1/2. - me2)*(1 + 2*rho)),-3*(me1 + 2*(-me1/2. - me2))*sqrtrho,-9*(me1 + 2*(-me1/2. - me2))*sqrtrho,3*me1*rho,9*me1*rho,(3*me1*sqrtrho)/2.,(9*me1*sqrtrho)/2.,6*(me1 + 4*(-me1/2. - me2))*sqrtrho,18*(me1 + 4*(-me1/2. - me2))*sqrtrho,-6*me3*rho,-18*me3*rho,-6*me4*sqrtrho,-18*me4*sqrtrho,(2*me3*(1 - rho) + 2*me4*(1 + 2*rho))/2.,(3*(2*me3*(1 - rho) + 2*me4*(1 + 2*rho)))/2.,(-3*(-me3 - 2*me4)*sqrtrho)/2.,(-9*(-me3 - 2*me4)*sqrtrho)/2.,6*(-me3 + 2*me4)*sqrtrho,18*(-me3 + 2*me4)*sqrtrho},
                {-3*(me1 + 2*me2)*sqrtrho,-3*(me1 + 2*(-me1/2. - me2))*sqrtrho,12*(me1 + 2*me2)*rho,12*(me1 + 2*(-me1/2. - me2))*rho,(-3*(me1 + 2*me2)*sqrtrho)/2.,(-3*(me1 + 2*(-me1/2. - me2))*sqrtrho)/2.,(-3*(me1 + 2*me2)*(1 - 2*rho))/2.,(-3*(me1 + 2*(-me1/2. - me2))*(1 - 2*rho))/2.,-18*(me1 + 2*me2)*(1 - 2*rho),-18*(me1 + 2*(-me1/2. - me2))*(1 - 2*rho),3*me3*sqrtrho,-6*me4*sqrtrho,-6*me3*(1 - 2*rho),12*me4*(1 - 2*rho),(3*me3*sqrtrho)/2.,-3*me4*sqrtrho,3*me3*rho,-6*me4*rho,36*me3*rho,-72*me4*rho},
                {-3*(me1 + 2*(-me1/2. - me2))*sqrtrho,-9*(me1 + 2*(-me1/2. - me2))*sqrtrho,12*(me1 + 2*(-me1/2. - me2))*rho,36*(me1 + 2*(-me1/2. - me2))*rho,(-3*(me1 + 2*(-me1/2. - me2))*sqrtrho)/2.,(-9*(me1 + 2*(-me1/2. - me2))*sqrtrho)/2.,(-3*(me1 + 2*(-me1/2. - me2))*(1 - 2*rho))/2.,(-9*(me1 + 2*(-me1/2. - me2))*(1 - 2*rho))/2.,-18*(me1 + 2*(-me1/2. - me2))*(1 - 2*rho),-54*(me1 + 2*(-me1/2. - me2))*(1 - 2*rho),-6*me4*sqrtrho,-18*me4*sqrtrho,12*me4*(1 - 2*rho),36*me4*(1 - 2*rho),-3*me4*sqrtrho,-9*me4*sqrtrho,-6*me4*rho,-18*me4*rho,-72*me4*rho,-216*me4*rho},
                {3*me1*rho,3*me1*rho,(-3*(me1 + 2*me2)*sqrtrho)/2.,(-3*(me1 + 2*(-me1/2. - me2))*sqrtrho)/2.,(-(me1*(1 - 4*rho)) + 2*me2*(1 + 2*rho))/4.,(-(me1*(1 - 4*rho)) + 2*(-me1/2. - me2)*(1 + 2*rho))/4.,(3*me2*sqrtrho)/2.,(3*(-me1/2. - me2)*sqrtrho)/2.,6*(me1 + me2)*sqrtrho,6*(me1/2. - me2)*sqrtrho,(-4*me4*(1 - rho) - me3*(1 + 2*rho))/2.,(2*me3*(1 - rho) + 2*me4*(1 + 2*rho))/2.,(3*me3*sqrtrho)/2.,-3*me4*sqrtrho,3*me4*rho,(-3*me3*rho)/2.,(3*me4*sqrtrho)/2.,(-3*me3*sqrtrho)/4.,-6*(me3 + me4)*sqrtrho,-6*(-me3/2. - 2*me4)*sqrtrho},
                {3*me1*rho,9*me1*rho,(-3*(me1 + 2*(-me1/2. - me2))*sqrtrho)/2.,(-9*(me1 + 2*(-me1/2. - me2))*sqrtrho)/2.,(-(me1*(1 - 4*rho)) + 2*(-me1/2. - me2)*(1 + 2*rho))/4.,(3*(-(me1*(1 - 4*rho)) + 2*(-me1/2. - me2)*(1 + 2*rho)))/4.,(3*(-me1/2. - me2)*sqrtrho)/2.,(9*(-me1/2. - me2)*sqrtrho)/2.,6*(me1/2. - me2)*sqrtrho,18*(me1/2. - me2)*sqrtrho,(2*me3*(1 - rho) + 2*me4*(1 + 2*rho))/2.,(3*(2*me3*(1 - rho) + 2*me4*(1 + 2*rho)))/2.,-3*me4*sqrtrho,-9*me4*sqrtrho,(-3*me3*rho)/2.,(-9*me3*rho)/2.,(-3*me3*sqrtrho)/4.,(-9*me3*sqrtrho)/4.,-6*(-me3/2. - 2*me4)*sqrtrho,-18*(-me3/2. - 2*me4)*sqrtrho},
                {(3*me1*sqrtrho)/2.,(3*me1*sqrtrho)/2.,(-3*(me1 + 2*me2)*(1 - 2*rho))/2.,(-3*(me1 + 2*(-me1/2. - me2))*(1 - 2*rho))/2.,(3*me2*sqrtrho)/2.,(3*(-me1/2. - me2)*sqrtrho)/2.,3*me2*rho,3*(-me1/2. - me2)*rho,12*(me1 + me2)*rho,12*(me1/2. - me2)*rho,(-3*(me3 + 2*me4)*sqrtrho)/2.,(-3*(-me3 - 2*me4)*sqrtrho)/2.,3*me3*rho,-6*me4*rho,(3*me4*sqrtrho)/2.,(-3*me3*sqrtrho)/4.,(-(me3*(1 - 4*rho)) + 2*me4*(1 + 2*rho))/4.,(2*me4*(1 - 4*rho) - me3*(1 + 2*rho))/4.,-(me3*(5 - 8*rho)) - 2*me4*(1 + 2*rho),2*me4*(5 - 8*rho) + me3*(1 + 2*rho)},
                {(3*me1*sqrtrho)/2.,(9*me1*sqrtrho)/2.,(-3*(me1 + 2*(-me1/2. - me2))*(1 - 2*rho))/2.,(-9*(me1 + 2*(-me1/2. - me2))*(1 - 2*rho))/2.,(3*(-me1/2. - me2)*sqrtrho)/2.,(9*(-me1/2. - me2)*sqrtrho)/2.,3*(-me1/2. - me2)*rho,9*(-me1/2. - me2)*rho,12*(me1/2. - me2)*rho,36*(me1/2. - me2)*rho,(-3*(-me3 - 2*me4)*sqrtrho)/2.,(-9*(-me3 - 2*me4)*sqrtrho)/2.,-6*me4*rho,-18*me4*rho,(-3*me3*sqrtrho)/4.,(-9*me3*sqrtrho)/4.,(2*me4*(1 - 4*rho) - me3*(1 + 2*rho))/4.,(3*(2*me4*(1 - 4*rho) - me3*(1 + 2*rho)))/4.,2*me4*(5 - 8*rho) + me3*(1 + 2*rho),3*(2*me4*(5 - 8*rho) + me3*(1 + 2*rho))},
                {6*(me1 + 4*me2)*sqrtrho,6*(me1 + 4*(-me1/2. - me2))*sqrtrho,-18*(me1 + 2*me2)*(1 - 2*rho),-18*(me1 + 2*(-me1/2. - me2))*(1 - 2*rho),6*(me1 + me2)*sqrtrho,6*(me1/2. - me2)*sqrtrho,12*(me1 + me2)*rho,12*(me1/2. - me2)*rho,48*(2*me1 + 5*me2)*rho,48*(2*me1 + 5*(-me1/2. - me2))*rho,6*(-me3 + 2*me4)*sqrtrho,6*(-me3 + 2*me4)*sqrtrho,36*me3*rho,-72*me4*rho,-6*(me3 + me4)*sqrtrho,-6*(-me3/2. - 2*me4)*sqrtrho,-(me3*(5 - 8*rho)) - 2*me4*(1 + 2*rho),2*me4*(5 - 8*rho) + me3*(1 + 2*rho),-4*(me3*(13 - 28*rho) - 2*me4*(1 + 2*rho)),-4*(-2*me4*(13 - 28*rho) + me3*(1 + 2*rho))},
                {6*(me1 + 4*(-me1/2. - me2))*sqrtrho,18*(me1 + 4*(-me1/2. - me2))*sqrtrho,-18*(me1 + 2*(-me1/2. - me2))*(1 - 2*rho),-54*(me1 + 2*(-me1/2. - me2))*(1 - 2*rho),6*(me1/2. - me2)*sqrtrho,18*(me1/2. - me2)*sqrtrho,12*(me1/2. - me2)*rho,36*(me1/2. - me2)*rho,48*(2*me1 + 5*(-me1/2. - me2))*rho,144*(2*me1 + 5*(-me1/2. - me2))*rho,6*(-me3 + 2*me4)*sqrtrho,18*(-me3 + 2*me4)*sqrtrho,-72*me4*rho,-216*me4*rho,-6*(-me3/2. - 2*me4)*sqrtrho,-18*(-me3/2. - 2*me4)*sqrtrho,2*me4*(5 - 8*rho) + me3*(1 + 2*rho),3*(2*me4*(5 - 8*rho) + me3*(1 + 2*rho)),-4*(-2*me4*(13 - 28*rho) + me3*(1 + 2*rho)),-12*(-2*me4*(13 - 28*rho) + me3*(1 + 2*rho))},
                {12*me4*rho,-6*me3*rho,3*me3*sqrtrho,-6*me4*sqrtrho,(-4*me4*(1 - rho) - me3*(1 + 2*rho))/2.,(2*me3*(1 - rho) + 2*me4*(1 + 2*rho))/2.,(-3*(me3 + 2*me4)*sqrtrho)/2.,(-3*(-me3 - 2*me4)*sqrtrho)/2.,6*(-me3 + 2*me4)*sqrtrho,6*(-me3 + 2*me4)*sqrtrho,-(me1*(1 - 4*rho)) + 2*me2*(1 + 2*rho),-(me1*(1 - 4*rho)) + 2*(-me1/2. - me2)*(1 + 2*rho),-3*(me1 + 2*me2)*sqrtrho,-3*(me1 + 2*(-me1/2. - me2))*sqrtrho,3*me1*rho,3*me1*rho,(3*me1*sqrtrho)/2.,(3*me1*sqrtrho)/2.,6*(me1 + 4*me2)*sqrtrho,6*(me1 + 4*(-me1/2. - me2))*sqrtrho},
                {-6*me3*rho,-18*me3*rho,-6*me4*sqrtrho,-18*me4*sqrtrho,(2*me3*(1 - rho) + 2*me4*(1 + 2*rho))/2.,(3*(2*me3*(1 - rho) + 2*me4*(1 + 2*rho)))/2.,(-3*(-me3 - 2*me4)*sqrtrho)/2.,(-9*(-me3 - 2*me4)*sqrtrho)/2.,6*(-me3 + 2*me4)*sqrtrho,18*(-me3 + 2*me4)*sqrtrho,-(me1*(1 - 4*rho)) + 2*(-me1/2. - me2)*(1 + 2*rho),3*(-(me1*(1 - 4*rho)) + 2*(-me1/2. - me2)*(1 + 2*rho)),-3*(me1 + 2*(-me1/2. - me2))*sqrtrho,-9*(me1 + 2*(-me1/2. - me2))*sqrtrho,3*me1*rho,9*me1*rho,(3*me1*sqrtrho)/2.,(9*me1*sqrtrho)/2.,6*(me1 + 4*(-me1/2. - me2))*sqrtrho,18*(me1 + 4*(-me1/2. - me2))*sqrtrho},
                {3*me3*sqrtrho,-6*me4*sqrtrho,-6*me3*(1 - 2*rho),12*me4*(1 - 2*rho),(3*me3*sqrtrho)/2.,-3*me4*sqrtrho,3*me3*rho,-6*me4*rho,36*me3*rho,-72*me4*rho,-3*(me1 + 2*me2)*sqrtrho,-3*(me1 + 2*(-me1/2. - me2))*sqrtrho,12*(me1 + 2*me2)*rho,12*(me1 + 2*(-me1/2. - me2))*rho,(-3*(me1 + 2*me2)*sqrtrho)/2.,(-3*(me1 + 2*(-me1/2. - me2))*sqrtrho)/2.,(-3*(me1 + 2*me2)*(1 - 2*rho))/2.,(-3*(me1 + 2*(-me1/2. - me2))*(1 - 2*rho))/2.,-18*(me1 + 2*me2)*(1 - 2*rho),-18*(me1 + 2*(-me1/2. - me2))*(1 - 2*rho)},
                {-6*me4*sqrtrho,-18*me4*sqrtrho,12*me4*(1 - 2*rho),36*me4*(1 - 2*rho),-3*me4*sqrtrho,-9*me4*sqrtrho,-6*me4*rho,-18*me4*rho,-72*me4*rho,-216*me4*rho,-3*(me1 + 2*(-me1/2. - me2))*sqrtrho,-9*(me1 + 2*(-me1/2. - me2))*sqrtrho,12*(me1 + 2*(-me1/2. - me2))*rho,36*(me1 + 2*(-me1/2. - me2))*rho,(-3*(me1 + 2*(-me1/2. - me2))*sqrtrho)/2.,(-9*(me1 + 2*(-me1/2. - me2))*sqrtrho)/2.,(-3*(me1 + 2*(-me1/2. - me2))*(1 - 2*rho))/2.,(-9*(me1 + 2*(-me1/2. - me2))*(1 - 2*rho))/2.,-18*(me1 + 2*(-me1/2. - me2))*(1 - 2*rho),-54*(me1 + 2*(-me1/2. - me2))*(1 - 2*rho)},
                {(-4*me4*(1 - rho) - me3*(1 + 2*rho))/2.,(2*me3*(1 - rho) + 2*me4*(1 + 2*rho))/2.,(3*me3*sqrtrho)/2.,-3*me4*sqrtrho,3*me4*rho,(-3*me3*rho)/2.,(3*me4*sqrtrho)/2.,(-3*me3*sqrtrho)/4.,-6*(me3 + me4)*sqrtrho,-6*(-me3/2. - 2*me4)*sqrtrho,3*me1*rho,3*me1*rho,(-3*(me1 + 2*me2)*sqrtrho)/2.,(-3*(me1 + 2*(-me1/2. - me2))*sqrtrho)/2.,(-(me1*(1 - 4*rho)) + 2*me2*(1 + 2*rho))/4.,(-(me1*(1 - 4*rho)) + 2*(-me1/2. - me2)*(1 + 2*rho))/4.,(3*me2*sqrtrho)/2.,(3*(-me1/2. - me2)*sqrtrho)/2.,6*(me1 + me2)*sqrtrho,6*(me1/2. - me2)*sqrtrho},
                {(2*me3*(1 - rho) + 2*me4*(1 + 2*rho))/2.,(3*(2*me3*(1 - rho) + 2*me4*(1 + 2*rho)))/2.,-3*me4*sqrtrho,-9*me4*sqrtrho,(-3*me3*rho)/2.,(-9*me3*rho)/2.,(-3*me3*sqrtrho)/4.,(-9*me3*sqrtrho)/4.,-6*(-me3/2. - 2*me4)*sqrtrho,-18*(-me3/2. - 2*me4)*sqrtrho,3*me1*rho,9*me1*rho,(-3*(me1 + 2*(-me1/2. - me2))*sqrtrho)/2.,(-9*(me1 + 2*(-me1/2. - me2))*sqrtrho)/2.,(-(me1*(1 - 4*rho)) + 2*(-me1/2. - me2)*(1 + 2*rho))/4.,(3*(-(me1*(1 - 4*rho)) + 2*(-me1/2. - me2)*(1 + 2*rho)))/4.,(3*(-me1/2. - me2)*sqrtrho)/2.,(9*(-me1/2. - me2)*sqrtrho)/2.,6*(me1/2. - me2)*sqrtrho,18*(me1/2. - me2)*sqrtrho},
                {(-3*(me3 + 2*me4)*sqrtrho)/2.,(-3*(-me3 - 2*me4)*sqrtrho)/2.,3*me3*rho,-6*me4*rho,(3*me4*sqrtrho)/2.,(-3*me3*sqrtrho)/4.,(-(me3*(1 - 4*rho)) + 2*me4*(1 + 2*rho))/4.,(2*me4*(1 - 4*rho) - me3*(1 + 2*rho))/4.,-(me3*(5 - 8*rho)) - 2*me4*(1 + 2*rho),2*me4*(5 - 8*rho) + me3*(1 + 2*rho),(3*me1*sqrtrho)/2.,(3*me1*sqrtrho)/2.,(-3*(me1 + 2*me2)*(1 - 2*rho))/2.,(-3*(me1 + 2*(-me1/2. - me2))*(1 - 2*rho))/2.,(3*me2*sqrtrho)/2.,(3*(-me1/2. - me2)*sqrtrho)/2.,3*me2*rho,3*(-me1/2. - me2)*rho,12*(me1 + me2)*rho,12*(me1/2. - me2)*rho},
                {(-3*(-me3 - 2*me4)*sqrtrho)/2.,(-9*(-me3 - 2*me4)*sqrtrho)/2.,-6*me4*rho,-18*me4*rho,(-3*me3*sqrtrho)/4.,(-9*me3*sqrtrho)/4.,(2*me4*(1 - 4*rho) - me3*(1 + 2*rho))/4.,(3*(2*me4*(1 - 4*rho) - me3*(1 + 2*rho)))/4.,2*me4*(5 - 8*rho) + me3*(1 + 2*rho),3*(2*me4*(5 - 8*rho) + me3*(1 + 2*rho)),(3*me1*sqrtrho)/2.,(9*me1*sqrtrho)/2.,(-3*(me1 + 2*(-me1/2. - me2))*(1 - 2*rho))/2.,(-9*(me1 + 2*(-me1/2. - me2))*(1 - 2*rho))/2.,(3*(-me1/2. - me2)*sqrtrho)/2.,(9*(-me1/2. - me2)*sqrtrho)/2.,3*(-me1/2. - me2)*rho,9*(-me1/2. - me2)*rho,12*(me1/2. - me2)*rho,36*(me1/2. - me2)*rho},
                {6*(-me3 + 2*me4)*sqrtrho,6*(-me3 + 2*me4)*sqrtrho,36*me3*rho,-72*me4*rho,-6*(me3 + me4)*sqrtrho,-6*(-me3/2. - 2*me4)*sqrtrho,-(me3*(5 - 8*rho)) - 2*me4*(1 + 2*rho),2*me4*(5 - 8*rho) + me3*(1 + 2*rho),-4*(me3*(13 - 28*rho) - 2*me4*(1 + 2*rho)),-4*(-2*me4*(13 - 28*rho) + me3*(1 + 2*rho)),6*(me1 + 4*me2)*sqrtrho,6*(me1 + 4*(-me1/2. - me2))*sqrtrho,-18*(me1 + 2*me2)*(1 - 2*rho),-18*(me1 + 2*(-me1/2. - me2))*(1 - 2*rho),6*(me1 + me2)*sqrtrho,6*(me1/2. - me2)*sqrtrho,12*(me1 + me2)*rho,12*(me1/2. - me2)*rho,48*(2*me1 + 5*me2)*rho,48*(2*me1 + 5*(-me1/2. - me2))*rho},
                {6*(-me3 + 2*me4)*sqrtrho,18*(-me3 + 2*me4)*sqrtrho,-72*me4*rho,-216*me4*rho,-6*(-me3/2. - 2*me4)*sqrtrho,-18*(-me3/2. - 2*me4)*sqrtrho,2*me4*(5 - 8*rho) + me3*(1 + 2*rho),3*(2*me4*(5 - 8*rho) + me3*(1 + 2*rho)),-4*(-2*me4*(13 - 28*rho) + me3*(1 + 2*rho)),-12*(-2*me4*(13 - 28*rho) + me3*(1 + 2*rho)),6*(me1 + 4*(-me1/2. - me2))*sqrtrho,18*(me1 + 4*(-me1/2. - me2))*sqrtrho,-18*(me1 + 2*(-me1/2. - me2))*(1 - 2*rho),-54*(me1 + 2*(-me1/2. - me2))*(1 - 2*rho),6*(me1/2. - me2)*sqrtrho,18*(me1/2. - me2)*sqrtrho,12*(me1/2. - me2)*rho,36*(me1/2. - me2)*rho,48*(2*me1 + 5*(-me1/2. - me2))*rho,144*(2*me1 + 5*(-me1/2. - me2))*rho}
            }};

            return result;
        }

        complex<double> G_12_cc() const
        {
            const auto wc = model->wet_sbcc();

            const double mu      = mu_sb();
            const double m_b     = model->m_b_msbar(mu);
            const double sqrtrho = model->m_c_msbar(mu) / m_b;
            const double rho     = sqrtrho * sqrtrho;
            const auto   M       = M_cc(mu, sqrtrho);

            // transforming the Wilson coefficients from the EOS basis
            // to the basis used in [LMPR:2022A], eqs. (2.1) to (2.6) (with the replacement u->c, d->s)
            const std::array<std::complex<double>, 20u> C
            {
                wc.c2() / 2.0 + 8.0 * wc.c4(),
                wc.c1() - wc.c2() / 6.0 + 16.0 * wc.c3() - (8.0 * wc.c4()) / 3.0,
                -4.0 * wc.c10() - wc.c6() / 4.0,
                1.0 / 12.0 * (16.0 * wc.c10() - 6.0 * wc.c5() + wc.c6() - 96.0 * wc.c9()),
                -wc.c2() - 4.0 * wc.c4(),
                1.0 / 3.0 * (-6.0 * wc.c1() + wc.c2() + 4.0 * (wc.c4() - 6.0 * wc.c3())),
                32.0 * wc.c10() - wc.c6() / 4.0 - 3.0 * wc.c8(),
                -((32.0 * wc.c10()) / 3.0) - wc.c5() / 2.0 + wc.c6() / 12.0 - 6.0 * wc.c7() + wc.c8() + 64.0 * wc.c9(),
                -8.0 * wc.c10() - wc.c6() / 16.0 + wc.c8() / 4.0,
                1.0 / 48.0 * (128.0 * wc.c10() - 6.0 * wc.c5() + wc.c6() + 24.0 * wc.c7() - 4.0 * wc.c8() - 768.0 * wc.c9()),
                wc.c2p() / 2.0 + 8.0 * wc.c4p(), wc.c1p() - wc.c2p() / 6.0 + 16.0 * wc.c3p() - (8.0 * wc.c4p()) / 3.0,
                -4.0 * wc.c10p() - wc.c6p() / 4.0, 1.0 / 12.0 * (16.0 * wc.c10p() - 6.0 * wc.c5p() + wc.c6p() - 96.0 * wc.c9p()),
                -wc.c2p() - 4.0 * wc.c4p(),
                1.0 / 3.0 * (-6.0 * wc.c1p() + wc.c2p() + 4.0 * (wc.c4p() - 6.0 * wc.c3p())), 32.0 * wc.c10p() - wc.c6p() / 4.0 - 3.0 * wc.c8p(),
                -((32.0 * wc.c10p()) / 3.0) - wc.c5p() / 2.0 + wc.c6p() / 12.0 - 6.0 * wc.c7p() + wc.c8p() + 64.0 * wc.c9p(),
                -8.0 * wc.c10p() - wc.c6p() / 16.0 + wc.c8p() / 4.0, 1.0 / 48.0 * (128.0 * wc.c10p() - 6.0 * wc.c5p() + wc.c6p() + 24.0 * wc.c7p() - 4.0 * wc.c8p() - 768.0 * wc.c9p())
            };
            std::array<std::complex<double>, 20u> Cconj;
            complex<double> (*conj)(const std::complex<double> &) = &std::conj<double>;
            std::transform(C.cbegin(), C.cend(), Cconj.begin(), conj);

            const double ckm = abs(model->ckm_cs() * model->ckm_cb());

            // Extra factor of 2 mB in the denominator takes account of the matrix element normalisation <O> = <Bbar|O|B> / 2 mB
            return power_of<2>(g_fermi * m_b * ckm * (1.0 - rho)) / (48.0 * m_B * M_PI) * real(dot(Cconj, (M * C))) * 1.0e-12;
        }
    };

    const std::vector<OptionSpecification>
    Implementation<BMixing>::options
    {
    };

    BMixing::BMixing(const Parameters & parameters, const Options & options) :
        PrivateImplementationPattern<BMixing>(new Implementation<BMixing>(parameters, options, *this))
    {
    }

    BMixing::~BMixing()
    {
    }

    double
    BMixing::delta_m() const
    {
        // cf. [BBL:1995A], eq. (XVIII.16), p. 153
        return 2.0 * std::abs(_imp->M_12()) / _imp->hbar() * 1.0e-12; // return value in ps^-1
    }

    double
    BMixing::delta_gamma() const
    {
        // cf. [ABL:2015A], eq. (12), p. 3
        // Note this is just the leading order contributions from bscc operators, not the full SM result
        const auto G_12_cc = _imp->G_12_cc();
        const auto phi_12 = std::arg(-_imp->M_12() / G_12_cc);
        return 2.0 * std::abs(G_12_cc) * std::cos(phi_12) / _imp->hbar() * 1.0e-12; // return value in ps^-1
    }

    const std::set<ReferenceName>
    BMixing::references
    {
    };

    std::vector<OptionSpecification>::const_iterator
    BMixing::begin_options()
    {
        return Implementation<BMixing>::options.cbegin();
    }

    std::vector<OptionSpecification>::const_iterator
    BMixing::end_options()
    {
        return Implementation<BMixing>::options.cend();
    }
}
