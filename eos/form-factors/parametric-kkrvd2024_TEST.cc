/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2024 Matthew Kirk
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

#include <test/test.hh>
#include <eos/form-factors/parametric-kkrvd2024.hh>

#include <cmath>
#include <limits>
#include <vector>

using namespace test;
using namespace eos;

class ParametricKKRvD2024Test :
    public TestCase
{
    public:
        ParametricKKRvD2024Test() :
            TestCase("parametric_KKRvD2024_test")
        {
        }

        virtual void run() const
        {
            static const double eps = 1e-7;

            Parameters p = Parameters::Defaults();
            p["mass::pi^+"]                  = 0.13957;
            p["pi->pi::t_0@KKRvD2024"]       = 0.0;
            p["pi->pi::b_(+,1)^1@KKRvD2024"] = 0.01101;
            p["pi->pi::b_(+,1)^2@KKRvD2024"] = 0.01342;
            p["pi->pi::b_(+,1)^3@KKRvD2024"] = 0.01119;
            p["pi->pi::b_(+,1)^4@KKRvD2024"] = 0.02332;
            p["pi->pi::b_(+,1)^5@KKRvD2024"] = 0.004137;
            p["pi->pi::b_(+,1)^6@KKRvD2024"] = 0.01153;
            p["pi->pi::b_(+,1)^7@KKRvD2024"] = 0.001415;
            p["pi->pi::b_(+,1)^8@KKRvD2024"] = 0.009508;
            p["pi->pi::b_(+,1)^9@KKRvD2024"] = 0.0009328;
            p["pi->pi::M_(+,1)@KKRvD2024"] = 0.7616;
            p["pi->pi::Gamma_(+,1)@KKRvD2024"] = 0.1480;
            p["pi->pi::Re{c}_(+,1)@KKRvD2024"] = 0.08049;
            p["pi->pi::Im{c}_(+,1)@KKRvD2024"] = -0.2830;

            /* P->P factory */
            {
                std::shared_ptr<FormFactors<PToP>> ff = FormFactorFactory<PToP>::create("pi->pi::KKRvD2024", p, Options{ });

                TEST_CHECK(nullptr != ff);
            }

            /* f_+ and its auxiliary functions at spacelike and lightlike q2 <= 0.0 */
            {
                KKRvD2024FormFactors<PiToPi> ff(p, Options{ });

                const double chi = 3.52e-3;

                TEST_CHECK_NEARLY_EQUAL(ff.z(-1.0), 0.5762159, eps);
                TEST_CHECK_NEARLY_EQUAL(ff.z( 0.0), 0.0,      eps);

                TEST_CHECK_NEARLY_EQUAL(ff.phi_p(ff.z(-1.0), chi), 1.303305e-1, eps);
                TEST_CHECK_NEARLY_EQUAL(ff.phi_p(ff.z( 0.0), chi), 2.969088e-2, eps);

                TEST_CHECK_NEARLY_EQUAL(ff.f_p(-1.0), 0.77095895,  eps);
                TEST_CHECK_NEARLY_EQUAL(ff.f_p( 0.0), 1.0, eps);
            }

            /* 0->PP factory */
            {
                std::shared_ptr<FormFactors<VacuumToPP>> ff = FormFactorFactory<VacuumToPP>::create("0->pipi::KKRvD2024", p, Options{ });

                TEST_CHECK(nullptr != ff);
            }

            /* f_+ at timelike q2 > 0.0 */
            {
                KKRvD2024FormFactors<VacuumToPiPi> ff(p, Options{ });

                const double chi = 3.52e-3;

                TEST_CHECK_NEARLY_EQUAL(real(ff.z( 0.0)),  0.0,      eps);
                TEST_CHECK_NEARLY_EQUAL(imag(ff.z( 0.0)),  0.0,      eps);
                TEST_CHECK_NEARLY_EQUAL(real(ff.z(+0.1)), -0.5583828, eps);
                TEST_CHECK_NEARLY_EQUAL(imag(ff.z(+0.1)),  0.8295834, eps);
                TEST_CHECK_NEARLY_EQUAL(real(ff.z(+0.5)),  0.6883234, eps);
                TEST_CHECK_NEARLY_EQUAL(imag(ff.z(+0.5)),  0.7254039, eps);

                TEST_CHECK_NEARLY_EQUAL(real(ff.phi_p(ff.z( 0.0), chi)),  0.02969088,  eps);
                TEST_CHECK_NEARLY_EQUAL(imag(ff.phi_p(ff.z( 0.0), chi)),  0.0,        eps);
                TEST_CHECK_NEARLY_EQUAL(real(ff.phi_p(ff.z(+0.1), chi)),  0.00361219, eps);
                TEST_CHECK_NEARLY_EQUAL(imag(ff.phi_p(ff.z(+0.1), chi)),  0.00827876, eps);
                TEST_CHECK_NEARLY_EQUAL(real(ff.phi_p(ff.z(+0.5), chi)), -0.04728994,  eps);
                TEST_CHECK_NEARLY_EQUAL(imag(ff.phi_p(ff.z(+0.5), chi)),  0.06171049,  eps);

                TEST_CHECK_NEARLY_EQUAL(real(ff.f_p( 0.0)),  1.0, eps);
                TEST_CHECK_NEARLY_EQUAL(imag(ff.f_p( 0.0)),  0.0,         eps);
                TEST_CHECK_NEARLY_EQUAL(real(ff.f_p(+0.1)),  0.48037094,   eps);
                TEST_CHECK_NEARLY_EQUAL(imag(ff.f_p(+0.1)), -1.92660137,  eps);
                TEST_CHECK_NEARLY_EQUAL(real(ff.f_p(+0.5)),  4.80894779,   eps);
                TEST_CHECK_NEARLY_EQUAL(imag(ff.f_p(+0.5)),  2.02035446,   eps);
            }
        }
} parametric_KKRvD2024_test;
