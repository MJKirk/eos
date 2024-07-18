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
            p["pi->pi::b_(+,1)^0@KKRvD2024"] = -0.4523;
            p["pi->pi::b_(+,1)^1@KKRvD2024"] = -0.03156;
            p["pi->pi::b_(+,1)^2@KKRvD2024"] = 0.007017;
            p["pi->pi::b_(+,1)^3@KKRvD2024"] = 0.06082;
            p["pi->pi::b_(+,1)^4@KKRvD2024"] = 0.05651;
            p["pi->pi::b_(+,1)^5@KKRvD2024"] = 0.05912;
            p["pi->pi::b_(+,1)^6@KKRvD2024"] = 0.04208;
            p["pi->pi::b_(+,1)^7@KKRvD2024"] = 0.03461;
            p["pi->pi::b_(+,1)^8@KKRvD2024"] = 0.008998;
            p["pi->pi::b_(+,1)^9@KKRvD2024"] = 0.01179;
            p["pi->pi::M_(+,1)@KKRvD2024"] = 0.7614;
            p["pi->pi::Gamma_(+,1)@KKRvD2024"] = 0.1475;
            p["pi->pi::Re{c}_(+,1)@KKRvD2024"] = -0.2101;
            p["pi->pi::Im{c}_(+,1)@KKRvD2024"] = -0.01010;

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

                TEST_CHECK_NEARLY_EQUAL(ff.f_p(-1.0), 0.13369753,  eps);
                TEST_CHECK_NEARLY_EQUAL(ff.f_p( 0.0), 0.024614699, eps);
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

                TEST_CHECK_NEARLY_EQUAL(real(ff.f_p( 0.0)),  0.024614699, eps);
                TEST_CHECK_NEARLY_EQUAL(imag(ff.f_p( 0.0)),  0.0,         eps);
                TEST_CHECK_NEARLY_EQUAL(real(ff.f_p(+0.1)),  1.8373454,   eps);
                TEST_CHECK_NEARLY_EQUAL(imag(ff.f_p(+0.1)), -0.43932373,  eps);
                TEST_CHECK_NEARLY_EQUAL(real(ff.f_p(+0.5)),  2.0599461,   eps);
                TEST_CHECK_NEARLY_EQUAL(imag(ff.f_p(+0.5)), -4.7897676,   eps);
            }
        }
} parametric_KKRvD2024_test;
