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

            // t0 = 0
            {

                Parameters p = Parameters::Defaults();
                p["mass::pi^+"]                  = 0.13957;
                p["pi->pi::t_0@KKRvD2024"]       = 0.0;
                p["pi->pi::b_(+,1)^1@KKRvD2024"] = -0.027048;
                p["pi->pi::b_(+,1)^2@KKRvD2024"] = -0.010680;
                p["pi->pi::b_(+,1)^3@KKRvD2024"] = -0.0070125;
                p["pi->pi::b_(+,1)^4@KKRvD2024"] = 0.0;
                p["pi->pi::b_(+,1)^5@KKRvD2024"] = 0.0;
                p["pi->pi::b_(+,1)^6@KKRvD2024"] = 0.0;
                p["pi->pi::b_(+,1)^7@KKRvD2024"] = 0.0;
                p["pi->pi::b_(+,1)^8@KKRvD2024"] = 0.0;
                p["pi->pi::b_(+,1)^9@KKRvD2024"] = 0.0;
                p["pi->pi::M_(+,1)@KKRvD2024"] = 0.76173;
                p["pi->pi::Gamma_(+,1)@KKRvD2024"] = 0.14620;
                p["pi->pi::Re{c}_(+,1)@KKRvD2024"] = -0.21691;
                p["pi->pi::Im{c}_(+,1)@KKRvD2024"] = 0.0030072;

                /* P->P factory */
                {
                    std::shared_ptr<FormFactors<PToP>> ff = FormFactorFactory<PToP>::create("pi->pi::KKRvD2024", p, Options{ });

                    TEST_CHECK(nullptr != ff);
                }

                /* f_+ and its auxiliary functions at spacelike and lightlike q2 <= 0.0 */
                {
                    KKRvD2024FormFactors<PiToPi> ff(p, Options{ });

                    const double chi = 0.0258815;

                    TEST_CHECK_NEARLY_EQUAL(ff.z(-1.0), 0.5762159, eps);
                    TEST_CHECK_NEARLY_EQUAL(ff.z( 0.0), 0.0,      eps);

                    TEST_CHECK_NEARLY_EQUAL(ff.phi_p(ff.z(-1.0), chi), 0.0818588362, eps);
                    TEST_CHECK_NEARLY_EQUAL(ff.phi_p(ff.z( 0.0), chi), 0.0243975773, eps);

                    TEST_CHECK_NEARLY_EQUAL(ff.f_p(-1.0), 0.168686490,  eps);
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

                    const double chi = 0.0258815;

                    TEST_CHECK_NEARLY_EQUAL(real(ff.z( 0.0)),  0.0,      eps);
                    TEST_CHECK_NEARLY_EQUAL(imag(ff.z( 0.0)),  0.0,      eps);
                    TEST_CHECK_NEARLY_EQUAL(real(ff.z(+0.1)), -0.5583828, eps);
                    TEST_CHECK_NEARLY_EQUAL(imag(ff.z(+0.1)), -0.8295834, eps);
                    TEST_CHECK_NEARLY_EQUAL(real(ff.z(+0.5)),  0.6883234, eps);
                    TEST_CHECK_NEARLY_EQUAL(imag(ff.z(+0.5)), -0.7254039, eps);

                    TEST_CHECK_NEARLY_EQUAL(real(ff.phi_p(ff.z( 0.0), chi)),  0.0243975773,  eps);
                    TEST_CHECK_NEARLY_EQUAL(imag(ff.phi_p(ff.z( 0.0), chi)),  0.0,        eps);
                    TEST_CHECK_NEARLY_EQUAL(real(ff.phi_p(ff.z(+0.1), chi)),  0.0044399384, eps);
                    TEST_CHECK_NEARLY_EQUAL(imag(ff.phi_p(ff.z(+0.1), chi)), -0.0075823849, eps);
                    TEST_CHECK_NEARLY_EQUAL(real(ff.phi_p(ff.z(+0.5), chi)), -0.0150675279,  eps);
                    TEST_CHECK_NEARLY_EQUAL(imag(ff.phi_p(ff.z(+0.5), chi)), -0.0598239508,  eps);

                    TEST_CHECK_NEARLY_EQUAL(real(ff.f_p( 0.0)),  1.0, eps);
                    TEST_CHECK_NEARLY_EQUAL(imag(ff.f_p( 0.0)),  0.0,         eps);
                    TEST_CHECK_NEARLY_EQUAL(real(ff.f_p(+0.1)),  1.000788609,   eps);
                    TEST_CHECK_NEARLY_EQUAL(imag(ff.f_p(+0.1)),  0.328870780,  eps);
                    TEST_CHECK_NEARLY_EQUAL(real(ff.f_p(+0.5)),  2.45351138,   eps);
                    TEST_CHECK_NEARLY_EQUAL(imag(ff.f_p(+0.5)),  4.58789522,   eps);
                }
            }

            // t0 = -2
            {

                Parameters p = Parameters::Defaults();
                p["mass::pi^+"]                  = 0.13957;
                p["pi->pi::t_0@KKRvD2024"]       = -2.0;
                p["pi->pi::b_(+,1)^1@KKRvD2024"] = -0.35095;
                p["pi->pi::b_(+,1)^2@KKRvD2024"] = -0.29473;
                p["pi->pi::b_(+,1)^3@KKRvD2024"] = -0.17205;
                p["pi->pi::b_(+,1)^4@KKRvD2024"] = 0.0;
                p["pi->pi::b_(+,1)^5@KKRvD2024"] = 0.0;
                p["pi->pi::b_(+,1)^6@KKRvD2024"] = 0.0;
                p["pi->pi::b_(+,1)^7@KKRvD2024"] = 0.0;
                p["pi->pi::b_(+,1)^8@KKRvD2024"] = 0.0;
                p["pi->pi::b_(+,1)^9@KKRvD2024"] = 0.0;
                p["pi->pi::M_(+,1)@KKRvD2024"] = 0.76032;
                p["pi->pi::Gamma_(+,1)@KKRvD2024"] = 0.14638;
                p["pi->pi::Re{c}_(+,1)@KKRvD2024"] = 0.32927;
                p["pi->pi::Im{c}_(+,1)@KKRvD2024"] = 0.43654;

                /* P->P factory */
                {
                    std::shared_ptr<FormFactors<PToP>> ff = FormFactorFactory<PToP>::create("pi->pi::KKRvD2024", p, Options{ });

                    TEST_CHECK(nullptr != ff);
                }

                /* f_+ and its auxiliary functions at spacelike and lightlike q2 <= 0.0 */
                {
                    KKRvD2024FormFactors<PiToPi> ff(p, Options{ });

                    const double chi = 0.0258815;

                    TEST_CHECK_NEARLY_EQUAL(ff.z(-1.0), -0.162626752, eps);
                    TEST_CHECK_NEARLY_EQUAL(ff.z( 0.0), -0.675539131, eps);

                    TEST_CHECK_NEARLY_EQUAL(ff.phi_p(ff.z(-1.0), chi), 0.145047910, eps);
                    TEST_CHECK_NEARLY_EQUAL(ff.phi_p(ff.z( 0.0), chi), 0.242821037, eps);

                    TEST_CHECK_NEARLY_EQUAL(ff.f_p(-1.0), 0.200915869,  eps);
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

                    const double chi = 0.0258815;

                    TEST_CHECK_NEARLY_EQUAL(real(ff.z( 0.0)),  -0.6755391, eps);
                    TEST_CHECK_NEARLY_EQUAL(imag(ff.z( 0.0)),  0.0,      eps);
                    TEST_CHECK_NEARLY_EQUAL(real(ff.z(+0.1)), -0.97897061, eps);
                    TEST_CHECK_NEARLY_EQUAL(imag(ff.z(+0.1)), -0.20400134, eps);
                    TEST_CHECK_NEARLY_EQUAL(real(ff.z(+0.5)), -0.66233531, eps);
                    TEST_CHECK_NEARLY_EQUAL(imag(ff.z(+0.5)), -0.74920754, eps);

                    TEST_CHECK_NEARLY_EQUAL(real(ff.phi_p(ff.z( 0.0), chi)),  0.242821037,  eps);
                    TEST_CHECK_NEARLY_EQUAL(imag(ff.phi_p(ff.z( 0.0), chi)),  0.0,        eps);
                    TEST_CHECK_NEARLY_EQUAL(real(ff.phi_p(ff.z(+0.1), chi)),  0.333648456, eps);
                    TEST_CHECK_NEARLY_EQUAL(imag(ff.phi_p(ff.z(+0.1), chi)),  0.107444433, eps);
                    TEST_CHECK_NEARLY_EQUAL(real(ff.phi_p(ff.z(+0.5), chi)),  0.142282074,  eps);
                    TEST_CHECK_NEARLY_EQUAL(imag(ff.phi_p(ff.z(+0.5), chi)),  0.140166151,  eps);

                    TEST_CHECK_NEARLY_EQUAL(real(ff.f_p( 0.0)),  1.0, eps);
                    TEST_CHECK_NEARLY_EQUAL(imag(ff.f_p( 0.0)),  0.0, eps);
                    TEST_CHECK_NEARLY_EQUAL(real(ff.f_p(+0.1)),  1.216681140,   eps);
                    TEST_CHECK_NEARLY_EQUAL(imag(ff.f_p(+0.1)),  -0.042359497,  eps);
                    TEST_CHECK_NEARLY_EQUAL(real(ff.f_p(+0.5)),  3.63653056,   eps);
                    TEST_CHECK_NEARLY_EQUAL(imag(ff.f_p(+0.5)),  3.73112339,   eps);
                }
            }
        }
} parametric_KKRvD2024_test;
