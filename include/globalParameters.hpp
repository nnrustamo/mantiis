/*
 * Copyright (c) 2024 Nijat Rustamov
 *
 * Author: Nijat Rustamov
 * Organization: University of Wyoming
 * Email: nrustamo@uwyo.edu
 *
 * Academic Supervisor: Saman Aryana
 * Email: saryana@uwyo.edu
 * Organization: University of Wyoming
 *
 * This file is a part of Lattice Boltzmann Simulation Software
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in all
 * copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 */

/*
    Globally Accessible Simulation Parameters
*/

#pragma once

#include <math.h>
#include <vector>
#include <limits>

#include "utils.hpp"

namespace _GLOBAL_
{   
    std::string system_directory = "";
    // ================================================= Simulation Conditions
    // these parameters can be modified from elsewhere
    double Cl = 1.0e-9;     // resolution
    double Fbody = 1e-8;  // body force for transport
    double T_phy = 300.0; // physical temperature, K;
    double P_phy = 2.0e6; // physical pressure, Pa;

    // ================================================= Fluid specific parameters
    // methane case, these parameters can be modified from elsewhere
    double Tc_phy = 190.564;        // methane bulk critical temperature, K
    double Pc_phy = 4.5992 * 1.0e6; // methane bulk critical pressure, Pa;
    double d_mol = 3.800e-10;       // methane molecular diameter, meters
    double omega = 0.011;           // methane accentric factor
    double Mw = 16.043 * 1.0e-3;    // methane molecular weight, kg/mol
    //

    // ================================================= Solving PR-EOS
    // the followig paremeters need not to be accessible outside this file
    // ================================================= Constants
    // these parameters cannot be modified from elsewhere
    const double NA = 6.02214076e23; // Avagadro constant
    const double R_EOS = 8.31446;    // Universal Gas Constant, m3*Pa/K/mol
    // lattice unit constants
    const double a_lu = 2.0 / 49.0; // PR-EOS constant
    const double b_lu = 2.0 / 21.0; // PR-EOS constant
    // Derived values that depend on runtime settings (omega, T_phy, P_phy, etc.)
    double k = 0.37464 + 1.54226 * omega - 0.26992 * pow(omega, 2);
    double alpha = pow((1 + k * (1 - sqrt(T_phy / Tc_phy))), 2);
    double a_alpha_lu = a_lu * alpha;
    double Rs = R_EOS / Mw;
    const double R_lu = 1.0;
    const double Mw_lu = 1.0;
    const double Rs_lu = R_lu / Mw_lu;
    //
    double a = 0.45724 * (R_EOS * Tc_phy) * (R_EOS * Tc_phy) / Pc_phy;
    double aG = a * alpha;
    double bG = 0.07780 * R_EOS * Tc_phy / Pc_phy;
    double AG = aG * P_phy / (R_EOS * R_EOS) / (T_phy * T_phy);
    double BG = bG * P_phy / R_EOS / T_phy;

    // Z3*Z^3 + Z2*Z^2 + Z1*Z + Z0 = 0
    const double ZG3 = 1;
    double ZG2 = BG - 1;
    double ZG1 = AG - 3 * (BG * BG) - 2 * BG;
    double ZG0 = BG * BG * BG + (BG * BG) - AG * BG;
    std::vector<double> ZG_solution = utils::solveCubic(ZG3, ZG2, ZG1, ZG0);
    double ZG = utils::findMaxValue(ZG_solution);
    double vG = ZG * R_EOS * T_phy / P_phy;
    double rho_phy_mol = 1 / vG;       //  mol/m^3 - molar density
    double rho_phy = rho_phy_mol * Mw; // kg/m^3 - density

    // conversion constants C_
    double C_Rs = Rs / Rs_lu;                               // (m2/s2/K)/(m_lu2/s_lu2/K_lu)      = Cl^2/Ct2/CT
    double b_Rs = 0.07780 * Rs * Tc_phy / Pc_phy;           // m3/kg
    double a_Rs = 0.45724 * pow((Rs * Tc_phy), 2) / Pc_phy; // m5/kg/s2
    double C_as = a_Rs / a_lu;                              // (m5/kg/s2)/(m_lu5/kg_lu/s_lu2)    = Cl^5/Ct2/Cm
    double C_bs = b_Rs / b_lu;                              //(m3/kg)/(m_lu3/kg_lu)
    double C_aob = C_as / C_bs;                             // (m2/s2)/(m_lu2/s_lu2)             = Cl^2/Ct2
    double Ct = Cl / sqrt(C_aob);                           // conversion factor of time, s/s_lu
    double CT = 1 / (C_Rs / pow(Cl, 2) * pow(Ct, 2));       // conversion factor of temperature, K/K_lu
    double T_lu = T_phy / CT;
    double C_Mw = Mw / Mw_lu;          //(kg/mol)/(kg_lu/mol_lu)           = Cm/Cn
    double Cm = pow(Cl, 3) / C_bs;     // conversion factor of mass, kg/kg_lu
    double Crho = Cm / (Cl * Cl * Cl); // conversion factor of density, kg/m3
    double Cn = Cm / C_Mw;             // conversion factor of amount, mol/mol_lu
    double Cu = Cl / Ct;               // conversion factor of velocity, m/s
    double Cf = Cm / (pow(Ct, 2) * pow(Cl, 2));               // conversion factor of force per unit volume, kg/s2/m2
    
    double Crho_mol = Cn / pow(Cl, 3); // conversion factor of mole density, mol/m3
    double rhoin = rho_phy / Crho;

    // ================================================= Interaction Force constants
    const double Gfs = -12.7; // solid-fluid interaction strength
    const double Gff = -1.0;  // fluid-fluid interaction strength
    const std::vector<double> wf{0.0, 1.0 / 3.0, 1.0 / 3.0, 1.0 / 3.0, 1.0 / 3.0, 1.0 / 12.0, 1.0 / 12.0, 1.0 / 12.0, 1.0 / 12.0};

    // mean free path
    double mfp = 1.0 / (sqrt(2.0) * Crho_mol * rhoin / Mw_lu * NA * M_PI * (2 * d_mol / 2.0 + 2 * d_mol / 2.0) * (2 * d_mol / 2.0 + 2 * d_mol / 2.0));

    // ================================================= Boundary parameters
    const double A1 = 1.2; // half-way bounce-back C1
    double r = 0.5; //2.0 * A1 / (sqrt(6.0 / M_PI) + A1);
    const double A2 = 0.9;        // 1.0 / M_PI + 0.5 * pow(A1, 2); //

    // Specular reflection boundary files
    int N_bnd_types = 52;
    std::vector<std::vector<int>> icsr(N_bnd_types, std::vector<int>(9, 0));
    std::vector<std::vector<bool>> ifStream(N_bnd_types, std::vector<bool>(9, false));
    std::vector<std::vector<double>> xy_norm(N_bnd_types, std::vector<double>(2, 0.0));
    std::vector<int> bnd_types(N_bnd_types, 0);
    
    // ================================================= Collision parameters
    double tau_SRT = 1.0;
    // Identity Matrix
    const std::vector<double> I{
        1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
        0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
        0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
        0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0,
        0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0,
        0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0,
        0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0,
        0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0,
        0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0};

    // MRT collision matrix
    const std::vector<double> M{
        1., 1., 1., 1., 1., 1., 1., 1., 1.,
        -4., -1., -1., -1., -1., 2., 2., 2., 2.,
        4., -2., -2., -2., -2., 1., 1., 1., 1.,
        0., 1., 0., -1., 0., 1., -1., -1., 1.,
        0., -2., 0., 2., 0., 1., -1., -1., 1.,
        0., 0., 1., 0., -1., 1., 1., -1., -1.,
        0., 0., -2., 0., 2., 1., 1., -1., -1.,
        0., 1., -1., 1., -1., 0., 0., 0., 0.,
        0., 0., 0., 0., 0., 1., -1., 1., -1.};

    // Inverse of MRT collision matrix
    const std::vector<double> Minv{
        0.111111111111111, -0.111111111111111, 0.111111111111111, 0.000000000000000, -0.000000000000000, -0.000000000000000, -0.000000000000000, 0.000000000000000, 0.000000000000000,
        0.111111111111111, -0.027777777777778, -0.055555555555556, 0.166666666666667, -0.166666666666667, 0.000000000000000, 0.000000000000000, 0.250000000000000, -0.000000000000000,
        0.111111111111111, -0.027777777777778, -0.055555555555556, -0.000000000000000, -0.000000000000000, 0.166666666666667, -0.166666666666667, -0.250000000000000, 0.000000000000000,
        0.111111111111111, -0.027777777777778, -0.055555555555556, -0.166666666666667, 0.166666666666667, -0.000000000000000, -0.000000000000000, 0.250000000000000, -0.000000000000000,
        0.111111111111111, -0.027777777777778, -0.055555555555556, 0.000000000000000, 0.000000000000000, -0.166666666666667, 0.166666666666667, -0.250000000000000, -0.000000000000000,
        0.111111111111111, 0.055555555555556, 0.027777777777778, 0.166666666666667, 0.083333333333333, 0.166666666666667, 0.083333333333333, 0.000000000000000, 0.250000000000000,
        0.111111111111111, 0.055555555555556, 0.027777777777778, -0.166666666666667, -0.083333333333333, 0.166666666666667, 0.083333333333333, 0.000000000000000, -0.250000000000000,
        0.111111111111111, 0.055555555555556, 0.027777777777778, -0.166666666666667, -0.083333333333333, -0.166666666666667, -0.083333333333333, 0.000000000000000, 0.250000000000000,
        0.111111111111111, 0.055555555555556, 0.027777777777778, 0.166666666666667, 0.083333333333333, -0.166666666666667, -0.083333333333333, 0.000000000000000, -0.250000000000000};
    // MRT relaxation matrix
    std::vector<double> R;

    // Recompute all derived globals after changing inputs like T_phy, P_phy, Cl, omega, Mw, etc.
    inline void recompute()
    {
        k = 0.37464 + 1.54226 * omega - 0.26992 * pow(omega, 2);
        alpha = pow((1 + k * (1 - sqrt(T_phy / Tc_phy))), 2);
        a_alpha_lu = a_lu * alpha;
        Rs = R_EOS / Mw;

        a = 0.45724 * (R_EOS * Tc_phy) * (R_EOS * Tc_phy) / Pc_phy;
        aG = a * alpha;
        bG = 0.07780 * R_EOS * Tc_phy / Pc_phy;
        AG = aG * P_phy / (R_EOS * R_EOS) / (T_phy * T_phy);
        BG = bG * P_phy / R_EOS / T_phy;

        ZG2 = BG - 1;
        ZG1 = AG - 3 * (BG * BG) - 2 * BG;
        ZG0 = BG * BG * BG + (BG * BG) - AG * BG;
        ZG_solution = utils::solveCubic(ZG3, ZG2, ZG1, ZG0);
        ZG = utils::findMaxValue(ZG_solution);
        vG = ZG * R_EOS * T_phy / P_phy;
        rho_phy_mol = 1 / vG;
        rho_phy = rho_phy_mol * Mw;

        C_Rs = Rs / Rs_lu;
        b_Rs = 0.07780 * Rs * Tc_phy / Pc_phy;
        a_Rs = 0.45724 * pow((Rs * Tc_phy), 2) / Pc_phy;
        C_as = a_Rs / a_lu;
        C_bs = b_Rs / b_lu;
        C_aob = C_as / C_bs;
        Ct = Cl / sqrt(C_aob);
        CT = 1 / (C_Rs / pow(Cl, 2) * pow(Ct, 2));
        T_lu = T_phy / CT;
        C_Mw = Mw / Mw_lu;
        Cm = pow(Cl, 3) / C_bs;
        Crho = Cm / (Cl * Cl * Cl);
        Cn = Cm / C_Mw;
        Cu = Cl / Ct;
        Cf = Cm / (pow(Ct, 2) * pow(Cl, 2));
        Crho_mol = Cn / pow(Cl, 3);
        rhoin = rho_phy / Crho;

        mfp = 1.0 / (sqrt(2.0) * Crho_mol * rhoin / Mw_lu * NA * M_PI * (2 * d_mol / 2.0 + 2 * d_mol / 2.0) * (2 * d_mol / 2.0 + 2 * d_mol / 2.0));
    }
}
