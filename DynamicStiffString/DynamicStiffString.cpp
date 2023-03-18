/*
  ==============================================================================

    DynamicStiffString.cpp
    Created: 7 Feb 2022 5:09:58pm
    Author:  Silvin Willemsen

  ==============================================================================
*/

#include "DynamicStiffString.h"

#include <cassert>


//==============================================================================
DynamicStiffString::DynamicStiffString(const SimulationParameters& parameters, double k) : k(k)
{

    // Initialise member variables using the parameter set
    L = parameters.L;
    rho = parameters.rho;
    r = parameters.r;
    A = static_cast<double>(M_PI) * r * r;
    T = parameters.T;
    E = parameters.E;
    I = static_cast<double>(M_PI) * r * r * r * r * 0.25;
    sigma0 = parameters.sigma0;
    sigma1 = parameters.sigma1;

    origR = r;
    origT = T;
    origE = E;
    origL = L;
    origRho = rho;

    parameterPtrs.reserve(8);
    parameterPtrs.push_back(&L);
    parameterPtrs.push_back(&rho);
    parameterPtrs.push_back(&r);
    parameterPtrs.push_back(&T);
    parameterPtrs.push_back(&E);
    parameterPtrs.push_back(&sigma0);
    parameterPtrs.push_back(&sigma1);

    parametersToGoTo.resize(parameterPtrs.size(), 0);
    for (int i = 0; i < parameterPtrs.size(); ++i)
        parametersToGoTo[i] = *parameterPtrs[i];

    parameterChanged.resize(parameterPtrs.size(), false);

    double cSqMin = 0.5 * T / (2.0 * rho * 2.0 * r * 2.0 * r * static_cast<double>(M_PI));

    double hMin = sqrt(cSqMin * k * k + 4.0 * Global::sig1min * k);
    Nmax = static_cast<int>(floor(2.0 * L / hMin));

    double rMax = 0.5 * r; //?
    double cSqMax = 2.0 * T / (0.5 * rho * rMax * rMax * static_cast<double>(M_PI));
    double kappaSqMax = 2.0 * E * static_cast<double>(M_PI) * rMax * rMax * rMax * rMax * 0.25 / (0.5 * rho * rMax * rMax * static_cast<double>(M_PI));

    double hMax = sqrt((cSqMax * k * k + 4.0 * sigma1 * 2.0 * k + sqrt(pow(cSqMax * k * k + 4.0 * sigma1 * 2.0 * k, 2) + 16 * kappaSqMax * k * k)) / 2.0);

    double Nmin = floor(0.5 * L / hMax);

    // only add to left system (v)
    int MvMax = Nmax - Mw;

    // Initialise vectors (excluding outer boundaries
    vStates = std::vector<std::vector<double>>(3,
        std::vector<double>(MvMax + 1, 0));

    wStates = std::vector<std::vector<double>>(3,
        std::vector<double>(Mw + 1, 0));

    /*  Make u pointers point to the first index of the state vectors.
        To use u (and obtain a vector from the state vectors) use indices like u[n][l] where,
             - n = 0 is u^{n+1},
             - n = 1 is u^n, and
             - n = 2 is u^{n-1}.
        Also see calculateScheme()
     */

     // Initialise pointer vector
    v.resize(3, nullptr);
    w.resize(3, nullptr);

    // Make set memory addresses to first index of the state vectors.
    for (int i = 0; i < 3; ++i)
    {
        v[i] = &vStates[i][0];
        w[i] = &wStates[i][0];
    }
    customIp.resize(4, 0);

    refreshCoefficients(true);

    Nprev = N;
    NfracPrev = Nfrac;
}

DynamicStiffString::~DynamicStiffString()
{
}

void DynamicStiffString::calculateScheme()
{

    // simply supported left boundary
    v[0][1] = Bss * v[1][1] + B1 * (v[1][2] + v[1][0]) + B2 * v[1][3]
        + C0 * v[2][1] + C1 * (v[2][2] + v[2][0]);

    // main left scheme
    for (int l = 2; l < Mv - 1; ++l)
        v[0][l] = B0 * v[1][l] + B1 * (v[1][l + 1] + v[1][l - 1]) + B2 * (v[1][l + 2] + v[1][l - 2])
        + C0 * v[2][l] + C1 * (v[2][l + 1] + v[2][l - 1]);

    // inner boundary calculations
    A0 = Iterm * Iterm - 4.0 * Iterm + 6.0;
    A1 = Iterm - 4.0;
    A2 = -(Iterm * Iterm) + 4.0 * Iterm + 1.0;
    A3 = -Iterm;
    AA = Iterm - 2.0;

    if (Mw == 1)
    {
        // next-to-boundary point
        v[0][Mv - 1] = (2.0 * v[1][Mv - 1] - v[2][Mv - 1]
            + lambdaSq * (v[1][Mv] - 2.0 * v[1][Mv - 1] + v[1][Mv - 2])
            - muSq * (v[1][Mv - 3] - 4 * v[1][Mv - 2] + 6 * v[1][Mv - 1] + A1 * v[1][Mv] + w[1][0])
            + S0 * v[2][Mv - 1]
            + S1 * (v[1][Mv] - 2.0 * v[1][Mv - 1] + v[1][Mv - 2])
            - S1 * (v[2][Mv] - 2.0 * v[2][Mv - 1] + v[2][Mv - 2])) / (1.0 + S0);

        // boundary points
        v[0][Mv] = (2.0 * v[1][Mv]
            + lambdaSq * (w[1][0] + AA * v[1][Mv] + v[1][Mv - 1])
            - muSq * (v[1][Mv - 2] - 4 * v[1][Mv - 1] + A0 * v[1][Mv] + (A1 - A3) * w[1][0])
            + (-1.0 + S0) * v[2][Mv]
            + S1 * (w[1][0] + AA * v[1][Mv] + v[1][Mv - 1])
            - S1 * (w[2][0] + AA * v[2][Mv] + v[2][Mv - 1])) / (1.0 + S0);

        // right system (single point now)
        w[0][0] = (2.0 * w[1][0]
            + lambdaSq * (A3 * v[1][Mv - 1] + v[1][Mv] + AA * w[1][0]) // w[1][1] is 0
            - muSq * (A3 * v[1][Mv - 2] + A2 * v[1][Mv - 1] + A1 * v[1][Mv]
                + (A0 - 1.0) * w[1][0])  // w[1][1] is 0
            + (-1.0 + S0) * w[2][0]
            + S1 * (A3 * v[1][Mv - 1] + v[1][Mv] + AA * w[1][0])
            - S1 * (A3 * v[2][Mv - 1] + v[2][Mv] + AA * w[2][0])) / (1.0 + S0);
    }
    else
    {
        // next-to-boundary point (left)
        v[0][Mv - 1] = ((2.0 - 2.0 * lambdaSq - 6.0 * muSq - 2.0 * S1) * v[1][Mv - 1]
            + (lambdaSq + 4.0 * muSq + S1) * v[1][Mv - 2]
            + (lambdaSq - A1 * muSq + S1) * v[1][Mv]
            - muSq * (v[1][Mv - 3] + w[1][0] + A3 * w[1][1])
            + (S0 + 2.0 * S1 - 1.0) * v[2][Mv - 1]
            - S1 * (v[2][Mv] + v[2][Mv - 2])) / (1.0 + S0);

        // left inner boundary
        v[0][Mv] = ((2.0 + AA * (lambdaSq + S1) - A0 * muSq) * v[1][Mv]
            + (lambdaSq + 4.0 * muSq + S1) * v[1][Mv - 1]
            + (lambdaSq - A1 * muSq + S1) * w[1][0]
            + (A3 * lambdaSq - A2 * muSq + A3 * S1) * w[1][1]
            - muSq * (v[1][Mv - 2] + A3 * w[1][2])
            + (S0 - AA * S1 - 1.0) * v[2][Mv]
            - S1 * (v[2][Mv - 1] + w[2][0] + A3 * w[2][1])) / (1.0 + S0);

        // right inner boundary
        w[0][0] = ((2.0 + AA * (lambdaSq + S1) - A0 * muSq) * w[1][0]
            + (lambdaSq + 4.0 * muSq + S1) * w[1][1]
            + (lambdaSq - A1 * muSq + S1) * v[1][Mv]
            + (A3 * lambdaSq - A2 * muSq + A3 * S1) * v[1][Mv - 1]
            - muSq * (w[1][2] + A3 * v[1][Mv - 2])
            + (S0 - AA * S1 - 1.0) * w[2][0]
            - S1 * (w[2][1] + v[2][Mv] + A3 * v[2][Mv - 1])) / (1.0 + S0);

        if (Mw == 2)
        {
            // next-to-boundary point (right) + simply supported
            w[0][1] = (2.0 * w[1][1]
                + lambdaSq * (w[1][0] - 2.0 * w[1][1] + w[1][2])
                - muSq * (A3 * v[1][Mv - 1] + v[1][Mv] + A1 * w[1][0] + 5.0 * w[1][1] - 4.0 * w[1][2])
                + (-1.0 + S0) * w[2][1]
                + S1 * (w[1][0] - 2.0 * w[1][1] + w[1][2])
                - S1 * (w[2][0] - 2.0 * w[2][1] + w[2][2])) / (1.0 + S0);
        }
        else
        {
            // next-to-boundary point (right)
            w[0][1] = ((2.0 - 2.0 * lambdaSq - 6.0 * muSq - 2.0 * S1) * w[1][1]
                + (lambdaSq + 4.0 * muSq + S1) * w[1][2]
                + (lambdaSq - A1 * muSq + S1) * w[1][0]
                - muSq * (w[1][3] + v[1][Mv] + A3 * v[1][Mv - 1])
                + (S0 + 2.0 * S1 - 1.0) * w[2][1]
                - S1 * (w[2][0] + w[2][2])) / (1.0 + S0);

            // main right scheme
            for (int l = 2; l < Mw - 1; ++l)
                w[0][l] = B0 * w[1][l] + B1 * (w[1][l + 1] + w[1][l - 1]) + B2 * (w[1][l + 2] + w[1][l - 2])
                + C0 * w[2][l] + C1 * (w[2][l + 1] + w[2][l - 1]);

            // simply supported right boundary
            w[0][Mw - 1] = Bss * w[1][Mw - 1] + B1 * (w[1][Mw - 2] + w[1][Mw]) + B2 * w[1][Mw - 3]
                + C0 * w[2][Mw - 1] + C1 * (w[2][Mw - 2] + w[2][Mw]);

        }
    }
}

void DynamicStiffString::updateStates()
{
    // Do a pointer-switch. MUCH quicker than copying two entire state vectors every time-step.
    double* vTmp = v[2];
    v[2] = v[1];
    v[1] = v[0];
    v[0] = vTmp;

    double* wTmp = w[2];
    w[2] = w[1];
    w[1] = w[0];
    w[0] = wTmp;

    NfracPrev = Nfrac;
    Nprev = N;
}

void DynamicStiffString::excite(int loc)
{
    //// Arbitrary excitation function (raised cosine) ////

    // width (in grid points) of the excitation
    double width = 10;

    // make sure we're not going out of bounds at the left boundary
    int start = (loc == -1) ? static_cast<int>(std::max(floor((N + 1) * excitationLoc) - floor(width * 0.5), 1.0)) : loc;

    for (int l = 0; l < width; ++l)
    {
        // make sure we're not going out of bounds at the right boundary of the left system(this does 'cut off' the raised cosine)
        if (l + start > Mv)
            break;

        v[1][l + start] += 0.5 * (1 - cos(2.0 * static_cast<double>(M_PI) * l / (width - 1.0)));
        v[2][l + start] += 0.5 * (1 - cos(2.0 * static_cast<double>(M_PI) * l / (width - 1.0)));
    }
    // Disable the excitation flag to only excite once
    excitationFlag = false;
}


void DynamicStiffString::refreshParameter(int changedParameterIdx, double changedParameterValue)
{
    parametersToGoTo[changedParameterIdx] = changedParameterValue;
    parameterChanged[changedParameterIdx] = true;
}

void DynamicStiffString::refreshCoefficients(bool init)
{
    double NmaxChange = Global::NmaxChange;
    double paramDiffMax = 0.0;
    double NfracNext;

    bool needsRefresh = false;
    for (int i = 0; i < parameterPtrs.size(); ++i)
    {
        // if parameter hasn't changed, continue to next
        if (!parameterChanged[i])
            continue;

        needsRefresh = true;

        if (&L == parameterPtrs[i])
            paramDiffMax = NmaxChange * h;
        else if (&rho == parameterPtrs[i])
        {
            // bigger rho means bigger N
            NfracNext = Nfrac + (*parameterPtrs[i] > parametersToGoTo[i] ? -1 : 1) * NmaxChange;

            paramDiffMax = abs((k * k * L * L * NfracNext * NfracNext * T + 4.0 * I * k * k * NfracNext * NfracNext * NfracNext * NfracNext * E) / (A * L * L * (L * L - 4.0 * k * NfracNext * NfracNext * sigma1)) - rho);
        }
        else if (&T == parameterPtrs[i])
        {
            // bigger T means smaller N
            NfracNext = Nfrac + (*parameterPtrs[i] < parametersToGoTo[i] ? -1 : 1) * NmaxChange;

            paramDiffMax = abs((-4.0 * A * k * L * L * NfracNext * NfracNext * rho * sigma1 + A * L * L * L * L * rho - 4.0 * I * k * k * NfracNext * NfracNext * NfracNext * NfracNext * E) / (k * k * L * L * NfracNext * NfracNext) - T);
        }
        else if (&r == parameterPtrs[i])
        {

            if (E != 0)
            {
                double NfracNextPlus = Nfrac + NmaxChange;
                double NfracNextMin = Nfrac - NmaxChange;


                double bCoeffPlus = (16.0 * L * L * sigma1 * k) / (NfracNextPlus * NfracNextPlus) - (4.0 * L * L * L * L) / (NfracNextPlus * NfracNextPlus * NfracNextPlus * NfracNextPlus);
                double bCoeffMin = (16.0 * L * L * sigma1 * k) / (NfracNextMin * NfracNextMin) - (4.0 * L * L * L * L) / (NfracNextMin * NfracNextMin * NfracNextMin * NfracNextMin);

                std::vector<double> rVals(4, 0);

                // The graph of N (y-axis) vs r (x-axis) is a negative parabola. For a change in N (either positive or negative, there are 4 possible r values. Here we're trying to find the one that corresponds to the one we're trying to find.

                // r right side of parabola, increasing N
                rVals[0] = sqrt((-bCoeffPlus + sqrt(bCoeffPlus * bCoeffPlus - 16.0 * E * k * k / rho * (4.0 * L * L * T * k * k) / (NfracNextPlus * NfracNextPlus * rho * static_cast<double>(M_PI)))) / (8.0 * (E * k * k / rho)));
                // r right side of parabola, decreasing N
                rVals[1] = sqrt((-bCoeffMin + sqrt(bCoeffMin * bCoeffMin - 16.0 * E * k * k / rho * (4.0 * L * L * T * k * k) / (NfracNextMin * NfracNextMin * rho * static_cast<double>(M_PI)))) / (8.0 * (E * k * k / rho)));

                // r left side of parabola, increasing N
                rVals[2] = sqrt((-bCoeffPlus - sqrt(bCoeffPlus * bCoeffPlus - 16.0 * E * k * k / rho * (4.0 * L * L * T * k * k) / (NfracNextPlus * NfracNextPlus * rho * static_cast<double>(M_PI)))) / (8.0 * (E * k * k / rho)));

                // r left side of parabola, decreasing N
                rVals[3] = sqrt((-bCoeffMin - sqrt(bCoeffMin * bCoeffMin - 16.0 * E * k * k / rho * (4.0 * L * L * T * k * k) / (NfracNextMin * NfracNextMin * rho * static_cast<double>(M_PI)))) / (8.0 * (E * k * k / rho)));

                double rDiff = 1;
                double rToGoTo = parametersToGoTo[i];
                int idxToChoose = -1;
                for (int i = 0; i < rVals.size(); ++i)
                {
                    if (isnan(rVals[i]))
                        continue;
                    // if r is decreased, don't choose larger r values
                    if (rToGoTo < r && rVals[i] > r)
                        continue;

                    // if r is increased, don't choose smaller r values
                    if (rToGoTo > r && rVals[i] < r)
                        continue;

                    if (abs(rVals[i] - r) < rDiff)
                    {
                        rDiff = rVals[i] - r;
                        idxToChoose = i;
                    }
                }
                paramDiffMax = abs(rVals[idxToChoose] - r);

            }
            else
            {
                // if E = 0, bigger r means bigger N
                NfracNext = Nfrac + (*parameterPtrs[i] > parametersToGoTo[i] ? -1 : 1) * NmaxChange;

                paramDiffMax = abs((k * NfracNext * sqrt(T)) / (sqrt(rho) * sqrt(static_cast<double>(M_PI) * L * L - 4.0 * static_cast<double>(M_PI) * k * NfracNext * NfracNext * sigma1)) - r);
            }

        }
        else if (&E == parameterPtrs[i])
        {
            // bigger E means smaller N
            NfracNext = Nfrac + (*parameterPtrs[i] < parametersToGoTo[i] ? -1 : 1) * NmaxChange;

            paramDiffMax = abs((-4.0 * A * k * L * L * NfracNext * NfracNext * rho * sigma1 + A * L * L * L * L * rho - k * k * L * L * NfracNext * NfracNext * T) / (4.0 * I * k * k * NfracNext * NfracNext * NfracNext * NfracNext) - E);

        }
        else if (&sigma1 == parameterPtrs[i])
        {
            // bigger sigma1 means smaller N
            NfracNext = Nfrac + (*parameterPtrs[i] < parametersToGoTo[i] ? -1 : 1) * NmaxChange;

            paramDiffMax = abs((A * L * L * L * L * rho - k * k * L * L * NfracNext * NfracNext * T - 4.0 * I * k * k * NfracNext * NfracNext * NfracNext * NfracNext * E) / (4.0 * A * k * L * L * NfracNext * NfracNext * rho) - sigma1);

        }
        else if (&sigma0 == parameterPtrs[i])
        {
            //            *parameterPtrs[i] = parametersToGoTo[i];
            paramDiffMax = 100;
        }
        //    L = (1-LfilterCoeff) * LtoGoTo + LfilterCoeff * Lpre
        if (abs(*parameterPtrs[i] - parametersToGoTo[i]) < paramDiffMax)
        {
            *parameterPtrs[i] = parametersToGoTo[i];
            parameterChanged[i] = false;
        }
        else if (*parameterPtrs[i] < parametersToGoTo[i])
            *parameterPtrs[i] += paramDiffMax;
        else if (*parameterPtrs[i] > parametersToGoTo[i])
            *parameterPtrs[i] -= paramDiffMax;
        //        if (parameterPtrs[i] == &r || parameterPtrs[i] == &E)
        //        {
        //            *parameterPtrs[i] = 0.9995 * (*parameterPtrs[i]) + 0.0005 * parametersToGoTo[i];
        //
        //        } else {
        //            *parameterPtrs[i] = 0.999 * (*parameterPtrs[i]) + 0.001 * parametersToGoTo[i];
        //        }
    }

    // if the parameters don't need refresh, return
    if (!needsRefresh && !init)
        return;

    A = static_cast<double>(M_PI) * r * r;
    I = static_cast<double>(M_PI) * r * r * r * r * 0.25;

    // Calculate wave speed (squared)
    cSq = T / (rho * A);

    // Calculate stiffness coefficient (squared)
    kappaSq = E * I / (rho * A);

    double stabilityTerm = cSq * k * k + 4.0 * sigma1 * k; // just easier to write down below

    h = sqrt(0.5 * (stabilityTerm + sqrt((stabilityTerm * stabilityTerm) + 16.0 * kappaSq * k * k)));
    Nfrac = L / h;

    // check if the change does not surpass a limit
    N = static_cast<int>(floor(Nfrac));
    alf = Nfrac - N;
    if (init)
        Nprev = N;

    // Check whether a grid point needs to be added or removed
    if (Nprev != N)
        addRemovePoint();

    Mv = N - Mw;

    Iterm = (alf - 1.0) / (alf + 1.0);

    lambdaSq = cSq * k * k / (h * h);
    muSq = kappaSq * k * k / (h * h * h * h);

    // Coefficients used for damping
    S0 = sigma0 * k;
    S1 = (2.0 * sigma1 * k) / (h * h);

    // Scheme coefficients
    B0 = 2.0 - 2.0 * lambdaSq - 6.0 * muSq - 2.0 * S1; // u_l^n
    Bss = 2.0 - 2.0 * lambdaSq - 5.0 * muSq - 2.0 * S1;
    B1 = lambdaSq + 4.0 * muSq + S1;                   // u_{l+-1}^n
    B2 = -muSq;                                        // u_{l+-2}^n
    C0 = -1.0 + S0 + 2.0 * S1;                         // u_l^{n-1}
    C1 = -S1;                                          // u_{l+-1}^{n-1}

    Adiv = 1.0 / (1.0 + S0);                           // u_l^{n+1}

    // Divide by u_l^{n+1} term
    B0 *= Adiv;
    Bss *= Adiv;
    B1 *= Adiv;
    B2 *= Adiv;
    C0 *= Adiv;
    C1 *= Adiv;
}

void DynamicStiffString::addRemovePoint()
{
    assert(abs(N - Nprev) <= 1);
    refreshCustomIp();
    if (N > Nprev)
    {
        // possibly unnecessary to update up[0]
        v[0][Mv + 1] = customIp[0] * v[0][Mv - 1]
            + customIp[1] * v[0][Mv]
            + customIp[2] * w[0][0]
            + customIp[3] * w[0][1];

        v[1][Mv + 1] = customIp[0] * v[1][Mv - 1]
            + customIp[1] * v[1][Mv]
            + customIp[2] * w[1][0]
            + customIp[3] * w[1][1];

        v[2][Mv + 1] = customIp[0] * v[2][Mv - 1]
            + customIp[1] * v[2][Mv]
            + customIp[2] * w[2][0]
            + customIp[3] * w[2][1];

    }
    else {
        v[0][Mv] = 0;
        v[1][Mv] = 0;
        v[2][Mv] = 0;
    }
}

void DynamicStiffString::refreshCustomIp()
{
    customIp[0] = -alf * (alf + 1.0) / ((alf + 2.0) * (alf + 3.0));
    customIp[1] = 2.0 * alf / (alf + 2.0);
    customIp[2] = 2.0 / (alf + 2.0);
    customIp[3] = -2.0 * alf / ((alf + 3.0) * (alf + 2.0));
}
