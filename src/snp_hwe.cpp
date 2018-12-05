/*
// This code implements an exact SNP test of Hardy-Weinberg Equilibrium as
// described in Wigginton, JE, Cutler, DJ, and Abecasis, GR (2005) A Note on
// Exact Tests of Hardy-Weinberg Equilibrium. AJHG 76: 887-893
//
// Written by Jan Wigginton
*/

#include <vector>
#include <algorithm>

#include "snp_hwe.h"

double SNPHWE(int obs_hets, int obs_hom1, int obs_hom2) {
    if (obs_hom1 < 0 || obs_hom2 < 0 || obs_hets < 0) {
        return -1.0;
    }

    int obs_homr = std::min(obs_hom1, obs_hom2);
    int obs_homc = std::max(obs_hom1, obs_hom2);

    int rare = 2 * obs_homr + obs_hets;
    int genotypes = obs_hets + obs_homc + obs_homr;

    if (genotypes == 0) {
        return -1.0;
    }

    std::vector<double> probs(rare + 1, 0.0);

    // get distribution midpoint, but ensure midpoint and rare alleles have
    // same parity
    int mid = rare * (2 * genotypes - rare) / (2 * genotypes);
    if (mid % 2 != rare % 2) {
        mid += 1;
    }

    probs[mid] = 1.0;
    double sum = probs[mid];
    
    int curr_homr = (rare - mid) / 2;
    int curr_homc = genotypes - mid - curr_homr;
    for (auto curr_hets = mid; curr_hets > 1; curr_hets -= 2) {
        probs[curr_hets - 2] = probs[curr_hets] * curr_hets * (curr_hets - 1.0)
                               / (4.0 * (curr_homr + 1.0) * (curr_homc + 1.0));
        sum += probs[curr_hets - 2];

        /* 2 fewer heterozygotes -> add one rare, one common homozygote */
        curr_homr += 1;
        curr_homc += 1;
    }

    // calculate probabilities from midpoint up
    curr_homr = (rare - mid) / 2;
    curr_homc = genotypes - mid - curr_homr;
    for (auto curr_hets = mid; curr_hets <= rare - 2; curr_hets += 2) {
        probs[curr_hets + 2] = probs[curr_hets] * 4.0 * curr_homr * curr_homc
                                / ((curr_hets + 2.0) * (curr_hets + 1.0));
        sum += probs[curr_hets + 2];

        /* add 2 heterozygotes -> subtract one rare, one common homozygote */
        curr_homr -= 1;
        curr_homc -= 1;
    }

    /*  p-value calculation for p_hwe  */
    double target = probs[obs_hets];
    double p_hwe = 0.0;
    for (auto p: probs) {
        if (p <= target) {
            p_hwe += p / sum;
        }
    }

    return std::min(1.0, p_hwe);
}
