#pragma once
#include <vector>
#include <utility>
#include <deque>
#include <tuple>
#include <random>
#include <algorithm>
#include <valarray>
#include <numeric>
#include <cassert>
#include <chrono>

namespace adrianyu {

class WalkerAlias {

public:

    WalkerAlias() : uni_01_dist(0, 1) {}

    /*
        prob_prop: need not to sum to one
    */
    void init(const std::vector<double> &prob_prop) {
        unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
        rnd_gen.seed(seed);
        probs.resize(prob_prop.size());
        assert(!prob_prop.empty());
        const double prob_p_sum = std::accumulate(prob_prop.begin(), prob_prop.end(), 0.0);
        probs.resize(prob_prop.size());
        for (size_t i = 0; i < probs.size(); ++i) {
            probs[i] = prob_prop[i] / prob_p_sum;
        }
        GenAlias(probs, Alias);
    }

    size_t operator()(void) {
        return SampleAlias(Alias);
    }

protected:
    std::vector<double> probs;
    std::mt19937_64 rnd_gen;
    std::uniform_real_distribution<double> uni_01_dist;
    std::vector<std::tuple<size_t, size_t, double> > Alias;

    static void GenAlias(const std::vector<double> & probs,
        std::vector<std::tuple<size_t, size_t, double> > & Alias)
    {
        const double uni_prob = 1.0 / static_cast<double>(probs.size());
        std::deque<std::pair<size_t, double> > L;
        std::deque<std::pair<size_t, double> > H;
        for (size_t i = 0; i < probs.size(); ++i) {
            if (probs[i] <= uni_prob) {
                L.push_back(std::make_pair(i, probs[i]));
            }
            else {
                H.push_back(std::make_pair(i, probs[i]));
            }
        }
        Alias.clear();
        while (!L.empty() && !H.empty()) {
            const std::pair<size_t, double> & L_item = L.front();
            const std::pair<size_t, double> & H_item = H.front();
            Alias.push_back(std::make_tuple(L_item.first, H_item.first, L_item.second));
            const double prob_rsd = H_item.second + L_item.second - uni_prob;
            if (prob_rsd > uni_prob) {
                H.push_back(std::make_pair(H_item.first, prob_rsd));
            }
            else {
                L.push_back(std::make_pair(H_item.first, prob_rsd));
            }
            H.pop_front();
            L.pop_front();
        }
        // if any of the H/L is not empty, we fill Alias with the same index.
        if (!L.empty()) {
            Alias.push_back(std::make_tuple(L.front().first, L.front().first, L.front().second));
        }
        else {
            if (!H.empty()) {
                Alias.push_back(std::make_tuple(H.front().first, H.front().first, H.front().second));
            }
        }
    }

    size_t SampleAlias(const std::vector<std::tuple<size_t, size_t, double> > & Alias) {
        uni_01_dist.reset();
        const double bin_p = uni_01_dist(rnd_gen);
        size_t bin = static_cast<size_t>(bin_p * static_cast<double>(Alias.size()));
        auto & tup = Alias[bin];
        uni_01_dist.reset();
        const double uni_p = uni_01_dist(rnd_gen);
        const double lp = static_cast<double>(Alias.size()) * std::get<2>(tup);
        // in real application, it is common that almost all prob in probs is near zero
        // hence, the value of lp is small.
        if (lp < uni_p) {
            return std::get<1>(tup);
        }
        else {
            return std::get<0>(tup);
        }
    }
};
}

