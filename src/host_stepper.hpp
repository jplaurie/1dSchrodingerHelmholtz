#pragma once

#include "integrator.hpp"
#include "parameters.hpp"
#include "types.hpp"

#include <array>
#include <stdexcept>

struct HostIntegrationWorkspace {
    std::array<SpectralField, 4> nonlinearStages;
    std::array<SpectralField, 3> stageStates;

    void initialize(std::size_t count, std::size_t nonlinearStageCount) {
        if (nonlinearStageCount < 2 || nonlinearStageCount > nonlinearStages.size())
            throw std::runtime_error("invalid host integration stage count");
        for (auto &field : nonlinearStages)
            field.clear();
        for (auto &field : stageStates)
            field.clear();
        for (std::size_t i = 0; i < nonlinearStageCount; ++i)
            nonlinearStages[i].resize(count);
        for (std::size_t i = 1; i < nonlinearStageCount; ++i)
            stageStates[i - 1].resize(count);
    }
};

template <class RightHandSide, class EnforceConstraints>
void advanceHostTimeStep(const Parameters &parameters, const IntegrationCoefficients &coefficients,
                         SpectralField &state, const SpectralField &noise,
                         HostIntegrationWorkspace &workspace, RightHandSide rightHandSide,
                         EnforceConstraints enforceConstraints) {
    if (!noise.empty() && noise.size() != state.size())
        throw std::runtime_error("invalid stochastic-noise field");

    auto &n1 = workspace.nonlinearStages[0];
    auto &n2 = workspace.nonlinearStages[1];
    auto &n3 = workspace.nonlinearStages[2];
    auto &n4 = workspace.nonlinearStages[3];
    auto &a = workspace.stageStates[0];
    auto &b = workspace.stageStates[1];
    auto &c = workspace.stageStates[2];

    rightHandSide(state, n1);
    for (std::size_t i = 0; i < state.size(); ++i)
        a[i] = integrationStageA(parameters.integrator, parameters.timeStep, i, coefficients,
                                 state[i], n1[i]);
    rightHandSide(a, n2);
    if (!n3.empty()) {
        for (std::size_t i = 0; i < state.size(); ++i)
            b[i] = integrationStageB(i, coefficients, state[i], n1[i], n2[i]);
        rightHandSide(b, n3);
    }
    if (!n4.empty()) {
        for (std::size_t i = 0; i < state.size(); ++i)
            c[i] = integrationStageC(i, coefficients, state[i], n1[i], n3[i]);
        rightHandSide(c, n4);
    }
    for (std::size_t i = 0; i < state.size(); ++i) {
        const IntegrationFinishValues values{.initial = state[i],
                                             .stageA = a[i],
                                             .nonlinear1 = n1[i],
                                             .nonlinear2 = n2[i],
                                             .nonlinear3 = n3.empty() ? Complex{} : n3[i],
                                             .nonlinear4 = n4.empty() ? Complex{} : n4[i]};
        state[i] =
            integrationFinish(parameters.integrator, parameters.timeStep, i, coefficients, values);
        if (!noise.empty())
            state[i] += noise[i];
    }
    enforceConstraints(state);
}
