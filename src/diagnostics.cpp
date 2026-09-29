#include "output.hpp"

#include "spectral.hpp"

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <stdexcept>

namespace {
bool finiteValues(const std::vector<double> &values) {
    return std::all_of(values.begin(), values.end(),
                       [](double value) { return std::isfinite(value); });
}

std::ofstream numericOutput(const std::filesystem::path &path,
                            std::ios::openmode mode = std::ios::out) {
    std::ofstream output(path, mode);
    if (!output)
        throw std::runtime_error("cannot open output file: " + path.string());
    output << std::scientific;
    return output;
}

void closeChecked(std::ofstream &output, const std::string &message) {
    output.close();
    if (!output)
        throw std::runtime_error(message);
}
} // namespace

double appendDiagnostics(const Parameters &parameters, NonlinearOperator &nonlinearOperator,
                         double time, std::uint64_t frame, const SpectralField &wavefunction,
                         DiagnosticsAverages &averages) {
    const std::size_t bins = parameters.gridPoints / 2;
    std::vector<double> spectrum(bins), quadraticSpectrum(bins);
    std::vector<double> waveTransfer(bins), energyTransfer(bins);
    SpectralField density, nonlinearRhs, wavefunctionDot, densityDot;
    nonlinearOperator.densitySpectrum(wavefunction, density);
    nonlinearOperator.evaluate(wavefunction, nonlinearRhs);
    wavefunctionDot.resize(parameters.gridPoints);

    double waveAction = 0.0, linearEnergy = 0.0, chemicalEnergy = 0.0;
    double waveHypo = 0.0, waveHyper = 0.0;
    double energyHypo = 0.0, energyHyper = 0.0;
    for (std::size_t i = 0; i < parameters.gridPoints; ++i) {
        const long mode = signedWave(i, parameters.gridPoints);
        if (std::abs(mode) >= static_cast<long>(bins))
            continue;
        const std::size_t bin = static_cast<std::size_t>(std::abs(mode));
        const double k = waveNumber(parameters, i), k2 = k * k;
        const double power = std::norm(wavefunction[i]);
        const double hamiltonian =
            -parameters.dispersionCoefficient * k2 + parameters.chemicalPotential;
        waveAction += parameters.domainLength * power;
        linearEnergy += parameters.domainLength * (-parameters.dispersionCoefficient * k2) * power;
        chemicalEnergy += parameters.domainLength * parameters.chemicalPotential * power;
        spectrum[bin] += power;
        quadraticSpectrum[bin] += hamiltonian * power;
        wavefunctionDot[i] = Complex(0.0, -hamiltonian) * wavefunction[i] + nonlinearRhs[i];
        if (parameters.hypoviscosity > 0.0 && (k2 > 0.0 || parameters.hypoviscosityOrder >= 0.0)) {
            const double rate =
                parameters.hypoviscosity * std::pow(k2, parameters.hypoviscosityOrder);
            waveHypo += 2.0 * parameters.domainLength * rate * power;
            energyHypo += 2.0 * parameters.domainLength * rate * hamiltonian * power;
        }
        if (parameters.hyperviscosity > 0.0) {
            const double rate =
                parameters.hyperviscosity * std::pow(k2, parameters.hyperviscosityOrder);
            waveHyper += 2.0 * parameters.domainLength * rate * power;
            energyHyper += 2.0 * parameters.domainLength * rate * hamiltonian * power;
        }
    }
    nonlinearOperator.densityTimeDerivativeSpectrum(wavefunction, wavefunctionDot, densityDot);

    double nonlinearEnergy = 0.0;
    for (std::size_t i = 0; i < density.size(); ++i) {
        const double k = waveNumber(parameters, i);
        nonlinearEnergy += 0.5 * parameters.nonlinearityCoefficient * parameters.domainLength *
                           densityResponseMultiplier(parameters, k * k) * std::norm(density[i]);
    }

    for (std::size_t i = 0; i < parameters.gridPoints; ++i) {
        const long mode = signedWave(i, parameters.gridPoints);
        if (std::abs(mode) >= static_cast<long>(bins))
            continue;
        const std::size_t bin = static_cast<std::size_t>(std::abs(mode));
        const double k = waveNumber(parameters, i);
        const double wave = -2.0 * parameters.domainLength *
                            std::real(std::conj(wavefunction[i]) * wavefunctionDot[i]);
        const Complex gradient = Complex(0.0, k) * wavefunction[i];
        const Complex gradientDot = Complex(0.0, k) * wavefunctionDot[i];
        const double energy = 2.0 * parameters.dispersionCoefficient * parameters.domainLength *
                                  std::real(gradient * std::conj(gradientDot)) +
                              parameters.chemicalPotential * wave -
                              parameters.nonlinearityCoefficient * parameters.domainLength *
                                  densityResponseMultiplier(parameters, k * k) *
                                  std::real(density[i] * std::conj(densityDot[i]));
        waveTransfer[bin] += wave;
        energyTransfer[bin] += energy;
    }
    for (std::size_t i = 1; i < bins; ++i) {
        waveTransfer[i] += waveTransfer[i - 1];
        energyTransfer[i] += energyTransfer[i - 1];
    }

    if (averages.spectrum.empty()) {
        averages.spectrum.assign(bins, 0.0);
        averages.waveFlux.assign(bins, 0.0);
        averages.energyFlux.assign(bins, 0.0);
    }
    ++averages.count;
    for (std::size_t i = 0; i < bins; ++i) {
        averages.spectrum[i] += spectrum[i];
        averages.waveFlux[i] += waveTransfer[i];
        averages.energyFlux[i] += energyTransfer[i];
    }

    const double totalEnergy = linearEnergy + chemicalEnergy + nonlinearEnergy;
    if (!std::isfinite(totalEnergy) || !std::isfinite(waveAction) || !std::isfinite(linearEnergy) ||
        !std::isfinite(chemicalEnergy) || !std::isfinite(nonlinearEnergy) ||
        !std::isfinite(waveHypo) || !std::isfinite(waveHyper) || !std::isfinite(energyHypo) ||
        !std::isfinite(energyHyper) || !finiteValues(spectrum) ||
        !finiteValues(quadraticSpectrum) || !finiteValues(waveTransfer) ||
        !finiteValues(energyTransfer) || !finiteValues(averages.spectrum) ||
        !finiteValues(averages.waveFlux) || !finiteValues(averages.energyFlux))
        throw std::runtime_error("diagnostics contain a non-finite value");
    {
        auto out = numericOutput(parameters.outputDirectory / "diagnostics.csv", std::ios::app);
        out << std::setprecision(17) << time << ',' << frame << ',' << waveAction << ','
            << linearEnergy << ',' << chemicalEnergy << ',' << nonlinearEnergy << ',' << totalEnergy
            << ',' << waveHypo << ',' << waveHyper << ',' << energyHypo << ',' << energyHyper
            << '\n';
        closeChecked(out, "failed while writing diagnostics.csv");
    }
    {
        auto out = numericOutput(parameters.outputDirectory / "spectra.csv", std::ios::app);
        const double dk = 2.0 * sh1dPi / parameters.domainLength;
        for (std::size_t i = 0; i < bins; ++i)
            out << time << ',' << frame << ',' << dk * static_cast<double>(i) << ',' << spectrum[i]
                << ',' << quadraticSpectrum[i] << ','
                << averages.spectrum[i] / static_cast<double>(averages.count) << '\n';
        closeChecked(out, "failed while writing spectra.csv");
    }
    {
        auto out = numericOutput(parameters.outputDirectory / "fluxes.csv", std::ios::app);
        const double dk = 2.0 * sh1dPi / parameters.domainLength;
        for (std::size_t i = 0; i < bins; ++i)
            out << time << ',' << frame << ',' << dk * static_cast<double>(i) << ','
                << waveTransfer[i] << ',' << energyTransfer[i] << ','
                << averages.waveFlux[i] / static_cast<double>(averages.count) << ','
                << averages.energyFlux[i] / static_cast<double>(averages.count) << '\n';
        closeChecked(out, "failed while writing fluxes.csv");
    }
    if (parameters.writeModeDiagnostics) {
        auto out = numericOutput(parameters.outputDirectory / "modes.csv", std::ios::app);
        for (std::size_t i = 0; i < parameters.gridPoints; ++i)
            out << time << ',' << frame << ',' << signedWave(i, parameters.gridPoints) << ','
                << waveNumber(parameters, i) << ',' << wavefunction[i].real() << ','
                << wavefunction[i].imag() << '\n';
        closeChecked(out, "failed while writing modes.csv");
    }
    return totalEnergy;
}
