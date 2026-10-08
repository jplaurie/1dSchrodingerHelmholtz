#pragma once

#include <filesystem>
#include <fstream>
#include <iomanip>
#include <stdexcept>
#include <string>
#include <system_error>

inline std::ofstream numericOutput(const std::filesystem::path &path,
                                   std::ios::openmode mode = std::ios::out) {
    std::ofstream output(path, mode);
    if (!output)
        throw std::runtime_error("cannot open output file: " + path.string());
    output << std::scientific << std::setprecision(17);
    return output;
}

inline void closeChecked(std::ofstream &stream, const std::string &error) {
    stream.close();
    if (!stream)
        throw std::runtime_error(error);
}

template <class Writer> void writeAtomic(const std::filesystem::path &path, Writer writer) {
    const std::filesystem::path temporary = path.string() + ".tmp";
    std::error_code ignored;
    std::filesystem::remove(temporary, ignored);
    try {
        writer(temporary);
        std::filesystem::rename(temporary, path);
    } catch (...) {
        std::filesystem::remove(temporary, ignored);
        throw;
    }
}
