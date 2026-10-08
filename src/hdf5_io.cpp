#include "hdf5_io.hpp"

#include <algorithm>
#include <cmath>
#include <cstring>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

#ifdef SH1D_HAVE_HDF5
#include <hdf5.h>

namespace {
template <herr_t (*Close)(hid_t)> class Hdf5Handle {
  public:
    Hdf5Handle() = default;
    explicit Hdf5Handle(hid_t id) : id_(id) {}
    ~Hdf5Handle() {
        if (id_ >= 0)
            Close(id_);
    }
    Hdf5Handle(const Hdf5Handle &) = delete;
    Hdf5Handle &operator=(const Hdf5Handle &) = delete;
    Hdf5Handle(Hdf5Handle &&other) noexcept : id_(other.id_) { other.id_ = -1; }
    Hdf5Handle &operator=(Hdf5Handle &&other) noexcept {
        if (this != &other) {
            if (id_ >= 0)
                Close(id_);
            id_ = other.id_;
            other.id_ = -1;
        }
        return *this;
    }
    [[nodiscard]] hid_t get() const { return id_; }

  private:
    hid_t id_ = -1;
};

using File = Hdf5Handle<H5Fclose>;
using DataSet = Hdf5Handle<H5Dclose>;
using DataSpace = Hdf5Handle<H5Sclose>;
using PropertyList = Hdf5Handle<H5Pclose>;
using Attribute = Hdf5Handle<H5Aclose>;
using DataType = Hdf5Handle<H5Tclose>;

hid_t requireId(hid_t id, const std::string &operation) {
    if (id < 0)
        throw std::runtime_error(operation);
    return id;
}

void requireSuccess(herr_t status, const std::string &operation) {
    if (status < 0)
        throw std::runtime_error(operation);
}

template <class T>
void writeScalarAttribute(hid_t object, const char *name, hid_t type, const T &value) {
    DataSpace space(requireId(H5Screate(H5S_SCALAR), "cannot create HDF5 attribute space"));
    Attribute attribute(
        requireId(H5Acreate2(object, name, type, space.get(), H5P_DEFAULT, H5P_DEFAULT),
                  std::string("cannot create HDF5 attribute: ") + name));
    requireSuccess(H5Awrite(attribute.get(), type, &value),
                   std::string("cannot write HDF5 attribute: ") + name);
}

void writeStringAttribute(hid_t object, const char *name, const char *value) {
    DataType type(requireId(H5Tcopy(H5T_C_S1), "cannot create HDF5 string type"));
    requireSuccess(H5Tset_size(type.get(), std::strlen(value) + 1), "cannot size HDF5 string type");
    requireSuccess(H5Tset_strpad(type.get(), H5T_STR_NULLTERM),
                   "cannot configure HDF5 string type");
    DataSpace space(requireId(H5Screate(H5S_SCALAR), "cannot create HDF5 attribute space"));
    Attribute attribute(
        requireId(H5Acreate2(object, name, type.get(), space.get(), H5P_DEFAULT, H5P_DEFAULT),
                  std::string("cannot create HDF5 attribute: ") + name));
    requireSuccess(H5Awrite(attribute.get(), type.get(), value),
                   std::string("cannot write HDF5 attribute: ") + name);
}

template <class T> T readScalarAttribute(hid_t object, const char *name, hid_t type) {
    Attribute attribute(requireId(H5Aopen(object, name, H5P_DEFAULT),
                                  std::string("missing HDF5 attribute: ") + name));
    T value{};
    requireSuccess(H5Aread(attribute.get(), type, &value),
                   std::string("cannot read HDF5 attribute: ") + name);
    return value;
}

std::string readStringAttribute(hid_t object, const char *name) {
    Attribute attribute(requireId(H5Aopen(object, name, H5P_DEFAULT),
                                  std::string("missing HDF5 attribute: ") + name));
    DataType type(requireId(H5Aget_type(attribute.get()), "cannot inspect HDF5 string attribute"));
    const std::size_t size = H5Tget_size(type.get());
    if (size == 0 || size > 1024)
        throw std::runtime_error("invalid HDF5 string attribute size");
    std::vector<char> value(size + 1, '\0');
    requireSuccess(H5Aread(attribute.get(), type.get(), value.data()),
                   std::string("cannot read HDF5 attribute: ") + name);
    return value.data();
}
} // namespace
#endif

bool hdf5Available() {
#ifdef SH1D_HAVE_HDF5
    return true;
#else
    return false;
#endif
}

void writeHdf5Field(const std::filesystem::path &path, const Parameters &parameters, double time,
                    std::uint64_t frame, const SpectralField &physical) {
#ifdef SH1D_HAVE_HDF5
    if (physical.size() != parameters.gridPoints)
        throw std::runtime_error("invalid HDF5 wavefunction size");
    if (!std::all_of(physical.begin(), physical.end(), [](Complex value) {
            return std::isfinite(value.real()) && std::isfinite(value.imag());
        }))
        throw std::runtime_error("refusing to write non-finite HDF5 wavefunction");

    File file(requireId(H5Fcreate(path.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT),
                        "cannot create HDF5 field file: " + path.string()));
    writeStringAttribute(file.get(), "format", "sh1d_wavefunction_v1");
    const auto count = static_cast<unsigned long long>(parameters.gridPoints);
    const auto savedFrame = static_cast<unsigned long long>(frame);
    writeScalarAttribute(file.get(), "grid_points", H5T_NATIVE_ULLONG, count);
    writeScalarAttribute(file.get(), "domain_length", H5T_NATIVE_DOUBLE, parameters.domainLength);
    writeScalarAttribute(file.get(), "time", H5T_NATIVE_DOUBLE, time);
    writeScalarAttribute(file.get(), "frame", H5T_NATIVE_ULLONG, savedFrame);

    const hsize_t dimensions[2] = {parameters.gridPoints, 2};
    DataSpace space(requireId(H5Screate_simple(2, dimensions, nullptr),
                              "cannot create HDF5 wavefunction space"));
    PropertyList properties;
    hid_t creationProperties = H5P_DEFAULT;
    if (parameters.hdf5CompressionLevel > 0) {
        properties = PropertyList(
            requireId(H5Pcreate(H5P_DATASET_CREATE), "cannot create HDF5 dataset properties"));
        const hsize_t chunks[2] = {std::min<hsize_t>(parameters.gridPoints, 4096), 2};
        requireSuccess(H5Pset_chunk(properties.get(), 2, chunks),
                       "cannot set HDF5 wavefunction chunks");
        requireSuccess(H5Pset_deflate(properties.get(), parameters.hdf5CompressionLevel),
                       "cannot enable HDF5 deflate compression");
        creationProperties = properties.get();
    }
    DataSet dataset(requireId(H5Dcreate2(file.get(), "/wavefunction", H5T_IEEE_F64LE, space.get(),
                                         H5P_DEFAULT, creationProperties, H5P_DEFAULT),
                              "cannot create HDF5 wavefunction dataset"));
    std::vector<double> values(2 * physical.size());
    for (std::size_t i = 0; i < physical.size(); ++i) {
        values[2 * i] = physical[i].real();
        values[2 * i + 1] = physical[i].imag();
    }
    requireSuccess(
        H5Dwrite(dataset.get(), H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, values.data()),
        "cannot write HDF5 wavefunction dataset");
#else
    (void)path;
    (void)parameters;
    (void)time;
    (void)frame;
    (void)physical;
    throw std::runtime_error("this build has no HDF5 support");
#endif
}

Hdf5Field readHdf5Field(const std::filesystem::path &path) {
#ifdef SH1D_HAVE_HDF5
    File file(requireId(H5Fopen(path.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT),
                        "cannot open HDF5 field file: " + path.string()));
    if (readStringAttribute(file.get(), "format") != "sh1d_wavefunction_v1")
        throw std::runtime_error("unsupported HDF5 wavefunction format: " + path.string());
    DataSet dataset(requireId(H5Dopen2(file.get(), "/wavefunction", H5P_DEFAULT),
                              "HDF5 file has no /wavefunction dataset: " + path.string()));
    DataSpace space(requireId(H5Dget_space(dataset.get()), "cannot inspect HDF5 wavefunction"));
    hsize_t dimensions[2]{};
    if (H5Sget_simple_extent_ndims(space.get()) != 2 ||
        H5Sget_simple_extent_dims(space.get(), dimensions, nullptr) < 0 || dimensions[1] != 2 ||
        dimensions[0] == 0 || dimensions[0] > std::numeric_limits<std::size_t>::max() / 2)
        throw std::runtime_error("HDF5 wavefunction must have shape [grid_points,2]");

    Hdf5Field result;
    result.gridPoints = static_cast<std::size_t>(dimensions[0]);
    const auto savedCount =
        readScalarAttribute<unsigned long long>(file.get(), "grid_points", H5T_NATIVE_ULLONG);
    if (savedCount != result.gridPoints)
        throw std::runtime_error("HDF5 dimensions disagree with file metadata");
    result.domainLength =
        readScalarAttribute<double>(file.get(), "domain_length", H5T_NATIVE_DOUBLE);
    result.time = readScalarAttribute<double>(file.get(), "time", H5T_NATIVE_DOUBLE);
    result.frame = static_cast<std::uint64_t>(
        readScalarAttribute<unsigned long long>(file.get(), "frame", H5T_NATIVE_ULLONG));
    if (!(result.domainLength > 0.0) || !std::isfinite(result.domainLength) ||
        !std::isfinite(result.time))
        throw std::runtime_error("HDF5 wavefunction metadata is invalid");

    std::vector<double> values(2 * result.gridPoints);
    requireSuccess(
        H5Dread(dataset.get(), H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, values.data()),
        "cannot read HDF5 wavefunction dataset");
    result.wavefunction.resize(result.gridPoints);
    for (std::size_t i = 0; i < result.gridPoints; ++i)
        result.wavefunction[i] = {values[2 * i], values[2 * i + 1]};
    if (!std::all_of(result.wavefunction.begin(), result.wavefunction.end(), [](Complex value) {
            return std::isfinite(value.real()) && std::isfinite(value.imag());
        }))
        throw std::runtime_error("HDF5 wavefunction contains non-finite values");
    return result;
#else
    (void)path;
    throw std::runtime_error("this build has no HDF5 support");
#endif
}
