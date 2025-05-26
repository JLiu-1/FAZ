//
// Four member functions are heavily based on QccPack:
//    http://qccpack.sourceforge.net/index.shtml
//  - void QccWAVCDF97AnalysisSymmetricEvenEven(double* signal, size_t signal_length);
//  - void QccWAVCDF97AnalysisSymmetricOddEven(double* signal, size_t signal_length);
//  - void QccWAVCDF97SynthesisSymmetricEvenEven(double* signal, size_t signal_length);
//  - void QccWAVCDF97SynthesisSymmetricOddEven(double* signal, size_t signal_length);
//

#ifndef CDF97_H
#define CDF97_H

#include "sperr_helper.h"

#include <cmath>
#include <algorithm>
#include <cassert>
#include <numeric>  // std::accumulate()
#include <type_traits>


namespace sperr {

class CDF97 {
 public:
  //
  // Input
  //
  // Note that copy_data() and take_data() effectively resets internal states of this class.
  template <typename T>
  auto copy_data(const T* buf, size_t len, dims_type dims) -> RTNType;
  auto take_data(vecd_type&& buf, dims_type dims) -> RTNType;

  //
  // Output
  //
  auto view_data() const -> const vecd_type&;
  auto release_data() -> vecd_type&&;
  auto get_dims() const -> std::array<size_t, 3>;  // In 2D case, the 3rd value equals 1.

  //
  // Action items
  //
  void dwt1d();
  void dwt2d();
  void dwt3d();
  void idwt2d();
  void idwt1d();
  void idwt3d();

  //
  // Multi-resolution reconstruction
  //

  // If multi-resolution is supported (determined by `sperr::available_resolutions()`), then
  //    it returns all the coarsened volumes, which are placed in the same order of resolutions
  //    returned by `sperr::available_resolutions()`. The native resolution reconstruction should
  //    still be retrieved by the `view_data()` or `release_data()` functions.
  //    If multi-resolution is not supported, then it simply returns an empty vector, with the
  //    decompression still performed, and the native resolution reconstruction ready.
  [[nodiscard]] auto idwt2d_multi_res() -> std::vector<vecd_type>;
  void idwt3d_multi_res(std::vector<vecd_type>&);

 private:
  using itd_type = vecd_type::iterator;
  using citd_type = vecd_type::const_iterator;

  //
  // Private methods helping DWT.
  //

  // Multiple levels of 1D DWT/IDWT on a given array of length array_len.
  void m_dwt1d(itd_type array, size_t array_len, size_t num_of_xforms);
  void m_idwt1d(itd_type array, size_t array_len, size_t num_of_xforms);

  // Multiple levels of 2D DWT/IDWT on a given plane by repeatedly invoking
  // m_dwt2d_one_level(). The plane has a dimension (len_xy[0], len_xy[1]).
  void m_dwt2d(itd_type plane, std::array<size_t, 2> len_xy, size_t num_of_xforms);
  void m_idwt2d(itd_type plane, std::array<size_t, 2> len_xy, size_t num_of_xforms);

  // Perform one level of interleaved 3D dwt/idwt on a given volume (m_dims),
  // specifically on its top left (len_xyz) subset.
  void m_dwt3d_one_level(itd_type vol, std::array<size_t, 3> len_xyz);
  void m_idwt3d_one_level(itd_type vol, std::array<size_t, 3> len_xyz);

  // Perform one level of 2D dwt/idwt on a given plane (m_dims),
  // specifically on its top left (len_xy) subset.
  void m_dwt2d_one_level(itd_type plane, std::array<size_t, 2> len_xy);
  void m_idwt2d_one_level(itd_type plane, std::array<size_t, 2> len_xy);

  // Perform one level of 1D dwt/idwt on a given array (array_len).
  // A buffer space (tmp_buf) should be passed in for
  // this method to work on with length at least 2*array_len.
  void m_dwt1d_one_level(itd_type array, size_t array_len);
  void m_idwt1d_one_level(itd_type array, size_t array_len);

  // Separate even and odd indexed elements to be at the front and back of the dest array.
  // Note 1: sufficient memory space should be allocated by the caller.
  // Note 2: two versions for even and odd length input.
  void m_gather_even(citd_type begin, citd_type end, itd_type dest) const;
  void m_gather_odd(citd_type begin, citd_type end, itd_type dest) const;

  // Interleave low and high pass elements to be at even and odd positions of the dest array.
  // Note 1: sufficient memory space should be allocated by the caller.
  // Note 2: two versions for even and odd length input.
  void m_scatter_even(citd_type begin, citd_type end, itd_type dest) const;
  void m_scatter_odd(citd_type begin, citd_type end, itd_type dest) const;

  // Two flavors of 3D transforms.
  // They should be invoked by the `dwt3d()` and `idwt3d()` public methods, not users, though.
  void m_dwt3d_wavelet_packet();
  void m_idwt3d_wavelet_packet();
  void m_dwt3d_dyadic(size_t num_xforms);
  void m_idwt3d_dyadic(size_t num_xforms);

  // Extract a sub-slice/sub-volume starting with the same origin of the full slice/volume.
  // It is UB if `subdims` exceeds the full dimension (`m_dims`).
  // It is UB if `dst` does not point to a big enough space.
  auto m_sub_slice(std::array<size_t, 2> subdims) const -> vecd_type;
  void m_sub_volume(dims_type subdims, itd_type dst) const;

  //
  // Methods from QccPack, so keep their original names, interface, and the use of raw pointers.
  //
  void QccWAVCDF97AnalysisSymmetricEvenEven(double* signal, size_t signal_length);
  void QccWAVCDF97AnalysisSymmetricOddEven(double* signal, size_t signal_length);
  void QccWAVCDF97SynthesisSymmetricEvenEven(double* signal, size_t signal_length);
  void QccWAVCDF97SynthesisSymmetricOddEven(double* signal, size_t signal_length);

  //
  // Private data members
  //
  vecd_type m_data_buf;          // Holds the entire input data.
  dims_type m_dims = {0, 0, 0};  // Dimension of the data volume

  // Temporary buffers that are big enough for any (1D column * 2) or any 2D
  // slice. Note: `m_qcc_buf` should be used by m_***_one_level() functions and
  // should not be used by higher-level functions. `m_slice_buf` is only used by
  // wavelet-packet transforms.
  vecd_type m_qcc_buf;
  vecd_type m_slice_buf;

  //
  // Note on the coefficients and constants:
  // The ones from QccPack are slightly different from what's described in the
  // lifting scheme paper: Pg19 of "FACTORING WAVELET TRANSFORMS INTO LIFTING STEPS,"
  // DAUBECHIES and SWELDEN.  (https://9p.io/who/wim/papers/factor/factor.pdf)
  // JasPer, OpenJPEG, and FFMpeg use coefficients closer to the paper.
  // The filter bank coefficients (h[] array) are from "Biorthogonal Bases of
  // Compactly Supported Wavelets," by Cohen et al., Page 551.
  // (https://services.math.duke.edu/~ingrid/publications/CPAM_1992_p485.pdf)
  //

  // Paper coefficients
  const std::array<double, 5> h = {0.602949018236, 0.266864118443, -0.078223266529, -0.016864118443,
                                   0.026748757411};
  const double r0 = h[0] - 2.0 * h[4] * h[1] / h[3];
  const double r1 = h[2] - h[4] - h[4] * h[1] / h[3];
  const double s0 = h[1] - h[3] - h[3] * r0 / r1;
  const double t0 = h[0] - 2.0 * (h[2] - h[4]);
  const double ALPHA = h[4] / h[3];
  const double BETA = h[3] / r1;
  const double GAMMA = r1 / s0;
  const double DELTA = s0 / t0;
  const double EPSILON = std::sqrt(2.0) * t0;
  const double INV_EPSILON = 1.0 / EPSILON;

  // QccPack coefficients
  //
  // const double ALPHA   = -1.58615986717275;
  // const double BETA    = -0.05297864003258;
  // const double GAMMA   =  0.88293362717904;
  // const double DELTA   =  0.44350482244527;
  // const double EPSILON =  1.14960430535816;
  //
};

};  // namespace sperr

template <typename T>
auto sperr::CDF97::copy_data(const T* data, size_t len, dims_type dims) -> RTNType
{
  static_assert(std::is_floating_point<T>::value, "!! Only floating point values are supported !!");
  if (len != dims[0] * dims[1] * dims[2])
    return RTNType::WrongLength;

  m_data_buf.resize(len);
  std::copy(data, data + len, m_data_buf.begin());

  m_dims = dims;

  auto max_col = std::max(std::max(dims[0], dims[1]), dims[2]);
  if (max_col * 2 > m_qcc_buf.size())
    m_qcc_buf.resize(std::max(m_qcc_buf.size(), max_col) * 2);

  auto max_slice = std::max(std::max(dims[0] * dims[1], dims[0] * dims[2]), dims[1] * dims[2]);
  if (max_slice > m_slice_buf.size())
    m_slice_buf.resize(std::max(m_slice_buf.size() * 2, max_slice));

  return RTNType::Good;
}
template auto sperr::CDF97::copy_data(const float*, size_t, dims_type) -> RTNType;
template auto sperr::CDF97::copy_data(const double*, size_t, dims_type) -> RTNType;

auto sperr::CDF97::take_data(vecd_type&& buf, dims_type dims) -> RTNType
{
  if (buf.size() != dims[0] * dims[1] * dims[2])
    return RTNType::WrongLength;

  m_data_buf = std::move(buf);
  m_dims = dims;

  auto max_col = std::max(std::max(dims[0], dims[1]), dims[2]);
  if (max_col * 2 > m_qcc_buf.size())
    m_qcc_buf.resize(std::max(m_qcc_buf.size(), max_col) * 2);

  auto max_slice = std::max(std::max(dims[0] * dims[1], dims[0] * dims[2]), dims[1] * dims[2]);
  if (max_slice > m_slice_buf.size())
    m_slice_buf.resize(std::max(m_slice_buf.size() * 2, max_slice));

  return RTNType::Good;
}

auto sperr::CDF97::view_data() const -> const vecd_type&
{
  return m_data_buf;
}

auto sperr::CDF97::release_data() -> vecd_type&&
{
  return std::move(m_data_buf);
}

auto sperr::CDF97::get_dims() const -> std::array<size_t, 3>
{
  return m_dims;
}

void sperr::CDF97::dwt1d()
{
  auto num_xforms = sperr::num_of_xforms(m_dims[0]);
  m_dwt1d(m_data_buf.begin(), m_data_buf.size(), num_xforms);
}

void sperr::CDF97::idwt1d()
{
  auto num_xforms = sperr::num_of_xforms(m_dims[0]);
  m_idwt1d(m_data_buf.begin(), m_data_buf.size(), num_xforms);
}

void sperr::CDF97::dwt2d()
{
  auto xy = sperr::num_of_xforms(std::min(m_dims[0], m_dims[1]));
  m_dwt2d(m_data_buf.begin(), {m_dims[0], m_dims[1]}, xy);
}

void sperr::CDF97::idwt2d()
{
  auto xy = sperr::num_of_xforms(std::min(m_dims[0], m_dims[1]));
  m_idwt2d(m_data_buf.begin(), {m_dims[0], m_dims[1]}, xy);
}

auto sperr::CDF97::idwt2d_multi_res() -> std::vector<vecd_type>
{
  const auto xy = sperr::num_of_xforms(std::min(m_dims[0], m_dims[1]));
  auto ret = std::vector<vecd_type>();

  if (xy > 0) {
    ret.reserve(xy);
    for (size_t lev = xy; lev > 0; lev--) {
      auto [x, xd] = sperr::calc_approx_detail_len(m_dims[0], lev);
      auto [y, yd] = sperr::calc_approx_detail_len(m_dims[1], lev);
      ret.emplace_back(m_sub_slice({x, y}));
      m_idwt2d_one_level(m_data_buf.begin(), {x + xd, y + yd});
    }
  }

  return ret;
}

void sperr::CDF97::dwt3d()
{
  auto dyadic = sperr::can_use_dyadic(m_dims);
  if (dyadic)
    m_dwt3d_dyadic(*dyadic);
  else
    m_dwt3d_wavelet_packet();
}

void sperr::CDF97::idwt3d()
{
  auto dyadic = sperr::can_use_dyadic(m_dims);
  if (dyadic)
    m_idwt3d_dyadic(*dyadic);
  else
    m_idwt3d_wavelet_packet();
}

void sperr::CDF97::idwt3d_multi_res(std::vector<vecd_type>& h)
{
  auto dyadic = sperr::can_use_dyadic(m_dims);

  if (dyadic) {
    h.resize(*dyadic);
    for (size_t lev = *dyadic; lev > 0; lev--) {
      auto [x, xd] = sperr::calc_approx_detail_len(m_dims[0], lev);
      auto [y, yd] = sperr::calc_approx_detail_len(m_dims[1], lev);
      auto [z, zd] = sperr::calc_approx_detail_len(m_dims[2], lev);
      auto& buf = h[*dyadic - lev];
      buf.resize(x * y * z);
      m_sub_volume({x, y, z}, buf.begin());
      m_idwt3d_one_level(m_data_buf.begin(), {x + xd, y + yd, z + zd});
    }
  }
  else
    m_idwt3d_wavelet_packet();
}

void sperr::CDF97::m_dwt3d_wavelet_packet()
{
  /*
   *             Z
   *            /
   *           /
   *          /________
   *         /       /|
   *        /       / |
   *     0 |-------/-------> X
   *       |       |  |
   *       |       |  /
   *       |       | /
   *       |_______|/
   *       |
   *       |
   *       Y
   */

  const size_t plane_size_xy = m_dims[0] * m_dims[1];

  // First transform along the Z dimension
  //
  const auto num_xforms_z = sperr::num_of_xforms(m_dims[2]);

  for (size_t y = 0; y < m_dims[1]; y++) {
    const auto y_offset = y * m_dims[0];

    // Re-arrange values of one XZ slice so that they form many z_columns
    for (size_t z = 0; z < m_dims[2]; z++) {
      const auto cube_start_idx = z * plane_size_xy + y_offset;
      for (size_t x = 0; x < m_dims[0]; x++)
        m_slice_buf[z + x * m_dims[2]] = m_data_buf[cube_start_idx + x];
    }

    // DWT1D on every z_column
    for (size_t x = 0; x < m_dims[0]; x++)
      m_dwt1d(m_slice_buf.begin() + x * m_dims[2], m_dims[2], num_xforms_z);

    // Put back values of the z_columns to the cube
    for (size_t z = 0; z < m_dims[2]; z++) {
      const auto cube_start_idx = z * plane_size_xy + y_offset;
      for (size_t x = 0; x < m_dims[0]; x++)
        m_data_buf[cube_start_idx + x] = m_slice_buf[z + x * m_dims[2]];
    }
  }

  // Second transform each plane
  //
  const auto num_xforms_xy = sperr::num_of_xforms(std::min(m_dims[0], m_dims[1]));

  for (size_t z = 0; z < m_dims[2]; z++) {
    const size_t offset = plane_size_xy * z;
    m_dwt2d(m_data_buf.begin() + offset, {m_dims[0], m_dims[1]}, num_xforms_xy);
  }
}

void sperr::CDF97::m_idwt3d_wavelet_packet()
{
  const size_t plane_size_xy = m_dims[0] * m_dims[1];

  // First, inverse transform each plane
  //
  auto num_xforms_xy = sperr::num_of_xforms(std::min(m_dims[0], m_dims[1]));
  for (size_t i = 0; i < m_dims[2]; i++) {
    const size_t offset = plane_size_xy * i;
    m_idwt2d(m_data_buf.begin() + offset, {m_dims[0], m_dims[1]}, num_xforms_xy);
  }

  /*
   * Second, inverse transform along the Z dimension
   *
   *             Z
   *            /
   *           /
   *          /________
   *         /       /|
   *        /       / |
   *     0 |-------/-------> X
   *       |       |  |
   *       |       |  /
   *       |       | /
   *       |_______|/
   *       |
   *       |
   *       Y
   */

  // Process one XZ slice at a time
  //
  const auto num_xforms_z = sperr::num_of_xforms(m_dims[2]);
  for (size_t y = 0; y < m_dims[1]; y++) {
    const auto y_offset = y * m_dims[0];

    // Re-arrange values on one slice so that they form many z_columns
    for (size_t z = 0; z < m_dims[2]; z++) {
      const auto cube_start_idx = z * plane_size_xy + y_offset;
      for (size_t x = 0; x < m_dims[0]; x++)
        m_slice_buf[z + x * m_dims[2]] = m_data_buf[cube_start_idx + x];
    }

    // IDWT1D on every z_column
    for (size_t x = 0; x < m_dims[0]; x++)
      m_idwt1d(m_slice_buf.begin() + x * m_dims[2], m_dims[2], num_xforms_z);

    // Put back values from the z_columns to the cube
    for (size_t z = 0; z < m_dims[2]; z++) {
      const auto cube_start_idx = z * plane_size_xy + y_offset;
      for (size_t x = 0; x < m_dims[0]; x++)
        m_data_buf[cube_start_idx + x] = m_slice_buf[z + x * m_dims[2]];
    }
  }
}

void sperr::CDF97::m_dwt3d_dyadic(size_t num_xforms)
{
  for (size_t lev = 0; lev < num_xforms; lev++) {
    auto [x, xd] = sperr::calc_approx_detail_len(m_dims[0], lev);
    auto [y, yd] = sperr::calc_approx_detail_len(m_dims[1], lev);
    auto [z, zd] = sperr::calc_approx_detail_len(m_dims[2], lev);
    m_dwt3d_one_level(m_data_buf.begin(), {x, y, z});
  }
}

void sperr::CDF97::m_idwt3d_dyadic(size_t num_xforms)
{
  for (size_t lev = num_xforms; lev > 0; lev--) {
    auto [x, xd] = sperr::calc_approx_detail_len(m_dims[0], lev - 1);
    auto [y, yd] = sperr::calc_approx_detail_len(m_dims[1], lev - 1);
    auto [z, zd] = sperr::calc_approx_detail_len(m_dims[2], lev - 1);
    m_idwt3d_one_level(m_data_buf.begin(), {x, y, z});
  }
}

//
// Private Methods
//
void sperr::CDF97::m_dwt1d(itd_type array, size_t array_len, size_t num_of_lev)
{
  for (size_t lev = 0; lev < num_of_lev; lev++) {
    auto [x, xd] = sperr::calc_approx_detail_len(array_len, lev);
    m_dwt1d_one_level(array, x);
  }
}

void sperr::CDF97::m_idwt1d(itd_type array, size_t array_len, size_t num_of_lev)
{
  for (size_t lev = num_of_lev; lev > 0; lev--) {
    auto [x, xd] = sperr::calc_approx_detail_len(array_len, lev - 1);
    m_idwt1d_one_level(array, x);
  }
}

void sperr::CDF97::m_dwt2d(itd_type plane, std::array<size_t, 2> len_xy, size_t num_of_lev)
{
  for (size_t lev = 0; lev < num_of_lev; lev++) {
    auto [x, xd] = sperr::calc_approx_detail_len(len_xy[0], lev);
    auto [y, yd] = sperr::calc_approx_detail_len(len_xy[1], lev);
    m_dwt2d_one_level(plane, {x, y});
  }
}

void sperr::CDF97::m_idwt2d(itd_type plane, std::array<size_t, 2> len_xy, size_t num_of_lev)
{
  for (size_t lev = num_of_lev; lev > 0; lev--) {
    auto [x, xd] = sperr::calc_approx_detail_len(len_xy[0], lev - 1);
    auto [y, yd] = sperr::calc_approx_detail_len(len_xy[1], lev - 1);
    m_idwt2d_one_level(plane, {x, y});
  }
}

void sperr::CDF97::m_dwt1d_one_level(itd_type array, size_t array_len)
{
  std::copy(array, array + array_len, m_qcc_buf.begin());
  if (array_len % 2 == 0) {
    this->QccWAVCDF97AnalysisSymmetricEvenEven(m_qcc_buf.data(), array_len);
    m_gather_even(m_qcc_buf.cbegin(), m_qcc_buf.cbegin() + array_len, array);
  }
  else {
    this->QccWAVCDF97AnalysisSymmetricOddEven(m_qcc_buf.data(), array_len);
    m_gather_odd(m_qcc_buf.cbegin(), m_qcc_buf.cbegin() + array_len, array);
  }
}

void sperr::CDF97::m_idwt1d_one_level(itd_type array, size_t array_len)
{
  if (array_len % 2 == 0) {
    m_scatter_even(array, array + array_len, m_qcc_buf.begin());
    this->QccWAVCDF97SynthesisSymmetricEvenEven(m_qcc_buf.data(), array_len);
  }
  else {
    m_scatter_odd(array, array + array_len, m_qcc_buf.begin());
    this->QccWAVCDF97SynthesisSymmetricOddEven(m_qcc_buf.data(), array_len);
  }
  std::copy(m_qcc_buf.cbegin(), m_qcc_buf.cbegin() + array_len, array);
}

void sperr::CDF97::m_dwt2d_one_level(itd_type plane, std::array<size_t, 2> len_xy)
{
  // Note: here we call low-level functions (Qcc*()) instead of
  // m_dwt1d_one_level() because we want to have only one even/odd test at the outer loop.

  const auto max_len = std::max(len_xy[0], len_xy[1]);
  const auto beg = m_qcc_buf.begin();
  const auto beg2 = beg + max_len;

  // First, perform DWT along X for every row
  if (len_xy[0] % 2 == 0) {
    for (size_t i = 0; i < len_xy[1]; i++) {
      auto pos = plane + i * m_dims[0];
      std::copy(pos, pos + len_xy[0], beg);
      this->QccWAVCDF97AnalysisSymmetricEvenEven(m_qcc_buf.data(), len_xy[0]);
      m_gather_even(beg, beg + len_xy[0], pos);
    }
  }
  else  // Odd length
  {
    for (size_t i = 0; i < len_xy[1]; i++) {
      auto pos = plane + i * m_dims[0];
      std::copy(pos, pos + len_xy[0], beg);
      this->QccWAVCDF97AnalysisSymmetricOddEven(m_qcc_buf.data(), len_xy[0]);
      m_gather_odd(beg, beg + len_xy[0], pos);
    }
  }

  // Second, perform DWT along Y for every column
  // Note, I've tested that up to 1024^2 planes it is actually slightly slower
  // to transpose the plane and then perform the transforms. This was consistent
  // on both a MacBook and a RaspberryPi 3. Note2, I've tested transpose again
  // on an X86 linux machine using gcc, clang, and pgi. Again the difference is
  // either indistinguishable, or the current implementation has a slight edge.

  if (len_xy[1] % 2 == 0) {
    for (size_t x = 0; x < len_xy[0]; x++) {
      for (size_t y = 0; y < len_xy[1]; y++)
        m_qcc_buf[y] = *(plane + y * m_dims[0] + x);
      this->QccWAVCDF97AnalysisSymmetricEvenEven(m_qcc_buf.data(), len_xy[1]);
      m_gather_even(beg, beg + len_xy[1], beg2);
      for (size_t y = 0; y < len_xy[1]; y++)
        *(plane + y * m_dims[0] + x) = *(beg2 + y);
    }
  }
  else  // Odd length
  {
    for (size_t x = 0; x < len_xy[0]; x++) {
      for (size_t y = 0; y < len_xy[1]; y++)
        m_qcc_buf[y] = *(plane + y * m_dims[0] + x);
      this->QccWAVCDF97AnalysisSymmetricOddEven(m_qcc_buf.data(), len_xy[1]);
      m_gather_odd(beg, beg + len_xy[1], beg2);
      for (size_t y = 0; y < len_xy[1]; y++)
        *(plane + y * m_dims[0] + x) = *(beg2 + y);
    }
  }
}

void sperr::CDF97::m_idwt2d_one_level(itd_type plane, std::array<size_t, 2> len_xy)
{
  const auto max_len = std::max(len_xy[0], len_xy[1]);
  const auto beg = m_qcc_buf.begin();  // First half of the buffer
  const auto beg2 = beg + max_len;     // Second half of the buffer

  // First, perform IDWT along Y for every column
  if (len_xy[1] % 2 == 0) {
    for (size_t x = 0; x < len_xy[0]; x++) {
      for (size_t y = 0; y < len_xy[1]; y++)
        m_qcc_buf[y] = *(plane + y * m_dims[0] + x);
      m_scatter_even(beg, beg + len_xy[1], beg2);
      this->QccWAVCDF97SynthesisSymmetricEvenEven(m_qcc_buf.data() + max_len, len_xy[1]);
      for (size_t y = 0; y < len_xy[1]; y++)
        *(plane + y * m_dims[0] + x) = *(beg2 + y);
    }
  }
  else  // Odd length
  {
    for (size_t x = 0; x < len_xy[0]; x++) {
      for (size_t y = 0; y < len_xy[1]; y++)
        m_qcc_buf[y] = *(plane + y * m_dims[0] + x);
      m_scatter_odd(beg, beg + len_xy[1], beg2);
      this->QccWAVCDF97SynthesisSymmetricOddEven(m_qcc_buf.data() + max_len, len_xy[1]);
      for (size_t y = 0; y < len_xy[1]; y++)
        *(plane + y * m_dims[0] + x) = *(beg2 + y);
    }
  }

  // Second, perform IDWT along X for every row
  if (len_xy[0] % 2 == 0) {
    for (size_t i = 0; i < len_xy[1]; i++) {
      auto pos = plane + i * m_dims[0];
      m_scatter_even(pos, pos + len_xy[0], beg);
      this->QccWAVCDF97SynthesisSymmetricEvenEven(m_qcc_buf.data(), len_xy[0]);
      std::copy(beg, beg + len_xy[0], pos);
    }
  }
  else  // Odd length
  {
    for (size_t i = 0; i < len_xy[1]; i++) {
      auto pos = plane + i * m_dims[0];
      m_scatter_odd(pos, pos + len_xy[0], beg);
      this->QccWAVCDF97SynthesisSymmetricOddEven(m_qcc_buf.data(), len_xy[0]);
      std::copy(beg, beg + len_xy[0], pos);
    }
  }
}

void sperr::CDF97::m_dwt3d_one_level(itd_type vol, std::array<size_t, 3> len_xyz)
{
  // First, do one level of transform on all XY planes.
  const auto plane_size_xy = m_dims[0] * m_dims[1];
  for (size_t z = 0; z < len_xyz[2]; z++) {
    const size_t offset = plane_size_xy * z;
    m_dwt2d_one_level(vol + offset, {len_xyz[0], len_xyz[1]});
  }

  const auto beg = m_qcc_buf.begin();  // First half of the buffer
  const auto beg2 = beg + len_xyz[2];  // Second half of the buffer

  // Second, do one level of transform on all Z columns.  Strategy:
  // 1) extract a Z column to buffer space `m_qcc_buf`
  // 2) use appropriate even/odd Qcc*** function to transform it
  // 3) gather coefficients from `m_qcc_buf` to the second half of `m_qcc_buf`
  // 4) put the Z column back to their locations as a Z column.

  if (len_xyz[2] % 2 == 0) {  // Even length
    for (size_t y = 0; y < len_xyz[1]; y++) {
      for (size_t x = 0; x < len_xyz[0]; x++) {
        const size_t xy_offset = y * m_dims[0] + x;
        // Step 1
        for (size_t z = 0; z < len_xyz[2]; z++)
          m_qcc_buf[z] = m_data_buf[z * plane_size_xy + xy_offset];
        // Step 2
        this->QccWAVCDF97AnalysisSymmetricEvenEven(m_qcc_buf.data(), len_xyz[2]);
        // Step 3
        m_gather_even(beg, beg2, beg2);
        // Step 4
        for (size_t z = 0; z < len_xyz[2]; z++)
          m_data_buf[z * plane_size_xy + xy_offset] = *(beg2 + z);
      }
    }
  }
  else {  // Odd length
    for (size_t y = 0; y < len_xyz[1]; y++) {
      for (size_t x = 0; x < len_xyz[0]; x++) {
        const size_t xy_offset = y * m_dims[0] + x;
        // Step 1
        for (size_t z = 0; z < len_xyz[2]; z++)
          m_qcc_buf[z] = m_data_buf[z * plane_size_xy + xy_offset];
        // Step 2
        this->QccWAVCDF97AnalysisSymmetricOddEven(m_qcc_buf.data(), len_xyz[2]);
        // Step 3
        m_gather_odd(beg, beg2, beg2);
        // Step 4
        for (size_t z = 0; z < len_xyz[2]; z++)
          m_data_buf[z * plane_size_xy + xy_offset] = *(beg2 + z);
      }
    }
  }
}

void sperr::CDF97::m_idwt3d_one_level(itd_type vol, std::array<size_t, 3> len_xyz)
{
  const auto plane_size_xy = m_dims[0] * m_dims[1];
  const auto beg = m_qcc_buf.begin();  // First half of the buffer
  const auto beg2 = beg + len_xyz[2];  // Second half of the buffer

  // First, do one level of inverse transform on all Z columns.  Strategy:
  // 1) extract a Z column to buffer space `m_qcc_buf`
  // 2) scatter coefficients from `m_qcc_buf` to the second half of `m_qcc_buf`
  // 3) use appropriate even/odd Qcc*** function to transform it
  // 4) put the Z column back to their locations as a Z column.

  if (len_xyz[2] % 2 == 0) {
    for (size_t y = 0; y < len_xyz[1]; y++) {
      for (size_t x = 0; x < len_xyz[0]; x++) {
        const size_t xy_offset = y * m_dims[0] + x;
        // Step 1
        for (size_t z = 0; z < len_xyz[2]; z++)
          m_qcc_buf[z] = m_data_buf[z * plane_size_xy + xy_offset];
        // Step 2
        m_scatter_even(beg, beg2, beg2);
        // Step 3
        this->QccWAVCDF97SynthesisSymmetricEvenEven(m_qcc_buf.data() + len_xyz[2], len_xyz[2]);
        // Step 4
        for (size_t z = 0; z < len_xyz[2]; z++)
          m_data_buf[z * plane_size_xy + xy_offset] = *(beg2 + z);
      }
    }
  }
  else {
    for (size_t y = 0; y < len_xyz[1]; y++) {
      for (size_t x = 0; x < len_xyz[0]; x++) {
        const size_t xy_offset = y * m_dims[0] + x;
        // Step 1
        for (size_t z = 0; z < len_xyz[2]; z++)
          m_qcc_buf[z] = m_data_buf[z * plane_size_xy + xy_offset];
        // Step 2
        m_scatter_odd(beg, beg2, beg2);
        // Step 3
        this->QccWAVCDF97SynthesisSymmetricOddEven(m_qcc_buf.data() + len_xyz[2], len_xyz[2]);
        // Step 4
        for (size_t z = 0; z < len_xyz[2]; z++)
          m_data_buf[z * plane_size_xy + xy_offset] = *(beg2 + z);
      }
    }
  }

  // Second, do one level of inverse transform on all XY planes.
  for (size_t z = 0; z < len_xyz[2]; z++) {
    const size_t offset = plane_size_xy * z;
    m_idwt2d_one_level(vol + offset, {len_xyz[0], len_xyz[1]});
  }
}

void sperr::CDF97::m_gather_even(citd_type begin, citd_type end, itd_type dest) const
{
  auto len = end - begin;
  assert(len % 2 == 0);  // This function specifically for even length input
  size_t low_count = len / 2, high_count = len / 2;
  for (size_t i = 0; i < low_count; i++) {
    *dest = *(begin + i * 2);
    ++dest;
  }
  for (size_t i = 0; i < high_count; i++) {
    *dest = *(begin + i * 2 + 1);
    ++dest;
  }
}

void sperr::CDF97::m_gather_odd(citd_type begin, citd_type end, itd_type dest) const
{
  auto len = end - begin;
  assert(len % 2 == 1);  // This function specifically for odd length input
  size_t low_count = len / 2 + 1, high_count = len / 2;
  for (size_t i = 0; i < low_count; i++) {
    *dest = *(begin + i * 2);
    ++dest;
  }
  for (size_t i = 0; i < high_count; i++) {
    *dest = *(begin + i * 2 + 1);
    ++dest;
  }
}

void sperr::CDF97::m_scatter_even(citd_type begin, citd_type end, itd_type dest) const
{
  auto len = end - begin;
  assert(len % 2 == 0);  // This function specifically for even length input
  size_t low_count = len / 2, high_count = len / 2;
  for (size_t i = 0; i < low_count; i++) {
    *(dest + i * 2) = *begin;
    ++begin;
  }
  for (size_t i = 0; i < high_count; i++) {
    *(dest + i * 2 + 1) = *begin;
    ++begin;
  }
}

void sperr::CDF97::m_scatter_odd(citd_type begin, citd_type end, itd_type dest) const
{
  auto len = end - begin;
  assert(len % 2 == 1);  // This function specifically for odd length input
  size_t low_count = len / 2 + 1, high_count = len / 2;
  for (size_t i = 0; i < low_count; i++) {
    *(dest + i * 2) = *begin;
    ++begin;
  }
  for (size_t i = 0; i < high_count; i++) {
    *(dest + i * 2 + 1) = *begin;
    ++begin;
  }
}

auto sperr::CDF97::m_sub_slice(std::array<size_t, 2> subdims) const -> vecd_type
{
  assert(subdims[0] <= m_dims[0] && subdims[1] <= m_dims[1]);

  auto ret = vecd_type(subdims[0] * subdims[1]);
  auto dst = ret.begin();
  for (size_t y = 0; y < subdims[1]; y++) {
    auto beg = m_data_buf.begin() + y * m_dims[0];
    std::copy(beg, beg + subdims[0], dst);
    dst += subdims[0];
  }

  return ret;
}

void sperr::CDF97::m_sub_volume(dims_type subdims, itd_type dst) const
{
  assert(subdims[0] <= m_dims[0] && subdims[1] <= m_dims[1] && subdims[2] <= m_dims[2]);

  const auto slice_len = m_dims[0] * m_dims[1];
  for (size_t z = 0; z < subdims[2]; z++) {
    for (size_t y = 0; y < subdims[1]; y++) {
      auto beg = m_data_buf.begin() + z * slice_len + y * m_dims[0];
      std::copy(beg, beg + subdims[0], dst);
      dst += subdims[0];
    }
  }
}

//
// Methods from QccPack
//
void sperr::CDF97::QccWAVCDF97AnalysisSymmetricEvenEven(double* signal, size_t signal_length)
{
  for (size_t i = 1; i < signal_length - 2; i += 2)
    signal[i] += ALPHA * (signal[i - 1] + signal[i + 1]);

  signal[signal_length - 1] += 2.0 * ALPHA * signal[signal_length - 2];

  signal[0] += 2.0 * BETA * signal[1];

  for (size_t i = 2; i < signal_length; i += 2)
    signal[i] += BETA * (signal[i + 1] + signal[i - 1]);

  for (size_t i = 1; i < signal_length - 2; i += 2)
    signal[i] += GAMMA * (signal[i - 1] + signal[i + 1]);

  signal[signal_length - 1] += 2.0 * GAMMA * signal[signal_length - 2];

  signal[0] = EPSILON * (signal[0] + 2.0 * DELTA * signal[1]);

  for (size_t i = 2; i < signal_length; i += 2)
    signal[i] = EPSILON * (signal[i] + DELTA * (signal[i + 1] + signal[i - 1]));

  for (size_t i = 1; i < signal_length; i += 2)
    signal[i] *= -INV_EPSILON;
}

void sperr::CDF97::QccWAVCDF97SynthesisSymmetricEvenEven(double* signal, size_t signal_length)
{
  for (size_t i = 1; i < signal_length; i += 2)
    signal[i] *= (-EPSILON);

  signal[0] = signal[0] * INV_EPSILON - 2.0 * DELTA * signal[1];

  for (size_t i = 2; i < signal_length; i += 2)
    signal[i] = signal[i] * INV_EPSILON - DELTA * (signal[i + 1] + signal[i - 1]);

  for (size_t i = 1; i < signal_length - 2; i += 2)
    signal[i] -= GAMMA * (signal[i - 1] + signal[i + 1]);

  signal[signal_length - 1] -= 2.0 * GAMMA * signal[signal_length - 2];

  signal[0] -= 2.0 * BETA * signal[1];

  for (size_t i = 2; i < signal_length; i += 2)
    signal[i] -= BETA * (signal[i + 1] + signal[i - 1]);

  for (size_t i = 1; i < signal_length - 2; i += 2)
    signal[i] -= ALPHA * (signal[i - 1] + signal[i + 1]);

  signal[signal_length - 1] -= 2.0 * ALPHA * signal[signal_length - 2];
}

void sperr::CDF97::QccWAVCDF97SynthesisSymmetricOddEven(double* signal, size_t signal_length)
{
  for (size_t i = 1; i < signal_length - 1; i += 2)
    signal[i] *= (-EPSILON);

  signal[0] = signal[0] * INV_EPSILON - 2.0 * DELTA * signal[1];

  for (size_t i = 2; i < signal_length - 2; i += 2)
    signal[i] = signal[i] * INV_EPSILON - DELTA * (signal[i + 1] + signal[i - 1]);

  signal[signal_length - 1] =
      signal[signal_length - 1] * INV_EPSILON - 2.0 * DELTA * signal[signal_length - 2];

  for (size_t i = 1; i < signal_length - 1; i += 2)
    signal[i] -= GAMMA * (signal[i - 1] + signal[i + 1]);

  signal[0] -= 2.0 * BETA * signal[1];

  for (size_t i = 2; i < signal_length - 2; i += 2)
    signal[i] -= BETA * (signal[i + 1] + signal[i - 1]);

  signal[signal_length - 1] -= 2.0 * BETA * signal[signal_length - 2];

  for (size_t i = 1; i < signal_length - 1; i += 2)
    signal[i] -= ALPHA * (signal[i - 1] + signal[i + 1]);
}

void sperr::CDF97::QccWAVCDF97AnalysisSymmetricOddEven(double* signal, size_t signal_length)
{
  for (size_t i = 1; i < signal_length - 1; i += 2)
    signal[i] += ALPHA * (signal[i - 1] + signal[i + 1]);

  signal[0] += 2.0 * BETA * signal[1];

  for (size_t i = 2; i < signal_length - 2; i += 2)
    signal[i] += BETA * (signal[i + 1] + signal[i - 1]);

  signal[signal_length - 1] += 2.0 * BETA * signal[signal_length - 2];

  for (size_t i = 1; i < signal_length - 1; i += 2)
    signal[i] += GAMMA * (signal[i - 1] + signal[i + 1]);

  signal[0] = EPSILON * (signal[0] + 2.0 * DELTA * signal[1]);

  for (size_t i = 2; i < signal_length - 2; i += 2)
    signal[i] = EPSILON * (signal[i] + DELTA * (signal[i + 1] + signal[i - 1]));

  signal[signal_length - 1] =
      EPSILON * (signal[signal_length - 1] + 2.0 * DELTA * signal[signal_length - 2]);

  for (size_t i = 1; i < signal_length - 1; i += 2)
    signal[i] *= (-INV_EPSILON);
}




#endif
