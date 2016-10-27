#include "Triode3Filter.h"

#include <boost/math/constants/constants.hpp>
#include <cassert>

namespace ATK {

template <typename DataType>
Triode3Filter<DataType>::Triode3Filter(Triode3Type type) : Parent(1, 1) {
  init_tube(type);
  init_state();
}

template <typename DataType>
Triode3Filter<DataType>::Triode3Filter(Triode3Type type, DataType Ri, DataType Rk, DataType Rg, DataType Ra,
                                       DataType Ci, DataType Ck, DataType Co, DataType Vk, DataType Vb, DataType Kgk)
    : Parent(1, 1) {
  init_tube(type);
  m_Ri = Ri;
  m_Rk = Rk;
  m_Rg = Rg;
  m_Ra = Ra;
  m_Ci = Ci;
  m_Ck = Ck;
  m_Co = Co;
  m_Vk = Vk;
  m_Vb = Vb;
  m_Kgk = Kgk;
  init_state();
}

template <typename DataType>
void Triode3Filter<DataType>::setup() {
  Parent::setup();
  m_T = 1. / input_sampling_rate;
  using namespace std::placeholders;
  if (m_use_newton) {
    m_Vak_func = std::bind(&Triode3Filter<DataType>::Vak_newton, this, _1, _2, _3);
  } else {
    m_Vak_func = std::bind(&Triode3Filter<DataType>::Vak_exact, this, _1, _2, _3);
  }
  m_out_func = std::bind(&Triode3Filter<DataType>::process_output_final, this, _1, _2, _3, _4, _5);
  Result r;
  process(0, m_Ri, r);
  m_Vo_n1 = 0;
  adjust_cathode_lp();
}

template <typename DataType>
void Triode3Filter<DataType>::adjust_cathode_lp() {
  if (m_T > 0. && m_auto_adjust_cathode_lp) {
    m_Rk = m_Vk / m_Ia_n1;
    m_Ck = 1 / (4 * boost::math::constants::pi<DataType>() * m_Rk);
  }
}

template <typename DataType>
DataType Triode3Filter<DataType>::Vg(DataType Vi_n, DataType Ri) const {
  DataType Vgk_n1 = m_Vgk_n1 < 0 ? 0 : m_Vgk_n1;

  // Derivatives
  DataType Vi_n_t = (Vi_n - m_Vi_n1) * input_sampling_rate;
  DataType Vk_n_t = (m_Vk_n - m_Vk_n1) * input_sampling_rate;
  DataType Vg_n_t =
      ((3 * m_Kgk * m_Rg * Ri * Vk_n_t * sqrt(Vgk_n1) + 2 * m_Rg * Vi_n_t) * m_Ci - 2 * m_Rg * m_Ig_n1 - 2 * m_Vg_n1) /
      ((2 * m_Kgk * m_Rg * Ri * sqrt(Vgk_n1) + 2 * Ri + 2 * m_Rg) * m_Ci);

  DataType ret = m_Vg_n1 + Vg_n_t * m_T;
  return ret;
}

template <typename DataType>
DataType Triode3Filter<DataType>::Vak_exact(DataType g, DataType Vgk_n, DataType Vb_n) const {
  DataType a = g / m_h3;
  DataType b = (g / m_h2) * Vgk_n * 3;
  DataType c = (g / m_h) * std::pow(Vgk_n, 2) * 3 + 1 / m_Ra;
  DataType d = g * std::pow(Vgk_n, 3) - Vb_n / m_Ra;

  DataType a_3 = std::pow(a, 3);
  DataType a_14 = std::pow(a, 14);
  DataType a_15 = a_14 * a;
  DataType b_2 = std::pow(b, 2);
  DataType b_3 = b_2 * b;

  DataType u = std::pow(27 * std::pow(a, 8) * d - 9 * std::pow(a, 7) * b * c + 2 * std::pow(a, 6) * b_3 +
                            3 * SQRT_OF_3 * sqrt(27 * std::pow(a, 16) * d * d - 18 * a_15 * b * c * d +
                                                 4 * a_15 * std::pow(c, 3) + 4 * a_14 * b_3 * d - a_14 * b_2 * c * c),
                        1. / 3);

  DataType ret =
      -(u / (3 * SQRT3_OF_2 * a_3) -
        (SQRT3_OF_2 * (2187 * std::pow(a, 5) * c - 729 * std::pow(a, 4) * b_2)) / (2187 * a_3 * u) + b / (3 * a));

  return ret;
}

template <typename DataType>
DataType Triode3Filter<DataType>::Vak_newton(DataType g, DataType Vgk_n, DataType Vb_n) const {
  DataType a = g / m_h3;
  DataType b = (g / m_h2) * Vgk_n * 3;
  DataType c = (g / m_h) * std::pow(Vgk_n, 2) * 3 + 1 / m_Ra;
  DataType d = g * std::pow(Vgk_n, 3) - Vb_n / m_Ra;

  DataType Vak_last = m_Vak_n1;
  DataType ret;
  for (int i = 0; i < mMaxIterations; ++i) {
    ret = Vak_last - (((a * Vak_last + b) * Vak_last + c) * Vak_last + d) / ((3 * a * Vak_last + 2 * b) * Vak_last + c);
    auto diff = std::abs(Vak_last - ret);
    if (diff < 1e-5) {
      break;
    }
    Vak_last = ret;
  }

  return ret;
}

template <typename DataType>
void Triode3Filter<DataType>::process(DataType Vi_n, DataType Ri, Result& r) const {
  // Grid/cathode voltage
  auto Vg_n = Vg(Vi_n, Ri);
  auto Vgk_n = Vg_n - m_Vk_n1;

  // Calculate the effective supply voltage
  auto Vb_n = m_Vb - m_Vk_n1;

  // Tube model: anode/cathode voltage
  auto g = m_gInf - (m_gInf - m_g0) * std::exp(Vgk_n / m_D);
  DataType Vak_n;
  Vak_n = m_Vak_func(g, Vgk_n, Vb_n);

  // No grid current for negative grid/vathode volatge
  DataType Ig_n = Vgk_n < 0 ? 0 : m_Kgk * std::pow(Vgk_n, 1.5);
  // No anode current for positive grid/cathode voltage
  DataType Ia_n = Vgk_n > 0 ? 0 : g * std::pow(Vak_n / m_h + Vgk_n, 3);
  DataType Ik_n = Ig_n + Ia_n;

  // Cathode low pass
  auto Vk_n = (m_T / (m_Rk * m_Ck)) * (m_Rk * Ik_n - m_Vk_n1) + m_Vk_n1;
  m_Vk_n1 = m_Vk_n;
  m_Vk_n = Vk_n;

  // Output
  auto Va_n = Vak_n + Vk_n;
  m_out_func(Va_n, g, Vgk_n, Vak_n, r);

  // Update state
  m_Vi_n1 = Vi_n;
  m_Vo_n1 = r.Vo;
  m_Vg_n1 = Vg_n;
  m_Va_n1 = Va_n;
  m_Vgk_n1 = Vgk_n;
  m_Vak_n1 = Vak_n;
  m_Ig_n1 = Ig_n;
  m_Ia_n1 = Ia_n;
}

template <typename DataType>
void Triode3Filter<DataType>::process_impl(int64_t size) const {
  assert(input_sampling_rate == output_sampling_rate);
  using namespace std::placeholders;
  Result r;
  if (get_nb_input_ports() == 2) {
    // Follower stage
    if (m_first_sample) {
      // As the input signal from a previous stage has a high amount of DC we have to set the last input sample to
      // the current one for the right steady state.
      m_Vi_n1 = converted_inputs[0][0];
      m_first_sample = false;
    }
    if (get_nb_output_ports() > 1) {
      m_out_func = std::bind(&Triode3Filter<DataType>::process_output_pre, this, _1, _2, _3, _4, _5);
      process_loop_pre_follower(size, r);
    } else {
      m_out_func = std::bind(&Triode3Filter<DataType>::process_output_final, this, _1, _2, _3, _4, _5);
      process_loop_follower(size, r);
    }
  } else {
    if (get_nb_output_ports() > 1) {
      m_out_func = std::bind(&Triode3Filter<DataType>::process_output_pre, this, _1, _2, _3, _4, _5);
      process_loop_pre(size, r);
    } else {
      m_out_func = std::bind(&Triode3Filter<DataType>::process_output_final, this, _1, _2, _3, _4, _5);
      process_loop_single(size, r);
    }
  }
}

template class Triode3Filter<float>;
template class Triode3Filter<double>;

}  // ATK
