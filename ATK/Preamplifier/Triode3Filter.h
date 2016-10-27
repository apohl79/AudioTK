/**
 * A triode preamp filter based on a modfied version of
 * http://www.hs-ulm.de/opus/frontdoor.php?source_opus=114.
 *
 * The modeled circuitry:
 *
 *                                        +Vb
 *                                         o
 *                                         |
 *                                       +---+
 *                                       |   |
 *                                       |   | Ra
 *                                       +---+
 *                                         |         ||
 *                                      Va o---------||------o-------o Vo
 *                                        _|_        ||      |
 *                                      / _|_ \      Co      |
 *         +-----+     ||           Vg |       |             |
 * Vi o----|     |-----||------o-------+-..... |             |
 *         +-----+     ||      |       | .---. | ECC82       |
 *            Ri       Ci      |        \|_ _ /            +---+
 *                             |         |                 |   |
 *                           +---+    Vk o---------+       |   | Rg
 *                           |   |       |         |       +---+
 *                           |   | Rg  +---+       |         |
 *                           +---+     |   |     -----       |
 *                             |       |   | Rk  ----- Ck    |
 *                             |       +---+       |         |
 *                             |         |         |         |
 *    o------------------------o---------o---------o---------o-------o
 *
 */
#ifndef ATK_PREAMPLIFIER_TRIODE3FILTER_H
#define ATK_PREAMPLIFIER_TRIODE3FILTER_H

#include <ATK/Core/TypedBaseFilter.h>
#include <ATK/Preamplifier/config.h>

#include <cmath>
#include <functional>

namespace ATK {

struct Triode3TypeParams {
  // Tube model
  double g0;
  double gInf;
  double D;
  double h;
  // Circuit parameters
  double Vb;  // Supply voltage
  double Ra;  // The model has the smallest error to the real tube at this value for the anode resistor.
  double Vk;  // Cathode voltage for setting the bias point.
  double Rk;  // Rk and Ck form low pass and will be set to a cutoff freq of 2Hz. Usualy you set Vk by choosing proper
  double Ck;  // values for Rk and Ck. To make things easier we turn things around an calculate Rk/Ck based on Vk.
              // See adjust_cathode_lp().
};

/// Each type needs to have a parameter set at the type's index in the TubeMap.
enum Triode3Type : uint8_t { ECC82 = 0, ECC83 };

const Triode3TypeParams TubeMap[] = {
    // ECC82
    {.g0 = 9.888576e-15,
     .gInf = 6.415385e-24,
     .D = 2.95290,
     .h = 1.021711e-02,
     .Vb = 250,
     .Ra = 10e3,
     .Vk = 6.0,
     .Rk = 861.652929,
     .Ck = 9.235444e-05},
    // ECC83
    {.g0 = 1.609683e-15,
     .gInf = 2.140844e-23,
     .D = 0.61750,
     .h = 1.000794e-02,
     .Vb = 250,
     .Ra = 50e3,
     .Vk = 1.5,
     .Rk = 1430.037044,
     .Ck = 5.564714e-05},
};

/// Triode Filter with the circuitry above. The filter has the following inputs/outputs:
///
/// Single stage mode:
///   In 0 - Vi, Out 0 - Vo
///
///   This mode is for using a single triode stage.
///
/// Pre stage mode:
///
///  In 0 - Vi, Out 0 - Va
///             Out 1 - Ri
///
///   If the stage is followed by another stage, we need to calculate Ri for the next stage, as it is dependent on the
///   input signal. We can also skip the output high pass as this is part of the RC network of the next stage, so we
///   output Va.
///
/// Follower stage mode:
///
///   In 0 - Va, Out 0 - Vo
///   In 1 - Ri
///
///   A follower stage follows another stage and requires Ri of the previous stage.
///
/// Pre/Follower stage mode:
///
///   In 0 - Va, Out 0 Va
///   In 1 - Ri, Out 1 Ri
///
///   A stage could also be connected to a pre and a follower stage.
template <typename DataType_>
class ATK_PREAMPLIFIER_EXPORT Triode3Filter : public TypedBaseFilter<DataType_> {
public:
  typedef TypedBaseFilter<DataType_> Parent;
  using typename Parent::DataType;
  using Parent::input_delay;
  using Parent::output_delay;
  using Parent::input_sampling_rate;
  using Parent::output_sampling_rate;
  using Parent::converted_inputs;
  using Parent::outputs;
  using Parent::get_nb_input_ports;
  using Parent::set_nb_input_ports;
  using Parent::get_nb_output_ports;
  using Parent::set_nb_output_ports;
  using Parent::set_input_port;

  explicit Triode3Filter(Triode3Type type = Triode3Type::ECC82);
  Triode3Filter(Triode3Type type, DataType Ri, DataType Rk, DataType Rg, DataType Ra, DataType Ci, DataType Ck,
                DataType Co, DataType Vk, DataType Vb, DataType Kgk);

  void setup() override final;

  void process_impl(int64_t size) const override final;

  /// Connect another stage
  inline void connect_stage(Triode3Filter<DataType>* stage) {
    stage->set_nb_output_ports(2);  // Pre stage mode
    set_nb_input_ports(2);          // Follower stage mode
    set_input_port(0, stage, 0);    // Va
    set_input_port(1, stage, 1);    // Ri
  }

  inline void Ri(DataType v) { m_Ri = v; }
  inline void Rk(DataType v) { m_Rk = v; }
  inline void Rg(DataType v) { m_Rg = v; }
  inline void Ra(DataType v) { m_Ra = v; }
  inline void Ci(DataType v) { m_Ci = v; }
  inline void Ck(DataType v) { m_Ck = v; }
  inline void Co(DataType v) { m_Co = v; }
  inline void Vk(DataType v) {
    m_Vk = m_Vk_n = m_Vk_n1 = v;
    adjust_cathode_lp();
  }
  inline void Vb(DataType v) { m_Vb = v; }
  inline void Kgk(DataType v) { m_Kgk = v; }

  /// Enable/disable Newton's method to calculate the tube model instead of the exact solution.
  inline void use_newton(bool b) {
    m_use_newton = b;
    using namespace std::placeholders;
    if (m_use_newton) {
      m_Vak_func = std::bind(&Triode3Filter<DataType>::Vak_newton, this, _1, _2, _3);
    } else {
      m_Vak_func = std::bind(&Triode3Filter<DataType>::Vak_exact, this, _1, _2, _3);
    }
  }

  /// Enable/disable automatic adjustment of Rk and Ck to match 2Hz cutoff based on the defined Vk.
  inline void auto_adjust_cathode_lp(bool b) { m_auto_adjust_cathode_lp = b; }

  inline DataType Ri() { return m_Ri; }
  inline DataType Rk() { return m_Rk; }
  inline DataType Rg() { return m_Rg; }
  inline DataType Ra() { return m_Ra; }
  inline DataType Ci() { return m_Ci; }
  inline DataType Ck() { return m_Ck; }
  inline DataType Co() { return m_Co; }
  inline DataType Vk() { return m_Vk; }
  inline DataType Vb() { return m_Vb; }
  inline DataType Kgk() { return m_Kgk; }
  inline bool use_newton() { return m_use_newton; }
  inline bool auto_adjust_cathode_lp() { return m_auto_adjust_cathode_lp; }

private:
  mutable bool m_first_sample = true;
  bool m_auto_adjust_cathode_lp = true;
  bool m_use_newton = true;
  constexpr static int mMaxIterations = 6;

  // Tube parameters
  DataType m_g0;
  DataType m_gInf;
  DataType m_D;
  DataType m_h;
  DataType m_h2;
  DataType m_h3;

  // Circuit parameters
  DataType m_Ri = 100e3;
  DataType m_Rg = 1000e3;
  DataType m_Ci = 10e-9;
  DataType m_Co = 10e-9;
  DataType m_Kgk = 10e-6;

  // Circuit parameters loaded from the TubeMap
  DataType m_Rk = 0;
  DataType m_Ra = 0;
  DataType m_Ck = 0;
  DataType m_Vk = 0;
  DataType m_Vb = 0;

  // Runtime state
  mutable DataType m_Vk_n = 0;
  mutable DataType m_Vk_n1 = 0;
  mutable DataType m_Vg_n1 = 0;
  mutable DataType m_Va_n1 = 0;
  mutable DataType m_Vo_n1 = 0;
  mutable DataType m_Vi_n1 = 0;
  mutable DataType m_Vgk_n1 = 0;
  mutable DataType m_Vak_n1 = 0;
  mutable DataType m_Ig_n1 = 0;
  mutable DataType m_Ia_n1 = 0;

  struct Result {
    DataType Vo;
    DataType Va;
    DataType Ri;
  };

  /// Helpers
  DataType SQRT_OF_3 = std::sqrt(3);
  DataType SQRT3_OF_2 = std::pow(2, 1. / 3);
  DataType m_T = 0;

  /// Load tube parameters from the TubeMap
  void init_tube(Triode3Type type) {
    auto& tube = TubeMap[type];
    m_g0 = tube.g0;
    m_gInf = tube.gInf;
    m_D = tube.D;
    m_h = tube.h;
    m_h2 = std::pow(m_h, 2);
    m_h3 = std::pow(m_h, 3);
    m_Vb = tube.Vb;
    m_Ra = tube.Ra;
    m_Vk = tube.Vk;
    m_Rk = tube.Rk;
    m_Ck = tube.Ck;
    input_delay = output_delay = 1;
  }

  /// Initialize the runtime state
  void init_state() {
    m_Vg_n1 = m_Vo_n1 = m_Vi_n1 = m_Ig_n1 = m_Ia_n1 = 0;
    m_Va_n1 = m_Vak_n1 = m_Vb;
    m_Vk_n = m_Vk_n1 = m_Vk;
    m_Vgk_n1 = m_Vg_n1 - m_Vk_n1;
  }

  /// Adjustment of Rk and Ck to match 2Hz cutoff
  void adjust_cathode_lp();

  /// Calculate grid voltage
  DataType_ Vg(DataType Vi_n, DataType Ri) const;

  /// Tube model function
  using VakFuncType = std::function<DataType(DataType, DataType, DataType)>;
  VakFuncType m_Vak_func;

  /// Caclulate anode/kathode voltage
  DataType_ Vak_exact(DataType g, DataType Vgk_n, DataType Vb_n) const;

  /// Calculate anode/kathode voltage via Newton's method
  DataType_ Vak_newton(DataType g, DataType Vgk_n, DataType Vb_n) const;

  /// Output function (depends on the connected outputs)
  using OutFuncType = std::function<void(DataType, DataType, DataType, DataType, Result&)>;
  mutable OutFuncType m_out_func;

  /// Final output with no follower stage.
  void process_output_final(DataType Va_n, DataType, DataType, DataType, Result& r) const {
    // Output DC blocker
    r.Vo = Va_n - m_Va_n1 + m_Vo_n1 * (1 - m_T / (m_Co * m_Rg));
    r.Va = 0;
    r.Ri = 0;
  }

  /// Calculate output for a follower stage.
  void process_output_pre(DataType Va_n, DataType g, DataType Vgk_n, DataType Vak_n, Result& r) const {
    // Pre stage mode, calculate internal resistance and pass anode voltage to next stage
    DataType Xi = 1 / (3 * g / m_h * std::pow(Vgk_n + Vak_n / m_h, 2));
    r.Ri = m_Ra * Xi / (m_Ra + Xi);
    r.Va = Va_n;
    r.Vo = 0;
  }

  /// Process next sample
  void process(DataType Vi_n, DataType Ri, Result& r) const;

  /// Process loop for single stage mode
  inline void process_loop_single(int64_t size, Result& r) const {
    for (int64_t s = 0; s < size; ++s) {
      DataType Vi = converted_inputs[0][s];
      process(Vi, m_Ri, r);
      outputs[0][s] = r.Vo;
    }
  }

  /// Process loop for pre stage, no follower
  inline void process_loop_pre(int64_t size, Result& r) const {
    for (int64_t s = 0; s < size; ++s) {
      DataType Vi = converted_inputs[0][s];
      process(Vi, m_Ri, r);
      outputs[0][s] = r.Va;
      outputs[1][s] = r.Ri;
    }
  }

  /// Process loop for pre stage and follower
  inline void process_loop_pre_follower(int64_t size, Result& r) const {
    for (int64_t s = 0; s < size; ++s) {
      DataType Vi = converted_inputs[0][s];
      DataType Ri = converted_inputs[1][s];
      process(Vi, Ri, r);
      outputs[0][s] = r.Va;
      outputs[1][s] = r.Ri;
    }
  }

  /// Process loop for follower stage, no pre
  inline void process_loop_follower(int64_t size, Result& r) const {
    for (int64_t s = 0; s < size; ++s) {
      DataType Vi = converted_inputs[0][s];
      DataType Ri = converted_inputs[1][s];
      process(Vi, Ri, r);
      outputs[0][s] = r.Vo;
    }
  }
};

}  // ATK

#endif  // ATK_PREAMPLIFIER_TRIODE3FILTER_H
