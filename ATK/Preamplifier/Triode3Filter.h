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
 *         +-----+     ||      |       | .---. | Tube        |
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

using namespace std::placeholders;

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

/// Each type has a parameter set at the type's index in the TubeMap.
enum Triode3Type : uint8_t { ECC82 = 0, ECC83 };

/// Tube models
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
                DataType Co, DataType Vk, DataType Vb);

  void setup() override final;

  void process_impl(int64_t size) const override final;

  /// Connect another stage
  inline void connect_stage(Triode3Filter<DataType>* stage) {
    stage->pre_stage(true);       // Pre stage mode
    follower_stage(true);         // Follower stage mode
    set_input_port(0, stage, 0);  // Va
    set_input_port(1, stage, 1);  // Ri
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

  /// Enable/disable automatic adjustment of Rk and Ck to match 2Hz cutoff based on the defined Vk.
  inline void auto_adjust_cathode_lp(bool b) { m_auto_adjust_cathode_lp = b; }

  /// Enable/disable the pre stage mode.
  inline void pre_stage(bool b) {
    m_pre_stage = b;
    int ports = 1;
    if (b) {
      ports = 2;
    }
    set_nb_output_ports(ports);
    setup_process_functions();
  }

  /// Enable/disable the follower stage mode.
  inline void follower_stage(bool b) {
    m_follower_stage = b;
    int ports = 1;
    if (b) {
      ports = 2;
    }
    set_nb_input_ports(ports);
    setup_process_functions();
  }

  inline DataType Ri() const { return m_Ri; }
  inline DataType Rk() const { return m_Rk; }
  inline DataType Rg() const { return m_Rg; }
  inline DataType Ra() const { return m_Ra; }
  inline DataType Ci() const { return m_Ci; }
  inline DataType Ck() const { return m_Ck; }
  inline DataType Co() const { return m_Co; }
  inline DataType Vk() const { return m_Vk; }
  inline DataType Vb() const { return m_Vb; }
  inline bool auto_adjust_cathode_lp() const { return m_auto_adjust_cathode_lp; }
  inline bool pre_stage() const { return m_pre_stage; }
  inline bool follower_stage() const { return m_follower_stage; }

private:
  mutable bool m_first_sample = true;
  bool m_auto_adjust_cathode_lp = true;
  bool m_pre_stage = false;
  bool m_follower_stage = false;
  constexpr static int m_max_iterations = 6;

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
  DataType m_Vgamma = 0.6;

  // Circuit parameters loaded from the TubeMap
  DataType m_Rk = 0;
  DataType m_Ra = 0;
  DataType m_Ck = 0;
  DataType m_Vk = 0;
  DataType m_Vb = 0;

  // Runtime state
  DataType m_Vk_n = 0;
  DataType m_Vk_n1 = 0;
  DataType m_Vg_n1 = 0;
  DataType m_Va_n1 = 0;
  DataType m_Vo_n1 = 0;
  DataType m_Vi_n1 = 0;
  DataType m_Vgk_n1 = 0;
  DataType m_Vak_n1 = 0;
  DataType m_Ig_n1 = 0;
  DataType m_Ia_n1 = 0;

  struct Result {
    DataType Vo;
    DataType Va;
    DataType Ri;
  };

  /// Time delta between two samples
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
  }

  /// Initialize the runtime state
  void init_state() {
    m_Vg_n1 = m_Vo_n1 = m_Vi_n1 = m_Ig_n1 = m_Ia_n1 = 0;
    m_Va_n1 = m_Vak_n1 = m_Vb;
    m_Vk_n = m_Vk_n1 = m_Vk;
    m_Vgk_n1 = m_Vg_n1 - m_Vk_n1;
    input_delay = output_delay = 1;
  }

  /// Adjustment of Rk and Ck to match 2Hz cutoff
  void adjust_cathode_lp();

  /// Calculate grid voltage
  DataType_ Vg(DataType Vi_n, DataType Ri) const;

  /// Tube model: Calculate anode/cathode voltage via Newton's method
  DataType_ Vak(DataType g, DataType Vgk_n, DataType Vb_n) const;

  /// Output function (depends on the connected outputs)
  using OutFuncType = std::function<void(DataType, DataType, DataType, DataType, Result&)>;
  OutFuncType m_out_func;

  /// Final output with no follower stage.
  void process_output_final(DataType Va_n, DataType, DataType, DataType, Result& r) {
    // Output DC blocker
    r.Vo = Va_n - m_Va_n1 + m_Vo_n1 * (1 - m_T / (m_Co * m_Rg));
    r.Va = 0;
    r.Ri = 0;
  }

  /// Calculate output for a follower stage.
  void process_output_pre(DataType Va_n, DataType g, DataType Vgk_n, DataType Vak_n, Result& r) {
    // Pre stage mode, calculate internal resistance and pass anode voltage to next stage
    DataType Xi = 1 / (3 * g / m_h * std::pow(Vgk_n + Vak_n / m_h, 2));
    r.Ri = m_Ra * Xi / (m_Ra + Xi);
    r.Va = Va_n;
    r.Vo = 0;
  }

  /// Process next sample
  void process(DataType Vi_n, DataType Ri, Result& r);

  /// Process loop function for the defined mode (pre/follower)
  using LoopFuncType = std::function<void(int64_t, Result&)>;
  LoopFuncType m_loop_func;

  /// Process loop for single stage mode
  void process_loop_single(int64_t size, Result& r) {
    for (int64_t s = 0; s < size; ++s) {
      DataType Vi = converted_inputs[0][s];
      process(Vi, m_Ri, r);
      outputs[0][s] = r.Vo;
    }
  }

  /// Process loop for pre stage, no follower
  void process_loop_pre(int64_t size, Result& r) {
    for (int64_t s = 0; s < size; ++s) {
      DataType Vi = converted_inputs[0][s];
      process(Vi, m_Ri, r);
      outputs[0][s] = r.Va;
      outputs[1][s] = r.Ri;
    }
  }

  /// As the input signal from a previous stage has a high amount of DC we have to set the last input sample to
  /// the current one for the right steady state.
  inline void prepare_follower() {
    if (m_first_sample) {
      m_Vi_n1 = converted_inputs[0][0];
      m_first_sample = false;
    }
  }

  /// Process loop for pre stage and follower
  void process_loop_pre_follower(int64_t size, Result& r) {
    prepare_follower();
    for (int64_t s = 0; s < size; ++s) {
      DataType Vi = converted_inputs[0][s];
      DataType Ri = converted_inputs[1][s];
      process(Vi, Ri, r);
      outputs[0][s] = r.Va;
      outputs[1][s] = r.Ri;
    }
  }

  /// Process loop for follower stage, no pre
  void process_loop_follower(int64_t size, Result& r) {
    prepare_follower();
    for (int64_t s = 0; s < size; ++s) {
      DataType Vi = converted_inputs[0][s];
      DataType Ri = converted_inputs[1][s];
      process(Vi, Ri, r);
      outputs[0][s] = r.Vo;
    }
  }

  /// Setup the loop function pointer based on the pre/follower flags
  inline void setup_process_functions() {
    if (follower_stage()) {
      if (pre_stage()) {
        m_out_func = std::bind(&Triode3Filter<DataType>::process_output_pre, this, _1, _2, _3, _4, _5);
        m_loop_func = std::bind(&Triode3Filter<DataType>::process_loop_pre_follower, this, _1, _2);
      } else {
        m_out_func = std::bind(&Triode3Filter<DataType>::process_output_final, this, _1, _2, _3, _4, _5);
        m_loop_func = std::bind(&Triode3Filter<DataType>::process_loop_follower, this, _1, _2);
      }
    } else {
      if (pre_stage()) {
        m_out_func = std::bind(&Triode3Filter<DataType>::process_output_pre, this, _1, _2, _3, _4, _5);
        m_loop_func = std::bind(&Triode3Filter<DataType>::process_loop_pre, this, _1, _2);
      } else {
        m_out_func = std::bind(&Triode3Filter<DataType>::process_output_final, this, _1, _2, _3, _4, _5);
        m_loop_func = std::bind(&Triode3Filter<DataType>::process_loop_single, this, _1, _2);
      }
    }
  }
};

}  // ATK

#endif  // ATK_PREAMPLIFIER_TRIODE3FILTER_H
