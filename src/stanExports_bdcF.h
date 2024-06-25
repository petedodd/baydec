// Generated by rstantools.  Do not edit by hand.

/*
    baydec is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    baydec is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with baydec.  If not, see <http://www.gnu.org/licenses/>.
*/
#ifndef MODELS_HPP
#define MODELS_HPP
#define STAN__SERVICES__COMMAND_HPP
#ifndef USE_STANC3
#define USE_STANC3
#endif
#include <rstan/rstaninc.hpp>
// Code generated by stanc v2.32.2
#include <stan/model/model_header.hpp>
namespace model_bdcF_namespace {
using stan::model::model_base_crtp;
using namespace stan::math;
stan::math::profile_map profiles__;
static constexpr std::array<const char*, 31> locations_array__ =
  {" (found before start of program)",
  " (in 'bdcF', line 15, column 2 to column 20)",
  " (in 'bdcF', line 18, column 2 to column 29)",
  " (in 'bdcF', line 39, column 2 to column 18)",
  " (in 'bdcF', line 40, column 17 to column 44)",
  " (in 'bdcF', line 40, column 2 to column 44)",
  " (in 'bdcF', line 22, column 2 to column 31)",
  " (in 'bdcF', line 31, column 2 to column 39)",
  " (in 'bdcF', line 32, column 2 to column 53)",
  " (in 'bdcF', line 2, column 2 to column 9)",
  " (in 'bdcF', line 3, column 2 to column 9)",
  " (in 'bdcF', line 4, column 2 to column 9)",
  " (in 'bdcF', line 5, column 9 to column 11)",
  " (in 'bdcF', line 5, column 12 to column 14)",
  " (in 'bdcF', line 5, column 2 to column 18)",
  " (in 'bdcF', line 6, column 9 to column 11)",
  " (in 'bdcF', line 6, column 12 to column 14)",
  " (in 'bdcF', line 6, column 2 to column 18)",
  " (in 'bdcF', line 7, column 2 to column 11)",
  " (in 'bdcF', line 11, column 9 to column 11)",
  " (in 'bdcF', line 11, column 12 to column 14)",
  " (in 'bdcF', line 11, column 2 to column 29)",
  " (in 'bdcF', line 12, column 9 to column 11)",
  " (in 'bdcF', line 12, column 12 to column 14)",
  " (in 'bdcF', line 12, column 2 to column 29)",
  " (in 'bdcF', line 15, column 9 to column 11)",
  " (in 'bdcF', line 15, column 12 to column 14)",
  " (in 'bdcF', line 18, column 9 to column 11)",
  " (in 'bdcF', line 18, column 12 to column 14)",
  " (in 'bdcF', line 39, column 9 to column 11)",
  " (in 'bdcF', line 39, column 12 to column 14)"};
#include <stan_meta_header.hpp>
class model_bdcF final : public model_base_crtp<model_bdcF> {
private:
  int NG;
  int NC;
  int NR;
  Eigen::Matrix<double,-1,-1> Y_data__;
  Eigen::Matrix<double,-1,-1> X_data__;
  double tol;
  Eigen::Matrix<double,-1,-1> YtX_data__;
  Eigen::Matrix<double,-1,-1> XtX_data__;
  Eigen::Map<Eigen::Matrix<double,-1,-1>> Y{nullptr, 0, 0};
  Eigen::Map<Eigen::Matrix<double,-1,-1>> X{nullptr, 0, 0};
  Eigen::Map<Eigen::Matrix<double,-1,-1>> YtX{nullptr, 0, 0};
  Eigen::Map<Eigen::Matrix<double,-1,-1>> XtX{nullptr, 0, 0};
public:
  ~model_bdcF() {}
  model_bdcF(stan::io::var_context& context__, unsigned int
             random_seed__ = 0, std::ostream* pstream__ = nullptr)
      : model_base_crtp(0) {
    int current_statement__ = 0;
    using local_scalar_t__ = double;
    boost::ecuyer1988 base_rng__ =
      stan::services::util::create_rng(random_seed__, 0);
    // suppress unused var warning
    (void) base_rng__;
    static constexpr const char* function__ =
      "model_bdcF_namespace::model_bdcF";
    // suppress unused var warning
    (void) function__;
    local_scalar_t__ DUMMY_VAR__(std::numeric_limits<double>::quiet_NaN());
    // suppress unused var warning
    (void) DUMMY_VAR__;
    try {
      int pos__ = std::numeric_limits<int>::min();
      pos__ = 1;
      current_statement__ = 9;
      context__.validate_dims("data initialization", "NG", "int",
        std::vector<size_t>{});
      NG = std::numeric_limits<int>::min();
      current_statement__ = 9;
      NG = context__.vals_i("NG")[(1 - 1)];
      current_statement__ = 10;
      context__.validate_dims("data initialization", "NC", "int",
        std::vector<size_t>{});
      NC = std::numeric_limits<int>::min();
      current_statement__ = 10;
      NC = context__.vals_i("NC")[(1 - 1)];
      current_statement__ = 11;
      context__.validate_dims("data initialization", "NR", "int",
        std::vector<size_t>{});
      NR = std::numeric_limits<int>::min();
      current_statement__ = 11;
      NR = context__.vals_i("NR")[(1 - 1)];
      current_statement__ = 12;
      stan::math::validate_non_negative_index("Y", "NG", NG);
      current_statement__ = 13;
      stan::math::validate_non_negative_index("Y", "NR", NR);
      current_statement__ = 14;
      context__.validate_dims("data initialization", "Y", "double",
        std::vector<size_t>{static_cast<size_t>(NG), static_cast<size_t>(NR)});
      Y_data__ = Eigen::Matrix<double,-1,-1>::Constant(NG, NR,
                   std::numeric_limits<double>::quiet_NaN());
      new (&Y) Eigen::Map<Eigen::Matrix<double,-1,-1>>(Y_data__.data(), NG,
        NR);
      {
        std::vector<local_scalar_t__> Y_flat__;
        current_statement__ = 14;
        Y_flat__ = context__.vals_r("Y");
        current_statement__ = 14;
        pos__ = 1;
        current_statement__ = 14;
        for (int sym1__ = 1; sym1__ <= NR; ++sym1__) {
          current_statement__ = 14;
          for (int sym2__ = 1; sym2__ <= NG; ++sym2__) {
            current_statement__ = 14;
            stan::model::assign(Y, Y_flat__[(pos__ - 1)],
              "assigning variable Y", stan::model::index_uni(sym2__),
              stan::model::index_uni(sym1__));
            current_statement__ = 14;
            pos__ = (pos__ + 1);
          }
        }
      }
      current_statement__ = 15;
      stan::math::validate_non_negative_index("X", "NG", NG);
      current_statement__ = 16;
      stan::math::validate_non_negative_index("X", "NC", NC);
      current_statement__ = 17;
      context__.validate_dims("data initialization", "X", "double",
        std::vector<size_t>{static_cast<size_t>(NG), static_cast<size_t>(NC)});
      X_data__ = Eigen::Matrix<double,-1,-1>::Constant(NG, NC,
                   std::numeric_limits<double>::quiet_NaN());
      new (&X) Eigen::Map<Eigen::Matrix<double,-1,-1>>(X_data__.data(), NG,
        NC);
      {
        std::vector<local_scalar_t__> X_flat__;
        current_statement__ = 17;
        X_flat__ = context__.vals_r("X");
        current_statement__ = 17;
        pos__ = 1;
        current_statement__ = 17;
        for (int sym1__ = 1; sym1__ <= NC; ++sym1__) {
          current_statement__ = 17;
          for (int sym2__ = 1; sym2__ <= NG; ++sym2__) {
            current_statement__ = 17;
            stan::model::assign(X, X_flat__[(pos__ - 1)],
              "assigning variable X", stan::model::index_uni(sym2__),
              stan::model::index_uni(sym1__));
            current_statement__ = 17;
            pos__ = (pos__ + 1);
          }
        }
      }
      current_statement__ = 18;
      context__.validate_dims("data initialization", "tol", "double",
        std::vector<size_t>{});
      tol = std::numeric_limits<double>::quiet_NaN();
      current_statement__ = 18;
      tol = context__.vals_r("tol")[(1 - 1)];
      current_statement__ = 19;
      stan::math::validate_non_negative_index("YtX", "NR", NR);
      current_statement__ = 20;
      stan::math::validate_non_negative_index("YtX", "NC", NC);
      current_statement__ = 21;
      YtX_data__ = Eigen::Matrix<double,-1,-1>::Constant(NR, NC,
                     std::numeric_limits<double>::quiet_NaN());
      new (&YtX) Eigen::Map<Eigen::Matrix<double,-1,-1>>(YtX_data__.data(),
        NR, NC);
      current_statement__ = 21;
      stan::model::assign(YtX,
        stan::math::multiply(stan::math::transpose(Y), X),
        "assigning variable YtX");
      current_statement__ = 22;
      stan::math::validate_non_negative_index("XtX", "NC", NC);
      current_statement__ = 23;
      stan::math::validate_non_negative_index("XtX", "NC", NC);
      current_statement__ = 24;
      XtX_data__ = Eigen::Matrix<double,-1,-1>::Constant(NC, NC,
                     std::numeric_limits<double>::quiet_NaN());
      new (&XtX) Eigen::Map<Eigen::Matrix<double,-1,-1>>(XtX_data__.data(),
        NC, NC);
      current_statement__ = 24;
      stan::model::assign(XtX,
        stan::math::multiply(stan::math::transpose(X), X),
        "assigning variable XtX");
      current_statement__ = 25;
      stan::math::validate_non_negative_index("bet", "NC", NC);
      current_statement__ = 26;
      stan::math::validate_non_negative_index("bet", "NR", NR);
      current_statement__ = 27;
      stan::math::validate_non_negative_index("W", "NC", NC);
      current_statement__ = 28;
      stan::math::validate_non_negative_index("W", "NR", NR);
      current_statement__ = 29;
      stan::math::validate_non_negative_index("P", "NC", NC);
      current_statement__ = 30;
      stan::math::validate_non_negative_index("P", "NR", NR);
    } catch (const std::exception& e) {
      stan::lang::rethrow_located(e, locations_array__[current_statement__]);
    }
    num_params_r__ = (NC * NR);
  }
  inline std::string model_name() const final {
    return "model_bdcF";
  }
  inline std::vector<std::string> model_compile_info() const noexcept {
    return std::vector<std::string>{"stanc_version = stanc3 v2.32.2",
             "stancflags = --allow-undefined"};
  }
  template <bool propto__, bool jacobian__, typename VecR, typename VecI,
            stan::require_vector_like_t<VecR>* = nullptr,
            stan::require_vector_like_vt<std::is_integral, VecI>* = nullptr>
  inline stan::scalar_type_t<VecR>
  log_prob_impl(VecR& params_r__, VecI& params_i__, std::ostream*
                pstream__ = nullptr) const {
    using T__ = stan::scalar_type_t<VecR>;
    using local_scalar_t__ = T__;
    T__ lp__(0.0);
    stan::math::accumulator<T__> lp_accum__;
    stan::io::deserializer<local_scalar_t__> in__(params_r__, params_i__);
    int current_statement__ = 0;
    local_scalar_t__ DUMMY_VAR__(std::numeric_limits<double>::quiet_NaN());
    // suppress unused var warning
    (void) DUMMY_VAR__;
    static constexpr const char* function__ =
      "model_bdcF_namespace::log_prob";
    // suppress unused var warning
    (void) function__;
    try {
      Eigen::Matrix<local_scalar_t__,-1,-1> bet =
        Eigen::Matrix<local_scalar_t__,-1,-1>::Constant(NC, NR, DUMMY_VAR__);
      current_statement__ = 1;
      bet = in__.template read<Eigen::Matrix<local_scalar_t__,-1,-1>>(NC, NR);
      Eigen::Matrix<local_scalar_t__,-1,-1> W =
        Eigen::Matrix<local_scalar_t__,-1,-1>::Constant(NC, NR, DUMMY_VAR__);
      current_statement__ = 2;
      stan::model::assign(W, stan::math::exp(bet), "assigning variable W");
      {
        current_statement__ = 6;
        lp_accum__.add(stan::math::normal_lpdf<propto__>(
                         stan::math::to_vector(bet), 0, 5));
        current_statement__ = 7;
        lp_accum__.add((stan::math::trace(stan::math::multiply(YtX, W)) /
          (tol * tol)));
        current_statement__ = 8;
        lp_accum__.add((-stan::math::trace_quad_form(XtX, W) / ((2 * tol) *
          tol)));
      }
    } catch (const std::exception& e) {
      stan::lang::rethrow_located(e, locations_array__[current_statement__]);
    }
    lp_accum__.add(lp__);
    return lp_accum__.sum();
  }
  template <typename RNG, typename VecR, typename VecI, typename VecVar,
            stan::require_vector_like_vt<std::is_floating_point,
            VecR>* = nullptr, stan::require_vector_like_vt<std::is_integral,
            VecI>* = nullptr, stan::require_vector_vt<std::is_floating_point,
            VecVar>* = nullptr>
  inline void
  write_array_impl(RNG& base_rng__, VecR& params_r__, VecI& params_i__,
                   VecVar& vars__, const bool
                   emit_transformed_parameters__ = true, const bool
                   emit_generated_quantities__ = true, std::ostream*
                   pstream__ = nullptr) const {
    using local_scalar_t__ = double;
    stan::io::deserializer<local_scalar_t__> in__(params_r__, params_i__);
    stan::io::serializer<local_scalar_t__> out__(vars__);
    static constexpr bool propto__ = true;
    // suppress unused var warning
    (void) propto__;
    double lp__ = 0.0;
    // suppress unused var warning
    (void) lp__;
    int current_statement__ = 0;
    stan::math::accumulator<double> lp_accum__;
    local_scalar_t__ DUMMY_VAR__(std::numeric_limits<double>::quiet_NaN());
    // suppress unused var warning
    (void) DUMMY_VAR__;
    constexpr bool jacobian__ = false;
    static constexpr const char* function__ =
      "model_bdcF_namespace::write_array";
    // suppress unused var warning
    (void) function__;
    try {
      Eigen::Matrix<double,-1,-1> bet =
        Eigen::Matrix<double,-1,-1>::Constant(NC, NR,
          std::numeric_limits<double>::quiet_NaN());
      current_statement__ = 1;
      bet = in__.template read<Eigen::Matrix<local_scalar_t__,-1,-1>>(NC, NR);
      Eigen::Matrix<double,-1,-1> W =
        Eigen::Matrix<double,-1,-1>::Constant(NC, NR,
          std::numeric_limits<double>::quiet_NaN());
      out__.write(bet);
      if (stan::math::logical_negation(
            (stan::math::primitive_value(emit_transformed_parameters__) ||
            stan::math::primitive_value(emit_generated_quantities__)))) {
        return ;
      }
      current_statement__ = 2;
      stan::model::assign(W, stan::math::exp(bet), "assigning variable W");
      if (emit_transformed_parameters__) {
        out__.write(W);
      }
      if (stan::math::logical_negation(emit_generated_quantities__)) {
        return ;
      }
      Eigen::Matrix<double,-1,-1> P =
        Eigen::Matrix<double,-1,-1>::Constant(NC, NR,
          std::numeric_limits<double>::quiet_NaN());
      current_statement__ = 5;
      for (int j = 1; j <= NR; ++j) {
        current_statement__ = 4;
        stan::model::assign(P,
          stan::math::softmax(
            stan::model::rvalue(bet, "bet", stan::model::index_omni(),
              stan::model::index_uni(j))), "assigning variable P",
          stan::model::index_omni(), stan::model::index_uni(j));
      }
      out__.write(P);
    } catch (const std::exception& e) {
      stan::lang::rethrow_located(e, locations_array__[current_statement__]);
    }
  }
  template <typename VecVar, typename VecI,
            stan::require_vector_t<VecVar>* = nullptr,
            stan::require_vector_like_vt<std::is_integral, VecI>* = nullptr>
  inline void
  unconstrain_array_impl(const VecVar& params_r__, const VecI& params_i__,
                         VecVar& vars__, std::ostream* pstream__ = nullptr) const {
    using local_scalar_t__ = double;
    stan::io::deserializer<local_scalar_t__> in__(params_r__, params_i__);
    stan::io::serializer<local_scalar_t__> out__(vars__);
    int current_statement__ = 0;
    local_scalar_t__ DUMMY_VAR__(std::numeric_limits<double>::quiet_NaN());
    // suppress unused var warning
    (void) DUMMY_VAR__;
    try {
      int pos__ = std::numeric_limits<int>::min();
      pos__ = 1;
      Eigen::Matrix<local_scalar_t__,-1,-1> bet =
        Eigen::Matrix<local_scalar_t__,-1,-1>::Constant(NC, NR, DUMMY_VAR__);
      current_statement__ = 1;
      stan::model::assign(bet,
        in__.read<Eigen::Matrix<local_scalar_t__,-1,-1>>(NC, NR),
        "assigning variable bet");
      out__.write(bet);
    } catch (const std::exception& e) {
      stan::lang::rethrow_located(e, locations_array__[current_statement__]);
    }
  }
  template <typename VecVar, stan::require_vector_t<VecVar>* = nullptr>
  inline void
  transform_inits_impl(const stan::io::var_context& context__, VecVar&
                       vars__, std::ostream* pstream__ = nullptr) const {
    using local_scalar_t__ = double;
    stan::io::serializer<local_scalar_t__> out__(vars__);
    int current_statement__ = 0;
    local_scalar_t__ DUMMY_VAR__(std::numeric_limits<double>::quiet_NaN());
    // suppress unused var warning
    (void) DUMMY_VAR__;
    try {
      current_statement__ = 1;
      context__.validate_dims("parameter initialization", "bet", "double",
        std::vector<size_t>{static_cast<size_t>(NC), static_cast<size_t>(NR)});
      int pos__ = std::numeric_limits<int>::min();
      pos__ = 1;
      Eigen::Matrix<local_scalar_t__,-1,-1> bet =
        Eigen::Matrix<local_scalar_t__,-1,-1>::Constant(NC, NR, DUMMY_VAR__);
      {
        std::vector<local_scalar_t__> bet_flat__;
        current_statement__ = 1;
        bet_flat__ = context__.vals_r("bet");
        current_statement__ = 1;
        pos__ = 1;
        current_statement__ = 1;
        for (int sym1__ = 1; sym1__ <= NR; ++sym1__) {
          current_statement__ = 1;
          for (int sym2__ = 1; sym2__ <= NC; ++sym2__) {
            current_statement__ = 1;
            stan::model::assign(bet, bet_flat__[(pos__ - 1)],
              "assigning variable bet", stan::model::index_uni(sym2__),
              stan::model::index_uni(sym1__));
            current_statement__ = 1;
            pos__ = (pos__ + 1);
          }
        }
      }
      out__.write(bet);
    } catch (const std::exception& e) {
      stan::lang::rethrow_located(e, locations_array__[current_statement__]);
    }
  }
  inline void
  get_param_names(std::vector<std::string>& names__, const bool
                  emit_transformed_parameters__ = true, const bool
                  emit_generated_quantities__ = true) const {
    names__ = std::vector<std::string>{"bet"};
    if (emit_transformed_parameters__) {
      std::vector<std::string> temp{"W"};
      names__.reserve(names__.size() + temp.size());
      names__.insert(names__.end(), temp.begin(), temp.end());
    }
    if (emit_generated_quantities__) {
      std::vector<std::string> temp{"P"};
      names__.reserve(names__.size() + temp.size());
      names__.insert(names__.end(), temp.begin(), temp.end());
    }
  }
  inline void
  get_dims(std::vector<std::vector<size_t>>& dimss__, const bool
           emit_transformed_parameters__ = true, const bool
           emit_generated_quantities__ = true) const {
    dimss__ = std::vector<std::vector<size_t>>{std::vector<size_t>{static_cast<
                                                                    size_t>(
                                                                    NC),
                                                 static_cast<size_t>(NR)}};
    if (emit_transformed_parameters__) {
      std::vector<std::vector<size_t>>
        temp{std::vector<size_t>{static_cast<size_t>(NC),
               static_cast<size_t>(NR)}};
      dimss__.reserve(dimss__.size() + temp.size());
      dimss__.insert(dimss__.end(), temp.begin(), temp.end());
    }
    if (emit_generated_quantities__) {
      std::vector<std::vector<size_t>>
        temp{std::vector<size_t>{static_cast<size_t>(NC),
               static_cast<size_t>(NR)}};
      dimss__.reserve(dimss__.size() + temp.size());
      dimss__.insert(dimss__.end(), temp.begin(), temp.end());
    }
  }
  inline void
  constrained_param_names(std::vector<std::string>& param_names__, bool
                          emit_transformed_parameters__ = true, bool
                          emit_generated_quantities__ = true) const final {
    for (int sym1__ = 1; sym1__ <= NR; ++sym1__) {
      for (int sym2__ = 1; sym2__ <= NC; ++sym2__) {
        param_names__.emplace_back(std::string() + "bet" + '.' +
          std::to_string(sym2__) + '.' + std::to_string(sym1__));
      }
    }
    if (emit_transformed_parameters__) {
      for (int sym1__ = 1; sym1__ <= NR; ++sym1__) {
        for (int sym2__ = 1; sym2__ <= NC; ++sym2__) {
          param_names__.emplace_back(std::string() + "W" + '.' +
            std::to_string(sym2__) + '.' + std::to_string(sym1__));
        }
      }
    }
    if (emit_generated_quantities__) {
      for (int sym1__ = 1; sym1__ <= NR; ++sym1__) {
        for (int sym2__ = 1; sym2__ <= NC; ++sym2__) {
          param_names__.emplace_back(std::string() + "P" + '.' +
            std::to_string(sym2__) + '.' + std::to_string(sym1__));
        }
      }
    }
  }
  inline void
  unconstrained_param_names(std::vector<std::string>& param_names__, bool
                            emit_transformed_parameters__ = true, bool
                            emit_generated_quantities__ = true) const final {
    for (int sym1__ = 1; sym1__ <= NR; ++sym1__) {
      for (int sym2__ = 1; sym2__ <= NC; ++sym2__) {
        param_names__.emplace_back(std::string() + "bet" + '.' +
          std::to_string(sym2__) + '.' + std::to_string(sym1__));
      }
    }
    if (emit_transformed_parameters__) {
      for (int sym1__ = 1; sym1__ <= NR; ++sym1__) {
        for (int sym2__ = 1; sym2__ <= NC; ++sym2__) {
          param_names__.emplace_back(std::string() + "W" + '.' +
            std::to_string(sym2__) + '.' + std::to_string(sym1__));
        }
      }
    }
    if (emit_generated_quantities__) {
      for (int sym1__ = 1; sym1__ <= NR; ++sym1__) {
        for (int sym2__ = 1; sym2__ <= NC; ++sym2__) {
          param_names__.emplace_back(std::string() + "P" + '.' +
            std::to_string(sym2__) + '.' + std::to_string(sym1__));
        }
      }
    }
  }
  inline std::string get_constrained_sizedtypes() const {
    return std::string("[{\"name\":\"bet\",\"type\":{\"name\":\"matrix\",\"rows\":" + std::to_string(NC) + ",\"cols\":" + std::to_string(NR) + "},\"block\":\"parameters\"},{\"name\":\"W\",\"type\":{\"name\":\"matrix\",\"rows\":" + std::to_string(NC) + ",\"cols\":" + std::to_string(NR) + "},\"block\":\"transformed_parameters\"},{\"name\":\"P\",\"type\":{\"name\":\"matrix\",\"rows\":" + std::to_string(NC) + ",\"cols\":" + std::to_string(NR) + "},\"block\":\"generated_quantities\"}]");
  }
  inline std::string get_unconstrained_sizedtypes() const {
    return std::string("[{\"name\":\"bet\",\"type\":{\"name\":\"matrix\",\"rows\":" + std::to_string(NC) + ",\"cols\":" + std::to_string(NR) + "},\"block\":\"parameters\"},{\"name\":\"W\",\"type\":{\"name\":\"matrix\",\"rows\":" + std::to_string(NC) + ",\"cols\":" + std::to_string(NR) + "},\"block\":\"transformed_parameters\"},{\"name\":\"P\",\"type\":{\"name\":\"matrix\",\"rows\":" + std::to_string(NC) + ",\"cols\":" + std::to_string(NR) + "},\"block\":\"generated_quantities\"}]");
  }
  // Begin method overload boilerplate
  template <typename RNG> inline void
  write_array(RNG& base_rng, Eigen::Matrix<double,-1,1>& params_r,
              Eigen::Matrix<double,-1,1>& vars, const bool
              emit_transformed_parameters = true, const bool
              emit_generated_quantities = true, std::ostream*
              pstream = nullptr) const {
    const size_t num_params__ = (NC * NR);
    const size_t num_transformed = emit_transformed_parameters * ((NC * NR));
    const size_t num_gen_quantities = emit_generated_quantities * ((NC * NR));
    const size_t num_to_write = num_params__ + num_transformed +
      num_gen_quantities;
    std::vector<int> params_i;
    vars = Eigen::Matrix<double,-1,1>::Constant(num_to_write,
             std::numeric_limits<double>::quiet_NaN());
    write_array_impl(base_rng, params_r, params_i, vars,
      emit_transformed_parameters, emit_generated_quantities, pstream);
  }
  template <typename RNG> inline void
  write_array(RNG& base_rng, std::vector<double>& params_r, std::vector<int>&
              params_i, std::vector<double>& vars, bool
              emit_transformed_parameters = true, bool
              emit_generated_quantities = true, std::ostream*
              pstream = nullptr) const {
    const size_t num_params__ = (NC * NR);
    const size_t num_transformed = emit_transformed_parameters * ((NC * NR));
    const size_t num_gen_quantities = emit_generated_quantities * ((NC * NR));
    const size_t num_to_write = num_params__ + num_transformed +
      num_gen_quantities;
    vars = std::vector<double>(num_to_write,
             std::numeric_limits<double>::quiet_NaN());
    write_array_impl(base_rng, params_r, params_i, vars,
      emit_transformed_parameters, emit_generated_quantities, pstream);
  }
  template <bool propto__, bool jacobian__, typename T_> inline T_
  log_prob(Eigen::Matrix<T_,-1,1>& params_r, std::ostream* pstream = nullptr) const {
    Eigen::Matrix<int,-1,1> params_i;
    return log_prob_impl<propto__, jacobian__>(params_r, params_i, pstream);
  }
  template <bool propto__, bool jacobian__, typename T_> inline T_
  log_prob(std::vector<T_>& params_r, std::vector<int>& params_i,
           std::ostream* pstream = nullptr) const {
    return log_prob_impl<propto__, jacobian__>(params_r, params_i, pstream);
  }
  inline void
  transform_inits(const stan::io::var_context& context,
                  Eigen::Matrix<double,-1,1>& params_r, std::ostream*
                  pstream = nullptr) const final {
    std::vector<double> params_r_vec(params_r.size());
    std::vector<int> params_i;
    transform_inits(context, params_i, params_r_vec, pstream);
    params_r = Eigen::Map<Eigen::Matrix<double,-1,1>>(params_r_vec.data(),
                 params_r_vec.size());
  }
  inline void
  transform_inits(const stan::io::var_context& context, std::vector<int>&
                  params_i, std::vector<double>& vars, std::ostream*
                  pstream__ = nullptr) const {
    vars.resize(num_params_r__);
    transform_inits_impl(context, vars, pstream__);
  }
  inline void
  unconstrain_array(const std::vector<double>& params_constrained,
                    std::vector<double>& params_unconstrained, std::ostream*
                    pstream = nullptr) const {
    const std::vector<int> params_i;
    params_unconstrained = std::vector<double>(num_params_r__,
                             std::numeric_limits<double>::quiet_NaN());
    unconstrain_array_impl(params_constrained, params_i,
      params_unconstrained, pstream);
  }
  inline void
  unconstrain_array(const Eigen::Matrix<double,-1,1>& params_constrained,
                    Eigen::Matrix<double,-1,1>& params_unconstrained,
                    std::ostream* pstream = nullptr) const {
    const std::vector<int> params_i;
    params_unconstrained = Eigen::Matrix<double,-1,1>::Constant(num_params_r__,
                             std::numeric_limits<double>::quiet_NaN());
    unconstrain_array_impl(params_constrained, params_i,
      params_unconstrained, pstream);
  }
};
}
using stan_model = model_bdcF_namespace::model_bdcF;
#ifndef USING_R
// Boilerplate
stan::model::model_base&
new_model(stan::io::var_context& data_context, unsigned int seed,
          std::ostream* msg_stream) {
  stan_model* m = new stan_model(data_context, seed, msg_stream);
  return *m;
}
stan::math::profile_map& get_stan_profile_data() {
  return model_bdcF_namespace::profiles__;
}
#endif
#endif
