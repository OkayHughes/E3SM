#include "physics/p3/eamxx_p3_process_interface.hpp"

#include <ekat_team_policy_utils.hpp>

#ifdef EAMXX_HAS_PYTHON
#include "share/atm_process/atmosphere_process_pyhelpers.hpp"
#include <pybind11/numpy.h>
#endif

namespace scream {

void P3Microphysics::run_impl (const double dt)
{
  using TPF = ekat::TeamPolicyFactory<KT::ExeSpace>;

  // Set the dt for p3 postprocessing
  p3_postproc.m_dt = dt;

#ifdef EAMXX_HAS_PYTHON
  if (has_py_module()) {
    // Swap the ENTIRE P3 step (pre-process + p3_main + post-process) with
    // the python implementation; this matches the validation unit of the
    // JAX port (see jax_port/TEST_HARNESS_DESIGN.md and
    // scream_jax/p3/process.py).
    EKAT_REQUIRE_MSG (m_params.get<std::string>("py_backend","host")=="host",
        "Error! The P3 python bridge only supports the host backend.\n");
    EKAT_REQUIRE_MSG (not infrastructure.prescribedCCN and
                      not runtime_options.use_hetfrz_classnuc and
                      not runtime_options.use_separate_ice_liq_frac and
                      not runtime_options.extra_p3_diags and
                      not has_column_conservation_check(),
        "Error! The P3 python bridge does not marshal the fields of "
        "do_prescribed_ccn/use_hetfrz_classnuc/use_separate_ice_liq_frac/"
        "extra_p3_diags/column conservation checks.\n");

    for (auto n : {"T_mid","qv","qc","nc","qr","nr","qi","ni","qm","bm",
                   "qv_prev_micro_step","T_prev_micro_step",
                   "precip_liq_surf_mass","precip_ice_surf_mass"}) {
      get_field_out(n).sync_to_host();
    }
    for (auto n : {"p_mid","p_dry_mid","pseudo_density","pseudo_density_dry",
                   "cldfrac_tot","nc_nuceat_tend","ni_activated",
                   "inv_qc_relvar"}) {
      get_field_in(n).sync_to_host();
    }

    py_module_call("main",
        dt,
        infrastructure.predictNc,
        runtime_options.do_ice_production,
        runtime_options.set_cld_frac_l_to_one,
        runtime_options.set_cld_frac_i_to_one,
        runtime_options.set_cld_frac_r_to_one,
        // P3Runtime scalars (order matches scream_jax.p3.DEFAULT_OPTS)
        runtime_options.max_total_ni,
        runtime_options.autoconversion_prefactor,
        runtime_options.autoconversion_qc_exponent,
        runtime_options.autoconversion_nc_exponent,
        runtime_options.autoconversion_radius,
        runtime_options.accretion_prefactor,
        runtime_options.accretion_qc_exponent,
        runtime_options.accretion_qr_exponent,
        runtime_options.rain_selfcollection_prefactor,
        runtime_options.rain_selfcollection_breakup_diameter,
        runtime_options.constant_mu_rain,
        runtime_options.spa_ccn_to_nc_factor,
        runtime_options.spa_ccn_to_nc_exponent,
        runtime_options.cldliq_to_ice_collection_factor,
        runtime_options.rain_to_ice_collection_factor,
        runtime_options.min_rime_rho,
        runtime_options.max_rime_rho,
        runtime_options.immersion_freezing_exponent,
        runtime_options.deposition_nucleation_exponent,
        runtime_options.ice_sedimentation_factor,
        // inputs
        get_py_field_host("p_mid"), get_py_field_host("p_dry_mid"),
        get_py_field_host("pseudo_density"),
        get_py_field_host("pseudo_density_dry"),
        get_py_field_host("cldfrac_tot"),
        get_py_field_host("nc_nuceat_tend"),
        get_py_field_host("ni_activated"),
        get_py_field_host("inv_qc_relvar"),
        // updated
        get_py_field_host("T_mid"), get_py_field_host("qv"),
        get_py_field_host("qc"), get_py_field_host("nc"),
        get_py_field_host("qr"), get_py_field_host("nr"),
        get_py_field_host("qi"), get_py_field_host("ni"),
        get_py_field_host("qm"), get_py_field_host("bm"),
        get_py_field_host("qv_prev_micro_step"),
        get_py_field_host("T_prev_micro_step"),
        get_py_field_host("precip_liq_surf_mass"),
        get_py_field_host("precip_ice_surf_mass"),
        // computed
        get_py_field_host("eff_radius_qc"),
        get_py_field_host("eff_radius_qi"),
        get_py_field_host("eff_radius_qr"),
        get_py_field_host("precip_total_tend"),
        get_py_field_host("nevapr"),
        get_py_field_host("diag_equiv_reflectivity"),
        get_py_field_host("micro_liq_ice_exchange"),
        get_py_field_host("micro_vap_liq_exchange"),
        get_py_field_host("micro_vap_ice_exchange"),
        get_py_field_host("rainfrac"));

    for (auto n : {"T_mid","qv","qc","nc","qr","nr","qi","ni","qm","bm",
                   "qv_prev_micro_step","T_prev_micro_step",
                   "precip_liq_surf_mass","precip_ice_surf_mass",
                   "eff_radius_qc","eff_radius_qi","eff_radius_qr",
                   "precip_total_tend","nevapr","diag_equiv_reflectivity",
                   "micro_liq_ice_exchange","micro_vap_liq_exchange",
                   "micro_vap_ice_exchange","rainfrac"}) {
      get_field_out(n).sync_to_dev();
    }
    infrastructure.it++;
    return;
  }
#endif

  // Create policy for pre and post process pfor
  const auto nlev_packs  = ekat::npack<Pack>(m_num_levs);
  const auto policy = TPF::get_default_team_policy(m_num_cols, nlev_packs);

  // Assign values to local arrays used by P3, these are now stored in p3_loc.
  Kokkos::parallel_for(
    "p3_pre_process",
    policy,
    p3_preproc
  );
  Kokkos::fence();

  // Update the variables in the p3 input structures with local values.

  infrastructure.dt = dt;
  infrastructure.it++;

  // Reset internal WSM variables.
  workspace_mgr.reset_internals();

  // Run p3 main
  get_field_out("micro_liq_ice_exchange").deep_copy(0.0);
  get_field_out("micro_vap_liq_exchange").deep_copy(0.0);
  get_field_out("micro_vap_ice_exchange").deep_copy(0.0);

  // Optional extra p3 diags
  if (runtime_options.extra_p3_diags) {
    get_field_out("qr2qv_evap").deep_copy(0.0);
    get_field_out("qi2qv_sublim").deep_copy(0.0);
    get_field_out("qc2qr_accret").deep_copy(0.0);
    get_field_out("qc2qr_autoconv").deep_copy(0.0);
    get_field_out("qv2qi_vapdep").deep_copy(0.0);
    get_field_out("qc2qi_berg").deep_copy(0.0);
    get_field_out("qc2qr_ice_shed").deep_copy(0.0);
    get_field_out("qc2qi_collect").deep_copy(0.0);
    get_field_out("qr2qi_collect").deep_copy(0.0);
    get_field_out("qc2qi_hetero_freeze").deep_copy(0.0);
    get_field_out("qr2qi_immers_freeze").deep_copy(0.0);
    get_field_out("qi2qr_melt").deep_copy(0.0);
    get_field_out("qr_sed").deep_copy(0.0);
    get_field_out("qc_sed").deep_copy(0.0);
    get_field_out("qi_sed").deep_copy(0.0);
  }

  P3F::p3_main(runtime_options, prog_state, diag_inputs, diag_outputs, infrastructure,
               history_only, lookup_tables,
#ifdef SCREAM_P3_SMALL_KERNELS
               temporaries,
#endif
               workspace_mgr, m_num_cols, m_num_levs);

  // Conduct the post-processing of the p3_main output.
  Kokkos::parallel_for(
    "p3_post_process",
    policy,
    p3_postproc
  );
  Kokkos::fence();
}

} // namespace scream
