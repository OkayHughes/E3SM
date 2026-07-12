// Dump C++ RRTMGP gas optics (SW + LW) for BFB validation of the JAX port.
// Input text file: ncol nlay, then play(ncol*nlay), plev(ncol*(nlay+1)),
// tlay, tlev, tsfc(ncol), vmr(ngas=8, ncol*nlay each, order:
// h2o co2 o3 n2o co ch4 o2 n2).
#include "physics/rrtmgp/eamxx_rrtmgp_interface.hpp"
#include <fstream>
#include <cstdio>

using namespace scream;
using interface_t = rrtmgp::rrtmgp_interface<>;
using pool_t = interface_t::pool_t;
using gas_concs_t = interface_t::gas_concs_t;
using real1dk = interface_t::real1dk;
using real2dk = interface_t::real2dk;
using real3dk = interface_t::real3dk;
using optical_props1_t = interface_t::optical_props1_t;
using optical_props2_t = interface_t::optical_props2_t;
using source_func_t = interface_t::source_func_t;

int main(int argc, char** argv) {
  Kokkos::initialize(argc, argv);
  scream::init_kls();
  {
    std::ifstream in(argv[1]);
    int ncol, nlay; in >> ncol >> nlay;
    auto rd2 = [&](int n1, int n2) {
      Kokkos::View<double**, Kokkos::LayoutRight, Kokkos::HostSpace> h("h", n1, n2);
      for (int i = 0; i < n1; ++i) for (int k = 0; k < n2; ++k) in >> h(i,k);
      real2dk d("d", n1, n2);
      Kokkos::deep_copy(d, h);
      return d;
    };
    auto play = rd2(ncol, nlay);
    auto plev = rd2(ncol, nlay+1);
    auto tlay = rd2(ncol, nlay);
    auto tlev = rd2(ncol, nlay+1);
    Kokkos::View<double*, Kokkos::HostSpace> tsfc_h("tsfc", ncol);
    for (int i = 0; i < ncol; ++i) in >> tsfc_h(i);
    real1dk tsfc("tsfc", ncol);
    Kokkos::deep_copy(tsfc, tsfc_h);

    std::vector<std::string> gases = {"h2o","co2","o3","n2o","co","ch4","o2","n2"};
    gas_concs_t gas_concs;
    string1dv gas_names_v(gases);
    gas_concs.init(gas_names_v, ncol, nlay);
    for (int ig = 0; ig < 8; ++ig) {
      auto vmr = rd2(ncol, nlay);
      gas_concs.set_vmr(gases[ig], vmr);
    }

    interface_t::rrtmgp_initialize(gas_concs,
      "/work/e3sm-inputdata/atm/scream/init/rrtmgp-data-sw-g112-210809.nc",
      "/work/e3sm-inputdata/atm/scream/init/rrtmgp-data-lw-g128-210809.nc",
      "/work/e3sm-inputdata/atm/scream/init/rrtmgp-cloud-optics-coeffs-sw.nc",
      "/work/e3sm-inputdata/atm/scream/init/rrtmgp-cloud-optics-coeffs-lw.nc",
      nullptr, 4.0);

    auto& ksw = *interface_t::k_dist_sw_k;
    auto& klw = *interface_t::k_dist_lw_k;
    const int ngpt_sw = ksw.get_ngpt(), nband_sw = ksw.get_nband();
    const int ngpt_lw = klw.get_ngpt(), nband_lw = klw.get_nband();

    bool top_at_1 = true; // p increases with index in our input

    // ---- SW ----
    {
      auto b2g = pool_t::alloc<int>(2, nband_sw);
      auto g2b = pool_t::alloc<int>(ngpt_sw);
      auto tau_m = pool_t::alloc<double>(ncol, nlay, ngpt_sw);
      auto ssa_m = pool_t::alloc<double>(ncol, nlay, ngpt_sw);
      auto g_m   = pool_t::alloc<double>(ncol, nlay, ngpt_sw);
      auto col_gas = pool_t::alloc<double>(ncol, nlay, ksw.get_ngas()+1);
      auto toa = pool_t::alloc<double>(ncol, ngpt_sw);
      optical_props2_t optics;
      optics.alloc_2str_no_alloc(ncol, nlay, ksw, b2g, g2b, tau_m, ssa_m, g_m);
      ksw.gas_optics(ncol, nlay, top_at_1, play, plev, tlay, gas_concs, col_gas, optics, toa);
      auto tau_h = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), optics.tau);
      auto ssa_h = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), optics.ssa);
      auto g_h   = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), optics.g);
      auto toa_h = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), toa);
      std::printf("SW\n");
      for (int i=0;i<ncol;i++) for (int k=0;k<nlay;k++) for (int g=0;g<ngpt_sw;g++)
        std::printf("%.17e %.17e %.17e\n", tau_h(i,k,g), ssa_h(i,k,g), g_h(i,k,g));
      for (int i=0;i<ncol;i++) for (int g=0;g<ngpt_sw;g++) std::printf("%.17e\n", toa_h(i,g));
    }

    // ---- LW ----
    {
      auto b2g = pool_t::alloc<int>(2, nband_lw);
      auto g2b = pool_t::alloc<int>(ngpt_lw);
      auto tau_m = pool_t::alloc<double>(ncol, nlay, ngpt_lw);
      auto col_gas = pool_t::alloc<double>(ncol, nlay, klw.get_ngas()+1);
      auto sb2g = pool_t::alloc<int>(2, nband_lw);
      auto sg2b = pool_t::alloc<int>(ngpt_lw);
      auto sfc_m = pool_t::alloc<double>(ncol, ngpt_lw);
      auto lay_m = pool_t::alloc<double>(ncol, nlay, ngpt_lw);
      auto inc_m = pool_t::alloc<double>(ncol, nlay, ngpt_lw);
      auto dec_m = pool_t::alloc<double>(ncol, nlay, ngpt_lw);
      optical_props1_t optics;
      optics.alloc_1scl_no_alloc(ncol, nlay, klw, b2g, g2b, tau_m);
      source_func_t sources;
      sources.alloc_no_alloc(ncol, nlay, klw, sb2g, sg2b, sfc_m, lay_m, inc_m, dec_m);
      klw.gas_optics(ncol, nlay, top_at_1, play, plev, tlay, tsfc, gas_concs, col_gas, optics,
                     sources, interface_t::view_t<double**>(), tlev);
      auto tau_h = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), optics.tau);
      auto sfc_h = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), sources.sfc_source);
      auto lay_h = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), sources.lay_source);
      auto inc_h = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), sources.lev_source_inc);
      auto dec_h = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), sources.lev_source_dec);
      std::printf("LW\n");
      for (int i=0;i<ncol;i++) for (int k=0;k<nlay;k++) for (int g=0;g<ngpt_lw;g++)
        std::printf("%.17e %.17e %.17e %.17e\n", tau_h(i,k,g), lay_h(i,k,g), inc_h(i,k,g), dec_h(i,k,g));
      for (int i=0;i<ncol;i++) for (int g=0;g<ngpt_lw;g++) std::printf("%.17e\n", sfc_h(i,g));
    }
  }
  scream::finalize_kls();
  Kokkos::finalize();
  return 0;
}
