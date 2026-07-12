// Dump C++ rrtmgp_main for BFB validation of the JAX port.
// Input: ncol nlay, then play, plev, tlay, tlev, tsfc(unused), vmr(8),
// then lwp, iwp, rel, rei, cldfrac (ncol*nlay each), sfc_alb_dir_vis,
// dir_nir, dif_vis, dif_nir, mu0 (ncol each). Aerosols zero. tsi=1.
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
    auto rd1 = [&](int n1) {
      Kokkos::View<double*, Kokkos::HostSpace> h("h", n1);
      for (int i = 0; i < n1; ++i) in >> h(i);
      real1dk d("d", n1);
      Kokkos::deep_copy(d, h);
      return d;
    };
    auto play = rd2(ncol, nlay);
    auto plev = rd2(ncol, nlay+1);
    auto tlay = rd2(ncol, nlay);
    auto tlev = rd2(ncol, nlay+1);
    auto tsfc = rd1(ncol);
    std::vector<std::string> gases = {"h2o","co2","o3","n2o","co","ch4","o2","n2"};
    gas_concs_t gas_concs;
    string1dv gas_names_v(gases);
    gas_concs.init(gas_names_v, ncol, nlay);
    for (int ig = 0; ig < 8; ++ig) {
      auto vmr = rd2(ncol, nlay);
      gas_concs.set_vmr(gases[ig], vmr);
    }
    auto lwp = rd2(ncol, nlay);
    auto iwp = rd2(ncol, nlay);
    auto rel = rd2(ncol, nlay);
    auto rei = rd2(ncol, nlay);
    auto cld = rd2(ncol, nlay);
    auto adv = rd1(ncol);
    auto anv = rd1(ncol);
    auto afv = rd1(ncol);
    auto afn = rd1(ncol);
    auto mu0 = rd1(ncol);

    interface_t::rrtmgp_initialize(gas_concs,
      "/work/e3sm-inputdata/atm/scream/init/rrtmgp-data-sw-g112-210809.nc",
      "/work/e3sm-inputdata/atm/scream/init/rrtmgp-data-lw-g128-210809.nc",
      "/work/e3sm-inputdata/atm/scream/init/rrtmgp-cloud-optics-coeffs-sw.nc",
      "/work/e3sm-inputdata/atm/scream/init/rrtmgp-cloud-optics-coeffs-lw.nc",
      nullptr, 8.0);

    const int nswb = interface_t::k_dist_sw_k->get_nband();
    const int nlwb = interface_t::k_dist_lw_k->get_nband();
    const int nswg = interface_t::k_dist_sw_k->get_ngpt();
    const int nlwg = interface_t::k_dist_lw_k->get_ngpt();

    auto sfc_alb_dir = pool_t::alloc<double>(ncol, nswb);
    auto sfc_alb_dif = pool_t::alloc<double>(ncol, nswb);
    interface_t::compute_band_by_band_surface_albedos(ncol, nswb,
      adv, anv, afv, afn, sfc_alb_dir, sfc_alb_dif);

    auto z3 = [&](int a, int b, int c) {
      auto v = pool_t::alloc<double>(a, b, c);
      Kokkos::deep_copy(v, 0.0);
      return v;
    };
    auto z2 = [&](int a, int b) {
      auto v = pool_t::alloc<double>(a, b);
      Kokkos::deep_copy(v, 0.0);
      return v;
    };
    auto aer_tau_sw = z3(ncol, nlay, nswb);
    auto aer_ssa_sw = z3(ncol, nlay, nswb);
    auto aer_g_sw   = z3(ncol, nlay, nswb);
    auto aer_tau_lw = z3(ncol, nlay, nlwb);
    auto cld_tau_sw_bnd = z3(ncol, nlay, nswb);
    auto cld_tau_lw_bnd = z3(ncol, nlay, nlwb);
    auto cld_tau_sw_gpt = z3(ncol, nlay, nswg);
    auto cld_tau_lw_gpt = z3(ncol, nlay, nlwg);
    auto swu = z2(ncol,nlay+1); auto swd = z2(ncol,nlay+1); auto swdd = z2(ncol,nlay+1);
    auto lwu = z2(ncol,nlay+1); auto lwd = z2(ncol,nlay+1);
    auto ccsu = z2(ncol,nlay+1); auto ccsd = z2(ncol,nlay+1); auto ccsdd = z2(ncol,nlay+1);
    auto csu = z2(ncol,nlay+1); auto csd = z2(ncol,nlay+1); auto csdd = z2(ncol,nlay+1);
    auto cnu = z2(ncol,nlay+1); auto cnd = z2(ncol,nlay+1); auto cndd = z2(ncol,nlay+1);
    auto lccu = z2(ncol,nlay+1); auto lccd = z2(ncol,nlay+1);
    auto lcu = z2(ncol,nlay+1); auto lcd = z2(ncol,nlay+1);
    auto lnu = z2(ncol,nlay+1); auto lnd = z2(ncol,nlay+1);
    auto swbu = z3(ncol,nlay+1,nswb); auto swbd = z3(ncol,nlay+1,nswb); auto swbdd = z3(ncol,nlay+1,nswb);
    auto lwbu = z3(ncol,nlay+1,nlwb); auto lwbd = z3(ncol,nlay+1,nlwb);

    interface_t::rrtmgp_main(ncol, nlay, play, tlay, plev, tlev, gas_concs,
      sfc_alb_dir, sfc_alb_dif, mu0, lwp, iwp, rel, rei, cld,
      aer_tau_sw, aer_ssa_sw, aer_g_sw, aer_tau_lw,
      cld_tau_sw_bnd, cld_tau_lw_bnd, cld_tau_sw_gpt, cld_tau_lw_gpt,
      swu, swd, swdd, lwu, lwd,
      ccsu, ccsd, ccsdd, csu, csd, csdd, cnu, cnd, cndd,
      lccu, lccd, lcu, lcd, lnu, lnd,
      swbu, swbd, swbdd, lwbu, lwbd,
      1.0, nullptr, true, true);

    auto pr2 = [&](const real2dk& v) {
      auto h = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), v);
      for (size_t i=0;i<h.extent(0);++i) for (size_t k=0;k<h.extent(1);++k)
        std::printf("%.17e\n", h(i,k));
    };
    auto pr3 = [&](const real3dk& v) {
      auto h = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), v);
      for (size_t i=0;i<h.extent(0);++i) for (size_t k=0;k<h.extent(1);++k)
        for (size_t g=0;g<h.extent(2);++g) std::printf("%.17e\n", h(i,k,g));
    };
    pr2(swu); pr2(swd); pr2(swdd); pr2(lwu); pr2(lwd);
    pr2(ccsu); pr2(ccsd); pr2(ccsdd); pr2(csu); pr2(csd); pr2(csdd);
    pr2(cnu); pr2(cnd); pr2(cndd);
    pr2(lccu); pr2(lccd); pr2(lcu); pr2(lcd); pr2(lnu); pr2(lnd);
    pr3(cld_tau_sw_bnd); pr3(cld_tau_lw_bnd);
    pr3(cld_tau_sw_gpt); pr3(cld_tau_lw_gpt);
    pr3(swbu); pr3(swbd); pr3(swbdd); pr3(lwbu); pr3(lwbd);
    pr2(sfc_alb_dir); pr2(sfc_alb_dif);
  }
  scream::finalize_kls();
  Kokkos::finalize();
  return 0;
}
