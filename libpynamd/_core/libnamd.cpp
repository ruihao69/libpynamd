#include <pybind11/pybind11.h>
#include <pybind11/eigen.h>
#include <pybind11/numpy.h>
#include "libnamd.hpp"

using namespace rhbi::libnamd;

// ----------------
// Python interface
// ----------------

namespace py = pybind11;

PYBIND11_MODULE(_core, m){
  // ----------------
  // module descriptions 
  // ----------------
  m.doc() = "libnamd: a simple c++ function library for Non-Adiabatic Molecular Dynamics";

  // ----------------
  // general operator functions
  // ----------------


}