/* Copyright (C) 2019 Pablo Hernandez-Cerdan
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

#include "perm_common_py.hpp"
#include "perm_common_types.hpp"
#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <sstream>

namespace py = pybind11;
using namespace perm;

void init_perm_common_types(py::module &m) {
    py::class_<vec3D_t<int>>(m, "vec3D")
            .def(py::init())
            .def(py::init([](const int &x, const int &y, const int &z) {
                return vec3D_t<int>{x, y, z};
            }))
            .def_readwrite("x", &vec3D_t<int>::x)
            .def_readwrite("y", &vec3D_t<int>::y)
            .def_readwrite("z", &vec3D_t<int>::z)
            .def("__add__",
                 [](const vec3D_t<int> &lhs, const vec3D_t<int> &rhs) {
                     return perm::plus(lhs, rhs);
                 })
            .def("__eq__", [](const vec3D_t<int> &lhs,
                              const vec3D_t<int> &rhs) { return lhs == rhs; })
            .def("__getitem__",
                 [](const vec3D_t<int> &vec, size_t i) {
                     if (i >= 3)
                         throw py::index_error();
                     return vec[i];
                 })
            .def("__setitem__",
                 [](vec3D_t<int> &vec, size_t i, int new_value) {
                     if (i >= 3)
                         throw py::index_error();
                     vec[i] = new_value;
                 })
            .def("__len__", &vec3D_t<int>::size)
            .def("__array__",
                 [](const vec3D_t<int> &vec) {
                     py::array_t<int> out(3);
                     auto r = out.mutable_unchecked<1>();
                     r(0) = vec.x;
                     r(1) = vec.y;
                     r(2) = vec.z;
                     return out;
                 })
            .def("__repr__", [](const vec3D_t<int> &vec) {
                std::stringstream os;
                vec.print(os);
                return os.str();
            });
    py::bind_vector<std::vector<vec3D_t<int>>>(m, "points_vector",
                                               py::module_local(false));
    py::class_<single_chain_t<int>>(m, "single_chain")
            .def(py::init())
            .def_readwrite("monomers", &single_chain_t<int>::monomers)
            .def_readwrite("points", &single_chain_t<int>::points)
            .def("__repr__", [](const single_chain_t<int> &chain) {
                std::stringstream os;
                chain.print(os);
                return os.str();
            });
}
