/* Copyright (C) 2019 Pablo Hernandez-Cerdan
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */
#ifndef DECLARE_TYPES_PY_HPP
#define DECLARE_TYPES_PY_HPP

#include "vec3D.hpp"
#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>

namespace perm {
template <typename T>
void declare_vec3D(pybind11::module &m, const std::string &typestr) {
    namespace py = pybind11;
    using TVec3D = vec3D_t<T>;
    const std::string pyclass_name = std::string("vec3D") + typestr;
    py::class_<TVec3D>(m, pyclass_name.c_str())
            .def(py::init())
            .def(py::init([](const T &x, const T &y, const T &z) {
                return TVec3D{x, y, z};
            }))
            .def_property(
                    "x", [](const TVec3D &lhs) { return lhs.x(); },
                    [](TVec3D &lhs, const T &value) {
                        return lhs.x() = value;
                    })
            .def_property(
                    "y", [](const TVec3D &lhs) { return lhs.y(); },
                    [](TVec3D &lhs, const T &value) {
                        return lhs.y() = value;
                    })
            .def_property(
                    "z", [](const TVec3D &lhs) { return lhs.z(); },
                    [](TVec3D &lhs, const T &value) {
                        return lhs.z() = value;
                    })
            .def("__add__",
                 [](const TVec3D &lhs, const TVec3D &rhs) {
                     return perm::plus(lhs, rhs);
                 })
            .def("__eq__", [](const TVec3D &lhs,
                              const TVec3D &rhs) { return lhs == rhs; })
            .def("__getitem__",
                 [](const TVec3D &vec, size_t i) {
                     if (i >= 3)
                         throw py::index_error();
                     return vec[i];
                 })
            .def("__setitem__",
                 [](TVec3D &vec, size_t i, T new_value) {
                     if (i >= 3)
                         throw py::index_error();
                     vec[i] = new_value;
                 })
            .def("__len__", &TVec3D::size)
            .def("__array__",
                 [](const TVec3D &vec) {
                     py::array_t<T> out(3);
                     auto r = out.template mutable_unchecked<1>();
                     r(0) = vec.e[0];
                     r(1) = vec.e[1];
                     r(2) = vec.e[2];
                     return out;
                 })
            .def("__repr__", [](const TVec3D &vec) {
                std::stringstream os;
                vec.print(os);
                return os.str();
            });
}
} // end namespace perm
#endif
