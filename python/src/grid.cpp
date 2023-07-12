#include <optional>
#include <variant>

#include "pybind11/pybind11.h"
#include "pybind11/stl.h"

#include "xtensor-python/pytensor.hpp"
#include "xtensor-python/pyarray.hpp"

#include "fastscapelib/grid/base.hpp"
#include "fastscapelib/grid/profile_grid.hpp"
#include "fastscapelib/grid/raster_grid.hpp"

#include "grid.hpp"
#include "pytensor_utils.hpp"


namespace py = pybind11;
namespace fs = fastscapelib;


void
add_grid_bindings(py::module& m)
{
    py::options options;
    options.disable_function_signatures();

    // ==== Binding of the node_status enumeration ==== //
    py::enum_<fs::node_status> node_status(m,
                                           "NodeStatus",
                                           py::arithmetic(),
                                           R"doc(
        Node status values.

        Node status is a label assigned to each grid/mesh node. It is used to
        define the structure and boundaries of the modeled domain.

        The description of each value is given for indicative pruposes only.
        Exact semantics may differ depending on the grid type or the objects that
        are consuming grids.

        For example, ``LOOPED`` is used only for structured grids.
        )doc");

    node_status.value("CORE", fs::node_status::core, "Inner grid node")
        .value("FIXED_VALUE", fs::node_status::fixed_value, "boundary node with fixed value")
        .value(
            "FIXED_GRADIENT", fs::node_status::fixed_gradient, "boundary node with fixed gradient")
        .value("LOOPED", fs::node_status::looped, "reflective (periodic) boundary node");

    // ==== Binding of the node structure ==== //
    py::class_<fs::node> node(m, "Node", "Represents a grid/mesh node.");

    node.def(py::init<std::size_t, fs::node_status>(),
             py::arg("idx"),
             py::arg("status"),
             R"doc(__init__(self, idx: int, status: NodeStatus) -> None

        Node initializer (internal use only).

        Parameters
        ----------
        idx : int
            Node index.
        status : :class:`NodeStatus`
            Node status.

        )doc");

    node.def_readwrite("idx", &fs::node::idx, "Node index.")
        .def_readwrite("status", &fs::node::status, "Node status.");

    // ==== Binding of the neighbor structure ==== //
    py::class_<fs::neighbor> neighbor(m, "Neighbor", "Represents a grid/mesh node neighbor.");

    neighbor.def(py::init<std::size_t, double, fs::node_status>(),
                 py::arg("idx"),
                 py::arg("distance"),
                 py::arg("status"),
                 R"doc(__init__(self, idx: int, distance: float, status: NodeStatus) -> None

        Neighbor initializer (internal use only).

        Parameters
        ----------
        idx : int
            Neighbor index.
        distance : float
            Neighbor distance.
        status : :class:`NodeStatus`
            Neighbor status.

        )doc");

    neighbor.def("__eq__", &fs::neighbor::operator==)
        .def_readwrite("idx", &fs::neighbor::idx, "Neighbor index.")
        .def_readwrite("distance", &fs::neighbor::distance, "Neighbor distance.")
        .def_readwrite("status", &fs::neighbor::status, "Neighbor status.");

    /*
    ** Triangular Mesh
    */

    // ==== Binding of the TriMesh class ==== //
    py::class_<fs::py_trimesh> tmesh(m, "TriMesh", "A 2-dimensional, triangular (irregular) mesh.");

    using trimesh_nstatus_type = std::optional<fs::py_trimesh::nodes_status_map_type>;

    tmesh.def(py::init(
                  [](const fs::py_trimesh::points_type& points,
                     const fs::py_trimesh::triangles_type& triangles,
                     const trimesh_nstatus_type& nodes_status)
                  {
                      if (!nodes_status.has_value())
                      {
                          return fs::py_trimesh(points, triangles);
                      }
                      else
                      {
                          const auto& nstatus = nodes_status.value();
                          return fs::py_trimesh(points, triangles, nstatus);
                      }
                  }),
              py::arg("points"),
              py::arg("triangles"),
              py::arg("nodes_status") = py::none(),
              R"doc(__init__(*args, **kwargs) -> None

        TriMesh initializer (from an existing triangulation).

        Overloaded initializer.

        1. ``__init__(self, points: numpy.ndarray, triangles: numpy.ndarray, nodes_status: Dict[int, NodeStatus] | None = None) -> None``

        Create a new mesh with a dictionary of node status (optional).

        2. ``__init__(self, points: numpy.ndarray, triangles: numpy.ndarray, nodes_status: np.ndarray) -> None``

        Create a new mesh with an array of node status.

        Parameters
        ----------
        points : array-like
            Mesh node x,y coordinates (array of shape [N, 2]).
        triangles : array-like
            Indices of triangle vertices (array of shape [K, 3]).
        nodes_status : dict or array-like, optional
            A dictionary where keys are node indices and values are
            :class:`NodeStatus` values for setting the status at any
            node on the mesh. If ``None`` (default), "fixed value" is set
            for all boundary nodes (i.e., end-points of all the edges that
            are not shared by more than one triangle). Alternatively, an
            array of shape [N] can be given for setting the status of all
            nodes at once.

        )doc");
    tmesh.def(py::init<const fs::py_trimesh::points_type&,
                       const fs::py_trimesh::triangles_type&,
                       const fs::py_trimesh::nodes_status_array_type&>(),
              py::arg("points"),
              py::arg("triangles"),
              py::arg("nodes_status"));

    fs::register_grid_static_properties(tmesh);
    fs::register_base_grid_properties(tmesh);
    fs::register_grid_methods(tmesh);

    /*
    ** Profile Grid
    */

    // ==== Binding of the profile_boundary_status class ==== //
    py::class_<fs::profile_boundary_status> profile_bstatus(
        m, "ProfileBoundaryStatus", "Status at profile grid boundary nodes.");

    profile_bstatus.def(py::init<const fs::node_status>(),
                        py::arg("status"),
                        R"doc(__init__(*args, **kwargs) -> None

        Overloaded initializer.

        1. ``__init__(self, status: NodeStatus) -> None``

        Set the same status for the left/right boundary nodes.

        status : :class:`NodeStatus`
            Boundary node status.

        )doc");
    profile_bstatus.def(py::init<const fs::node_status, const fs::node_status>(),
                        py::arg("left_status"),
                        py::arg("right_status"),
                        R"doc(

        2. ``__init__(self, left_status: NodeStatus, right_status: NodeStatus) -> None``

        Set status for the left and right boundary nodes.

        left_status : :class:`NodeStatus`
            Left boundary node status.
        right_status : :class:`NodeStatus`
            Right boundary node status.

        )doc");
    profile_bstatus.def(py::init<const std::array<fs::node_status, 2>&>(),
                        py::arg("status"),
                        R"doc(

        3. ``__init__(self, status: List[NodeStatus]) -> None``

        Set status for the left and right boundary nodes.

        status : list of :class:`NodeStatus`
            2-length list of left and right boundary node status.

        )doc");

    profile_bstatus
        .def_readwrite("left", &fs::profile_boundary_status::left, "Left boundary status.")
        .def_readwrite("right", &fs::profile_boundary_status::right, "Right boundary status.")
        .def_property_readonly("is_horizontal_looped",
                               &fs::profile_boundary_status::is_horizontal_looped,
                               "Horizontal looped status.");

    // ==== Binding of the ProfileGrid class ==== //
    py::class_<fs::py_profile_grid> pgrid(m, "ProfileGrid", "A 1-dimensional profile grid.");

    using profile_bstatus_type = std::
        variant<fs::node_status, std::array<fs::node_status, 2>, fs::profile_boundary_status>;
    using profile_nstatus_type = std::optional<fs::py_profile_grid::nodes_status_map_type>;

    pgrid.def(
        py::init(
            [](fs::py_profile_grid::size_type size,
               fs::py_profile_grid::spacing_type spacing,
               const profile_bstatus_type& bounds_status,
               const profile_nstatus_type& nodes_status)
            {
                auto bstatus = std::visit([](auto&& bs) { return fs::profile_boundary_status(bs); },
                                          bounds_status);

                if (!nodes_status.has_value())
                {
                    return fs::py_profile_grid(size, spacing, bstatus);
                }
                else
                {
                    const auto& nstatus = nodes_status.value();
                    return fs::py_profile_grid(size, spacing, bstatus, nstatus);
                }
            }),
        py::arg("size"),
        py::arg("spacing"),
        py::arg("bounds_status"),
        py::arg("nodes_status") = py::none(),
        R"doc(__init__(self, size: int, spacing: float, bounds_status: NodeStatus | List[NodeStatus] | ProfileBoundaryStatus, nodes_status: Dict[int, NodeStatus] | None) -> None

        Profile grid initializer (overloaded).

        Parameters
        ----------
        size : int
            Total number of grid nodes.
        spacing : float
            Distance between two adjacent grid nodes.
        bounds_status : :class:`NodeStatus` or list or :class:`ProfileBoundaryStatus`
            Status at boundary nodes (left/right grid edges).
        nodes_status : dict, optional
            A dictionary where keys are node indices and values are
            :class:`NodeStatus` values. If present (default is ``None``),
            it is used  to manually define the status at any node on the grid.

        )doc");

    pgrid.def_static(
        "from_length",
        [](fs::py_profile_grid::size_type size,
           fs::py_profile_grid::length_type length,
           const profile_bstatus_type& bounds_status,
           const profile_nstatus_type& nodes_status)
        {
            auto bstatus = std::visit([](auto&& bs) { return fs::profile_boundary_status(bs); },
                                      bounds_status);

            if (!nodes_status.has_value())
            {
                return fs::py_profile_grid::from_length(size, length, bstatus);
            }
            else
            {
                const auto& nstatus = nodes_status.value();
                return fs::py_profile_grid::from_length(size, length, bstatus, nstatus);
            }
        },
        py::arg("size"),
        py::arg("length"),
        py::arg("bounds_status"),
        py::arg("nodes_status") = py::none(),
        R"doc(from_length(self, size: int, length: float, bounds_status: NodeStatus | List[NodeStatus] | ProfileBoundaryStatus, nodes_status: Dict[int, NodeStatus] | None) -> None

        Profile grid initializer from a given total length.

        Parameters
        ----------
        size : int
            Total number of grid nodes.
        length : float
            Total physical length of the grid.
        bounds_status : :class:`NodeStatus` or list or :class:`ProfileBoundaryStatus`
            Status at boundary nodes (left/right grid edges).
        nodes_status : dict, optional
            A dictionary where keys are node indices and values are
            :class:`NodeStatus` values. If present (default is ``None``),
            it is used  to manually define the status at any node on the grid.

        )doc");

    fs::register_grid_static_properties(pgrid);
    fs::register_base_grid_properties(pgrid);
    fs::register_structured_grid_properties(pgrid);
    fs::register_grid_methods(pgrid);

    /*
    ** Raster Grid
    */

    // ==== Binding of the raster_node structure ==== //
    py::class_<fs::raster_node> raster_node(m, "RasterNode", "Represents a raster grid node.");

    raster_node.def(py::init<std::size_t, std::size_t, fs::node_status>(),
                    py::arg("row"),
                    py::arg("col"),
                    py::arg("status"),
                    R"doc(__init__(self, row: int, col: int, status: NodeStatus) -> None

        Raster node initializer (internal use only).

        Parameters
        ----------
        row : int
            Node row index.
        col : int
            Node column index.
        status : :class:`NodeStatus`
            Node status.

        )doc");

    raster_node.def_readwrite("row", &fs::raster_node::row, "Node row index.")
        .def_readwrite("col", &fs::raster_node::col, "Node column index.")
        .def_readwrite("status", &fs::raster_node::status, "Node status.");

    // ==== Binding of the raster_neighbor structure ==== //
    py::class_<fs::raster_neighbor> raster_neighbor(
        m, "RasterNeighbor", "Represents a raster grid node neighbor.");

    raster_neighbor.def(
        py::init<std::size_t, std::size_t, std::size_t, double, fs::node_status>(),
        py::arg("flatten_idx"),
        py::arg("row"),
        py::arg("col"),
        py::arg("distance"),
        py::arg("status"),
        R"doc(__init__(self, flatten_idx: int, row: int, col: int, distance: float, status: NodeStatus) -> None

        Neighbor initializer (internal use only).

        Parameters
        ----------
        flatten_idx : int
            Neighbor flat index.
        row : int
            Neighbor row index.
        col : int
            Neighbor column index.
        distance : float
            Neighbor distance.
        status : :class:`NodeStatus`
            Neighbor status.

        )doc");

    raster_neighbor.def("__eq__", &fs::raster_neighbor::operator==)
        .def_readwrite("flatten_idx", &fs::raster_neighbor::flatten_idx, "Neighbor flatten index.")
        .def_readwrite("row", &fs::raster_neighbor::row, "Neighbor row index.")
        .def_readwrite("col", &fs::raster_neighbor::col, "Neighbor column index.")
        .def_readwrite("distance", &fs::raster_neighbor::distance, "Neighbor distance.")
        .def_readwrite("status", &fs::raster_neighbor::status, "Neighbor status.");

    // ==== Binding of the raster_boundary_status structure ==== //
    py::class_<fs::raster_boundary_status> raster_bstatus(m,
                                                          "RasterBoundaryStatus",
                                                          R"doc(Status at grid boundary nodes.

        To disambiguate the cases at each of the four raster grid corners, node
        status is set from one of their two overlaping borders according to the
        following precedance order:

        fixed value > fixed gradient > looped > core

        )doc");

    raster_bstatus.def(py::init<const fs::node_status>(),
                       py::arg("status"),
                       R"doc(__init__(*args, **kwargs) -> None

        Overloaded initializer.

        1. ``__init__(self, status: NodeStatus) -> None``

        Set the same status for all the boundary nodes.

        status : :class:`NodeStatus`
            Boundary node status.

        )doc");
    raster_bstatus.def(py::init<const std::array<fs::node_status, 4>&>(),
                       py::arg("status"),
                       R"doc(

        2. ``__init__(self, status: List[NodeStatus]) -> None``

        Set status at the left, right, top and bottom border nodes.

        status : list of :class:`NodeStatus`
            4-length list of node status for each border.

        )doc");

    raster_bstatus.def_readwrite("left", &fs::raster_boundary_status::left, "Left boundary status.")
        .def_readwrite("right", &fs::raster_boundary_status::right, "Right boundary status.")
        .def_readwrite("bottom", &fs::raster_boundary_status::bottom, "Bottom boundary status.")
        .def_readwrite("top", &fs::raster_boundary_status::top, "Top boundary status.")
        .def_property_readonly("is_horizontal_looped",
                               &fs::raster_boundary_status::is_horizontal_looped,
                               "Horizontal periodicity.")
        .def_property_readonly("is_vertical_looped",
                               &fs::raster_boundary_status::is_vertical_looped,
                               "Vertical periodicity.");

    // ==== Binding of the RasterGrid class ==== //
    py::class_<fs::py_raster_grid> rgrid(m,
                                         "RasterGrid",
                                         R"doc(A 2-dimensional raster grid.

        This raster grid type has 8-node connectivity (i.e., including diagonals).

        )doc");


    using raster_bstatus_type
        = std::variant<fs::node_status, std::array<fs::node_status, 4>, fs::raster_boundary_status>;
    using raster_nstatus_type = std::optional<fs::py_raster_grid::nodes_status_map_type>;

    rgrid.def(
        py::init(
            [](const fs::py_raster_grid::shape_type& shape,
               const xt::pytensor<double, 1>& spacing,
               const raster_bstatus_type& bounds_status,
               const raster_nstatus_type& nodes_status)
            {
                auto bstatus = std::visit([](auto&& bs) { return fs::raster_boundary_status(bs); },
                                          bounds_status);

                if (!nodes_status.has_value())
                {
                    return fs::py_raster_grid(shape, spacing, bstatus);
                }
                else
                {
                    const auto& nstatus = nodes_status.value();
                    return fs::py_raster_grid(shape, spacing, bstatus, nstatus);
                }
            }),
        py::arg("shape"),
        py::arg("spacing"),
        py::arg("bounds_status"),
        py::arg("nodes_status") = py::none(),
        R"doc(__init__(self, size: List[int], spacing: array_like, bounds_status: NodeStatus | List[NodeStatus] | RasterBoundaryStatus, nodes_status: Dict[Tuple[int, int], NodeStatus] | None = None) -> None

        Raster grid initializer.

        Parameters
        ----------
        shape : tuple
            Shape of the grid (number of rows and cols).
        spacing : array_like
            Distance between two adjacent grid nodes (row, cols).
        bounds_status : :class:`NodeStatus` or list or :class:`RasterBoundaryStatus`
            Status at boundary borders (left/right/top/bottom grid borders).
        nodes_status : dict, optional
            A dictionary where keys are node (row, col) indices and values are
            :class:`NodeStatus` values. If present (default is ``None``),
            it is used to manually define the status at any node on the grid.

        )doc");

    rgrid.def_static(
        "from_length",
        [](const fs::py_raster_grid::shape_type& shape,
           const xt::pytensor<double, 1>& length,
           const raster_bstatus_type& bounds_status,
           const raster_nstatus_type& nodes_status)
        {
            auto bstatus = std::visit([](auto&& bs) { return fs::raster_boundary_status(bs); },
                                      bounds_status);

            if (!nodes_status.has_value())
            {
                return fs::py_raster_grid::from_length(shape, length, bstatus);
            }
            else
            {
                const auto& nstatus = nodes_status.value();
                return fs::py_raster_grid::from_length(shape, length, bstatus, nstatus);
            }
        },
        py::arg("shape"),
        py::arg("length"),
        py::arg("bounds_status"),
        py::arg("nodes_status") = py::none(),
        R"doc(from_length(self, shape: tuple, length: array_like, bounds_status: NodeStatus | List[NodeStatus] | RasterBoundaryStatus, nodes_status: Dict[Tuple[int, int], NodeStatus] | None = None) -> None

        Raster grid initializer from given total lengths.

        Parameters
        ----------
        shape : tuple
            Shape of the grid (number of rows and cols).
        length : array_like
            Total physical length of the grid in y (rows) and x (cols).
        bounds_status : :class:`NodeStatus` or list or :class:`RasterBoundaryStatus`
            Status at boundary borders (left/right/top/bottom grid borders).
        nodes_status : dict, optional
            A dictionary where keys are node (row, col) indices and values are
            :class:`NodeStatus` values. If present (default is ``None``),
            it is used to manually define the status at any node on the grid.

        )doc");

    fs::register_grid_static_properties(rgrid);
    fs::register_base_grid_properties(rgrid);
    fs::register_structured_grid_properties(rgrid);
    fs::register_grid_methods(rgrid);

    auto raster_grid_funcs = fs::py_raster_grid_funcs();

    rgrid.def(
        "neighbors_indices", raster_grid_funcs.m_neighbors_indices, py::arg("row"), py::arg("col"));

    rgrid.def("neighbors", raster_grid_funcs.m_neighbors, py::arg("row"), py::arg("col"));
}
