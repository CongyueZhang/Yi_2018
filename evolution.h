#pragma once

struct mesh;
class half_edge_connectivity;
struct reduction_metric;

mesh evolution_stream(mesh& obj, const size_t& target_num_vertices, const size_t& Np, const size_t& n_iter, const double& F, const double& Cr, const double& tao, const double& nc, reduction_metric& metric);
