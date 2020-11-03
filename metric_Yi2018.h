#pragma once

#include "metric.h"

struct metric_Yi2018 : public reduction_metric
{
	const mesh* obj;
	half_edge_connectivity* connectivity;
	bool lock_boundaries = false;

	void setup(const mesh& obj, half_edge_connectivity& connec) override;
	bool flip_valid(uint32_t half_edge) const override;
	double cost_split(uint32_t half_edge) const override;
	double cost_collapse(uint32_t half_edge ,uint32_t to_keep ,uint32_t to_remove) const override;
	bool remove_valid(uint32_t half_edge) override;
	bool delaunay_valid(uint32_t h_index) const override;
	std::pair<double, Eigen::Vector3d> split_position(const half_edge he) const override;

	unsigned int post_collapse(uint32_t, uint32_t, uint32_t)
	{
		return 1;
	}
};
