#pragma once

#include "half_edges.h"
#include <unordered_set>

struct k_ring
{
	template<typename F>
	k_ring(const half_edge_connectivity& connec, unsigned int k, std::string option, uint32_t v, F f) :
		connectivity(connec)
	{
		if(option == "halfedge")
			traverse_k_ring(k, v, f);
		if (option == "edge")
			traverse_k_ring_edge(k, v, f);
	}

	template<typename F>
	k_ring& extend(unsigned int k, uint32_t v, F f)				// 通过调用extend来遍历一个edge的两个端点
	{
		visited_vertices.clear();
		traverse_k_ring(k, v, f);
		return *this;
	}

private:
	/// 遍历vertex包含v的所有 *halfedge*
	template<typename F>
	void traverse_k_ring(unsigned int k, uint32_t v, F f)
	{
		if(k == 0) return;
		visited_vertices.insert(v);
		traverse_1_ring(v, [&](const half_edge& h)				// 下面的是传给tranverse_1_ring的函数
		{
			if(visited_edges.count(h.index) == 0) f(h);			// 防止重复访问某个Halfedge
			visited_edges.insert(h.index);

			const uint32_t vert = h.vertex();					
			if(k > 1 && visited_vertices.count(vert) == 0) traverse_k_ring(k - 1, vert, f);		// 同时迭代v的每一个对侧点
		});
	}

	/// 自己加的
	/// 遍历包含vertex v的所有 *edge*
	template<typename F>
	void traverse_k_ring_edge(unsigned int k, uint32_t v, F f)
	{
		if (k == 0) return;
		visited_vertices.insert(v);
		traverse_1_ring(v, [&](const half_edge& h)						// 下面的是传给tranverse_1_ring的函数
			{
				if (visited_edges.count(h.index) == 0) f(h);			// 防止重复访问某个edge
				visited_edges.insert(h.index);
				visited_edges.insert(h.opposite().index);

				const uint32_t vert = h.vertex();
				if (k > 1 && visited_vertices.count(vert) == 0) traverse_k_ring(k - 1, vert, f);		// 同时迭代v的每一个对侧点
			});
	}

	template<typename F>
	void traverse_1_ring(uint32_t v, F f)										// 遍历vertex中包含v的三角形的所有halfedge，并对每个halfedge执行 f 操作
	{
		half_edge h = connectivity.handle(connectivity.vertex_half_edge(v));
		uint32_t start = uint32_t(-1);
		while(h.index != start && h.is_valid())
		{
			if(start == uint32_t(-1)) start = h.index;
			f(h);
			f(h = h.next());
			f(h = h.next());
			h = h.opposite();
		}
	}

private:
	std::unordered_set<uint32_t> visited_edges;
	std::unordered_set<uint32_t> visited_vertices;
	const half_edge_connectivity& connectivity;
};
