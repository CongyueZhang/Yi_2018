#pragma once

#include "half_edges.h"
#include <unordered_set>

inline bool ci_equal(const char* a, const char* b)
{
	if (!a || !b) return false;
	while (*a && *b && std::tolower(*a) == std::tolower(*b))
		a++, b++;
	return !*a && !*b;
}

struct k_ring
{
	template<typename F>
	k_ring(const half_edge_connectivity& connec, unsigned int k, const char* option, uint32_t v, F f) :
		connectivity(connec)
	{
		if (ci_equal(option, "edge"))
			traverse_k_ring_edge(k, v, f);
		else if (ci_equal(option, "halfedge"))
			traverse_k_ring(k, v, f);
		REM = NULL;
		REM_NLD = NULL;
	}

	k_ring(const half_edge_connectivity& connec, const char* option, std::unordered_set<uint32_t>& REM2update, std::unordered_set<uint32_t>& REM_NLD2update, uint32_t v) :
		connectivity(connec)
	{
		REM = &REM2update;
		REM_NLD = &REM_NLD2update;

		if (ci_equal(option, "edge_ring"))
			store_edge_ring(v);
		else if (ci_equal(option, "vertex_ring"))
			store_2_ring(v);
	}

	template<typename F>
	k_ring& extend(unsigned int k, uint32_t v, F f)				// 通过调用extend来遍历一个edge的两个端点
	{
		visited_vertices.clear();
		traverse_k_ring(k, v, f);
		return *this;
	}

private:
	void store_2_ring(uint32_t v)
	{
		traverse_1_ring(v, [&](const half_edge& h)
			{
				if (visited_edges.count(h.index) == 0)
				{
					REM_NLD->insert(h.index);
					visited_edges.insert(h.index);
					visited_edges.insert(h.opposite().index);
				}
				visited_vertices.insert(h.vertex());
			});

		visited_vertices.erase(v);

		for (uint32_t v2visite : visited_vertices)
		{
			traverse_1_ring(v2visite, [&](const half_edge& h)
				{
					if (visited_edges.count(h.index) == 0 && h.opposite().vertex() == v2visite)				// 外圈需要更新的是射出的halfedge 
						REM->insert(h.index);
				});
		}
	}

	void store_edge_ring(uint32_t h)
	{
		half_edge he = connectivity.handle(h);

		REM_NLD->insert(he.index);
		visited_vertices.insert(he.vertex());
		REM_NLD->insert(he.next().index);
		visited_vertices.insert(he.next().vertex());
		REM_NLD->insert(he.next().next().index);
		visited_vertices.insert(he.next().next().vertex());
		REM_NLD->insert(he.opposite().next().index);
		visited_vertices.insert(he.opposite().next().vertex());
		REM_NLD->insert(he.opposite().next().next().index);

		visited_edges = *REM_NLD;

		for (uint32_t v : visited_vertices)
			traverse_1_ring(v, [&](const half_edge& h)						// 下面的是传给tranverse_1_ring的函数
				{
					if (visited_edges.count(h.index) == 0 && visited_edges.count(h.opposite().index) == 0)
					{
						if (h.opposite().vertex() == v)								// 外圈需要更新的是射出的halfedge 
							REM->insert(h.index);
					}
				});
	}

	/// 遍历vertex包含v的所有 *halfedge*
	template<typename F>
	void traverse_k_ring(unsigned int k, uint32_t v, F f)
	{
		if (k == 0) return;
		visited_vertices.insert(v);
		traverse_1_ring(v, [&](const half_edge& h)				// 下面的是传给tranverse_1_ring的函数
			{
				if (visited_edges.count(h.index) == 0) f(h);			// 防止重复访问某个Halfedge
				visited_edges.insert(h.index);

				const uint32_t vert = h.vertex();
				if (k > 1 && visited_vertices.count(vert) == 0) traverse_k_ring(k - 1, vert, f);		// 同时迭代v的每一个对侧点
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
		while (h.index != start && h.is_valid())
		{
			if (start == uint32_t(-1)) start = h.index;
			f(h);
			f(h = h.next());
			f(h = h.next());
			h = h.opposite();
		}
	}

private:
	std::unordered_set<uint32_t> visited_edges;
	std::unordered_set<uint32_t> visited_vertices;

	std::unordered_set<uint32_t>* REM;
	std::unordered_set<uint32_t>* REM_NLD;

	const half_edge_connectivity& connectivity;
};