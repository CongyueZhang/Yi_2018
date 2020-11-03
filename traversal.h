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
	k_ring& extend(unsigned int k, uint32_t v, F f)				// ͨ������extend������һ��edge�������˵�
	{
		visited_vertices.clear();
		traverse_k_ring(k, v, f);
		return *this;
	}

private:
	/// ����vertex����v������ *halfedge*
	template<typename F>
	void traverse_k_ring(unsigned int k, uint32_t v, F f)
	{
		if(k == 0) return;
		visited_vertices.insert(v);
		traverse_1_ring(v, [&](const half_edge& h)				// ������Ǵ���tranverse_1_ring�ĺ���
		{
			if(visited_edges.count(h.index) == 0) f(h);			// ��ֹ�ظ�����ĳ��Halfedge
			visited_edges.insert(h.index);

			const uint32_t vert = h.vertex();					
			if(k > 1 && visited_vertices.count(vert) == 0) traverse_k_ring(k - 1, vert, f);		// ͬʱ����v��ÿһ���Բ��
		});
	}

	/// �Լ��ӵ�
	/// ��������vertex v������ *edge*
	template<typename F>
	void traverse_k_ring_edge(unsigned int k, uint32_t v, F f)
	{
		if (k == 0) return;
		visited_vertices.insert(v);
		traverse_1_ring(v, [&](const half_edge& h)						// ������Ǵ���tranverse_1_ring�ĺ���
			{
				if (visited_edges.count(h.index) == 0) f(h);			// ��ֹ�ظ�����ĳ��edge
				visited_edges.insert(h.index);
				visited_edges.insert(h.opposite().index);

				const uint32_t vert = h.vertex();
				if (k > 1 && visited_vertices.count(vert) == 0) traverse_k_ring(k - 1, vert, f);		// ͬʱ����v��ÿһ���Բ��
			});
	}

	template<typename F>
	void traverse_1_ring(uint32_t v, F f)										// ����vertex�а���v�������ε�����halfedge������ÿ��halfedgeִ�� f ����
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
