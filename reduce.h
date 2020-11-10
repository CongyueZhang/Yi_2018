#pragma once

#include <stdint.h>
#include <unordered_set>
#include <Eigen/Sparse>
#include "min_priority_queue.h"
#include "metric.h"


struct mesh;
class half_edge_connectivity;
struct reduction_metric;

enum _Operation { Split, Collapse, Flip };

namespace detail
{
	struct candidate_operation
	{
		uint32_t index = uint32_t(-1);								//halfedge在half_edges vector中的index，每对halfedge好像只存了index小的那个 （本质上是edge的Index？
		double weight = std::numeric_limits<double>::infinity();

		candidate_operation() = default;
		candidate_operation(uint32_t h, double p) : weight(p), index(h) {}

		bool operator<(const candidate_operation& c) const { return weight < c.weight; }			//按target.weight的大小排
		void make_cost_infinite() { weight = std::numeric_limits<double>::infinity(); }
	};

#if 0			
	struct statistics
	{
		statistics(const size_t&) {}
		void on_operation(_Operation _lastOperation) {}
		void on_invalid() {}
		void on_outdated() {}
		void on_reduction_end(size_t) {}
	};
#else
	struct statistics
	{
		size_t num_total = 0;
		size_t num_collapse = 0, num_split = 0, num_flip = 0;
		size_t num_vertices = 0;

		// 更新了策略，num_nonDelaunay的大小其实就是NLD中elements的数量
		// size_t num_nonDelaunay = 0;			// non-delaunay的 *edge*的数量

		void on_operation(_Operation _lastOperation)
		{
			++num_total;
			switch (_lastOperation)
			{
			case Flip: ++num_flip; break;
			case Collapse: ++num_collapse; --num_vertices; break;
			case Split: ++num_vertices; ++num_split; break;
			}
			//print_progress<0>("(reduce) " + std::to_string(current_size));  
		}

		/*
		void on_reduction_end(size_t remaining)
		{

			printf("Edges statistics:\n"
				"\tprocessed    : %zu\n"
				"\t  + invalid  : %zu (%.2f%%)\n"
				"\t    + flips  : %zu (%.2f%%)\n"
				"\t  + outdated : %zu (%.2f%%)\n"
				"\t    + flips  : %zu (%.2f%%)\n"
				"\tunprocessed  : %zu\n",
				num_total,
				remaining);
		}		
		*/


		void clear()
		{
			num_total = 0;
			num_collapse = 0, num_split = 0, num_flip = 0;
			num_vertices = 0;
		}
	};
#endif;
}

void remove_standalone_vertices(mesh& obj, const half_edge_connectivity& h);

struct reduction {
	mesh obj;
	half_edge_connectivity connectivity;
	reduction_metric& metric;
	detail::statistics stats;

	size_t ini_num_vertices = 0;
	size_t target_num_vertices = 0;

	reduction(mesh &obj, reduction_metric &metric, const size_t &target_vertices):obj(obj), metric(metric), target_num_vertices(target_vertices)
	{
		connectivity = obj.half_edges();
		ini_num_vertices = obj.num_vertices();
		metric.setup((const mesh&)obj, (half_edge_connectivity&)connectivity);
	}

	void initialize();
	size_t get_nLiu();
	size_t get_sequence_d();
	std::pair<mesh, std::vector<size_t>> reduce_stream(Eigen::ArrayXf X);

private:
	size_t nLiu = 0;
	size_t sequence_d = 0;		// operation sequence的长度
	min_priority_queue<detail::candidate_operation, std::function<uint32_t(const detail::candidate_operation& c)>> candidatesNLD;	//T是candidate_operation，IndexFn是uint32_t；将candidate_index函数给index_of
	min_priority_queue<detail::candidate_operation, std::function<uint32_t(const detail::candidate_operation& c)>> candidatesREM;
	std::unordered_set<uint32_t> visited_edges;
	std::unordered_set<uint32_t> visited_vertices;
	void add_collapse(const half_edge& he);
	void add_one_collapse(const half_edge& he);
	void add_split(const half_edge& he);
	void add_operation(const half_edge& he);
	void perform_flip(const half_edge& he);
	void process_he(const half_edge& he);
	void perform_collapse(const detail::candidate_operation& c);
	void perform_split(const detail::candidate_operation& c);

	template<typename F>
	void traverse_k_ring_edge(unsigned int k, uint32_t v, F f);
	template<typename F>
	void traverse_1_ring(uint32_t v, F f);
};

inline size_t reduction::get_nLiu() { return nLiu; };
inline size_t reduction::get_sequence_d() { return sequence_d; };
//inline uint32_t reduction::candidate_index(const detail::candidate_operation& c) { return c.index; };	// 求index的lambda函数