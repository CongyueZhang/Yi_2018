#include "reduce.h"
#include "mesh.h"
#include "ranges.h"
#include "geometry.h"
#include "half_edges.h"

#include <execution>
#include<math.h>
#include <iostream>

void remove_standalone_vertices(mesh& obj, const half_edge_connectivity& connectivity)
{
	std::vector<uint32_t> map(obj.vertices.size(), uint32_t(-1));
	uint32_t u = 0;
	for (uint32_t v = 0; v < obj.vertices.size(); v++)
	{
		if (!connectivity.is_connected(v)) continue;
		map[v] = u;
		if (u != v) obj.vertices[u] = obj.vertices[v];
		u++;
	}
	obj.vertices.resize(u);
	for (size_t t = 0; t < obj.triangles.size(); t++)
		for (uint32_t& idx : obj.triangles[t].idx)
			idx = map[idx];
}

/// 从traversal.h里copy来的
template<typename F>
void reduction::traverse_k_ring_edge(unsigned int k, uint32_t v, F f)
{
	visited_vertices.insert(v);
	traverse_1_ring(v, [&](const half_edge& h)				// 下面的是传给tranverse_1_ring的函数
		{
			if (visited_edges.count(h.index) == 0) f(h);			// 防止重复访问某个Halfedge
			visited_edges.insert(h.index);
			visited_edges.insert(h.opposite().index);

			const uint32_t vert = h.vertex();
			if (k > 1 && visited_vertices.count(vert) == 0) traverse_k_ring_edge(k - 1, vert, f);		// 同时迭代v的每一个对侧点
		});
}

template<typename F>
void reduction::traverse_1_ring(uint32_t v, F f)										// 遍历vertex中包含v的三角形的所有halfedge，并对每个halfedge执行 f 操作
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

/// TODO: 考虑boundary的情况
void reduction::initialize()
{
	std::function<uint32_t(const detail::candidate_operation& c)> candidate_index = [](const detail::candidate_operation& c)->uint32_t {
		return c.index;
	};
	candidatesNLD.setup_indexof(candidate_index);
	candidatesREM.setup_indexof(candidate_index);

	double l_max = 0, l_min = std::numeric_limits<double>::infinity();
	double angle_min = std::numeric_limits<double>::infinity();

	/// 在进行操作前先把所有需要flip的边flip掉
	/// 统计non_delaunay_valid的【edge】的数量
	/// 寻找l_max, l_min, angle_min
	size_t _number = connectivity.num_half_edges();
	for (uint32_t hx = 0; hx < _number; hx++)	/// TODO: parallel
	{
		const half_edge he = connectivity.handle(hx);
		const uint32_t ho = he.opposite().index;
		if (!he.is_boundary() && he.index < ho) continue;						/// he.index < ho 是为了防止重复遍历一条edge

		//if (!metric.delaunay_valid(hx))
		//{
			//if (metric.flip_valid(hx))
			//	connectivity->flip_edge(hx);
			//else
			//connectivity->half_edges[hx].delaunay_valid = false;
			//++stats.num_nonDelaunay;
		//}
		Eigen::Vector3d v1, v2;													/// he的两个顶点
		v1 = obj.vertices[he.vertex()];
		v2 = obj.vertices[he.opposite().vertex()];
		double l_v1v2 = (v1 - v2).norm();
		if (l_v1v2 < l_min)	l_min = l_v1v2;
		if (l_v1v2 > l_max) l_max = l_v1v2;

		auto [angle1, angle2] = opp_angle(obj, connectivity, hx);
		if (angle1 < angle_min)	angle_min = angle1;
		if (angle2 < angle_min) angle_min = angle2;
	}

	nLiu = ceil(ini_num_vertices * l_max / (l_min * sin(angle_min) * sin(angle_min)));
	sequence_d = 2 * (nLiu - ini_num_vertices);									/// operation sequence的长度
}

/// 一次add一对 (因为每个edge变化后，其对应的两个halfedge的collapse情况都会变化)
void reduction::add_collapse(const half_edge& he)
{
	add_one_collapse(he);
	add_one_collapse(he.opposite());
}

void reduction::add_one_collapse(const half_edge& he)
{
	if (!metric.remove_valid(he.index))
	{
		candidatesREM._delete(he.index);
	}
	else
	{
		const uint32_t to_keep = he.vertex();
		const uint32_t to_remove = he.next().next().vertex();

		detail::candidate_operation c;
		c.index = he.index;
		c.weight = metric.cost_collapse(he.index, to_keep, to_remove);				/// 计算这个edge collapse的cost
		candidatesREM.push(c);
	}
}

void reduction::add_split(const half_edge& he)
{
	if (he.is_boundary()) return;
	if (metric.delaunay_valid(he.index))										/// 如果本身是delaunay的，则不需要split
	{
		candidatesNLD._delete(he.index);		/// 该边之前可能是nonDelaunay的，现在变成了delaunay，所以要把NLD中的删除
		return;
	}

	detail::candidate_operation c;
	c.index = he.index;
	c.weight = metric.cost_split(he.index);
	candidatesNLD.push(c);
}

/// 将该halfedge推入priority queue
void reduction::add_operation(const half_edge& he)
{
	add_collapse(he);
	add_split(he);
}

void reduction::perform_flip(const half_edge& he)								/// 更新flip后的边及其neighbour
{
	connectivity.flip_edge(he.index);

	stats.on_operation(Flip);								/// TODO: 考虑对临边delaunay情况的影响

	// he所在的两个三角形所有vertex的one ring都要Update一下（对于REM）
	visited_edges.clear();
	visited_vertices.clear();
	traverse_k_ring_edge(1, he.vertex(),										/// collapse和split现在都是针对edge而不是Halfedge，所以要防止重复遍历同一个edge
		[&](half_edge h) {
			process_he(h);
		});
	traverse_k_ring_edge(1, he.opposite().vertex(),										/// collapse和split现在都是针对edge而不是Halfedge，所以要防止重复遍历同一个edge
		[&](half_edge h) {
			process_he(h);
		});
	traverse_k_ring_edge(1, he.next().vertex(),										/// collapse和split现在都是针对edge而不是Halfedge，所以要防止重复遍历同一个edge
		[&](half_edge h) {
			process_he(h);
		});
	traverse_k_ring_edge(1, he.opposite().next().vertex(),										/// collapse和split现在都是针对edge而不是Halfedge，所以要防止重复遍历同一个edge
		[&](half_edge h) {
			process_he(h);
		});

	// 
	/*
	* 以下这些应该是检查NLD的
	process_he(he);
	process_he(he.next());
	process_he(he.next().next());
	process_he(he.opposite().next());
	process_he(he.opposite().next().next());	
	*/

}

void reduction::process_he(const half_edge& he)
{
	if (he.is_valid())
	{
		if (!metric.delaunay_valid(he.index) && metric.flip_valid(he.index))										/// 如果不是delaunay_valid，能flip先flip，不能进行collapse/split
			perform_flip(he);
		else 
			add_operation(he);
	}
}

void reduction::perform_collapse(const detail::candidate_operation& c)			/// 从pirority queue中移除掉删掉了的halfedges，并更新周围的neighbour
{
	half_edge he = connectivity.handle(c.index);						/// he是从要移除的点出射的
	uint32_t to_keep = he.vertex();

	std::vector<uint32_t> he2update;
	/// 将需要更新的halfedge的坐标存起来
	/// 需要更新的是以 to_remove 的 2 ring
	/// TO OPTIMIZE: REM需要更新2 ring，NLD只需要更新1 ring
	uint32_t to_remove = he.opposite().vertex();
	visited_edges.clear();
	visited_vertices.clear();
	traverse_k_ring_edge(2, to_remove,										/// collapse和split现在都是针对edge而不是Halfedge，所以要防止重复遍历同一个edge
		[&](half_edge h) {
				he2update.push_back(h.index);										
		});	

	stats.on_operation(Collapse);

	for (uint32_t h : connectivity.collapse_edge(c.index))				/// connectivity.collapse_edge()返回的是 index of removed half edges
	{
		candidatesREM._delete(h);										/// 从priority queue中移除这些删掉了的边
		candidatesNLD._delete(h);
	}

	for (uint32_t h : he2update)
	{
		if (connectivity.handle(h).is_valid())							/// 更新未删掉的Halfedge
			process_he(connectivity.handle(h));
	}	

/// 临时debug加的
/// ----------------------------------------
// Plot the mesh
	
	mesh Mesh;
	Mesh.vertices = obj.vertices;
	connectivity.on_triangles([&](const std::array<uint32_t, 3>& t) { Mesh.triangles.push_back(t); });			// 在reduce后更新了mesh的triangle？
	remove_standalone_vertices(Mesh, connectivity);

	Mesh.save("oneNLD.obj");	
	

}

/// 先从priority queue中删掉已有的Halfedges
/// 再将新加的halfedge推入priority queue
/// 再更新neighbourhood
void reduction::perform_split(const detail::candidate_operation& c)
{
	half_edge h = connectivity.handle(c.index);

	auto [nonDelaunay, vertexPosition] = metric.split_position(h);
	obj.vertices.push_back(vertexPosition);											/// 将新插入的点的坐标，推入 obj.vertices 这个 vector 的最后
	
	// 采用了统计nondelaunay的新方法后，应该用不到了
	// stats.num_nonDelaunay += nonDelaunay - 1;									/// -1是因为h不再是non-delaunay的了

	uint32_t vertex_index = obj.vertices.size() - 1;
	connectivity.split_edge(h.index, vertex_index);									/// vector 的 index 从0开始

	stats.on_operation(Split);

	/// update neightbourhood
	/// 问题: process_he是类里的函数，不能直接传，并且connectivity,stats等的只在这个类存在
	/// 采用的解决方法：把相关函数单拿出来放到这边的类里
	/// TO OPTIMIZE: REM需要更新2 ring，NLD只需要更新1 ring
	visited_edges.clear();
	visited_vertices.clear();
	std::vector<uint32_t> he2update;
	traverse_k_ring_edge(2, vertex_index, [&](half_edge he) {
		he2update.push_back(h.index);
		});

	for (uint32_t h : he2update)					// 不能再traverse_k_ring_edge中直接更新，因为flip会改变结构，导致回不到start edge
	{
		if (connectivity.handle(h).is_valid())						
			process_he(connectivity.handle(h));
	}

	mesh Mesh1;
	Mesh1.vertices = obj.vertices;
	connectivity.on_triangles([&](const std::array<uint32_t, 3>& t) { Mesh1.triangles.push_back(t); });			// 在reduce后更新了mesh的triangle？
	remove_standalone_vertices(Mesh1, connectivity);

	Mesh1.save("split_test.obj");
}

std::pair<mesh, std::vector<size_t>> reduction::reduce_stream(Eigen::ArrayXf X)
{
	stats.clear();
	stats.num_vertices = ini_num_vertices;

	/// delete connectivity;			///TODO: 确认一下是否需要主动释放空间
	connectivity = obj.half_edges();
	metric.setup(obj, connectivity);

	obj.vertices.resize(ini_num_vertices);

	candidatesNLD.clear();
	candidatesNLD.clear();

	// Fill the candidate queue
	candidatesREM.reserve(connectivity.num_half_edges(), connectivity.num_half_edges());
	candidatesNLD.reserve(std::ceil(connectivity.num_half_edges() / 2.0), std::ceil(connectivity.num_half_edges() / 2.0));

	for (uint32_t hx = 0; hx < connectivity.num_half_edges(); hx++) /// TODO: parallel
	{
		half_edge he = connectivity.handle(hx);
		if (!he.is_boundary() && he.index > he.opposite().index) continue;								/// he.index > ho 防止重复遍历一条edge

		process_he(he);																
	}

	/// Algorithm2
	/// 问题 S是做什么用的？返回后好像没用到？？
	/// ----------------------------------------
	std::vector<size_t> S;
	bool _skip = false;
	mesh Mesh;
	X = X.round();

	for (size_t i = 0; i <= X.size(); ++i)
	{
		_skip = false;
		size_t j = 0;
		if (i % 2 == 0)								/// odd, E = Es  (因为X坐标从0开始，所以i为even时，对应的是odd)
		{
			for (; j < X[i]; ++j)
			{
				// 问题： 做第一次operation前是否要判断？？
				// 解决方案：改成了先判断再操作

				if (candidatesNLD.empty())
				{
					_skip = true;
					break;
				}
				if (stats.num_split >= nLiu)
				{
					// set S(x) = Sdummy and return 空集
					S.resize(1, -1);
					return { Mesh, S };			// dummy
				}

				const detail::candidate_operation c = candidatesNLD.pop();
				perform_split(c);
				++stats.num_split;
				++stats.num_vertices;
			}
		}

		if (i % 2 != 0)							/// even, E = Ec
		{
			for (; j < X[i]; ++j)
			{
				if (candidatesNLD.empty())
				{
					_skip = true;
					break;
				}
				if (candidatesREM.empty())
				{
					_skip = true;
					// set S(x) = Sdummy and return 空集
					S.resize(1, -1);
					return { Mesh, S };			// dummy
				}

				const detail::candidate_operation c = candidatesREM.pop();
				perform_collapse(c);
				++stats.num_collapse;
				--stats.num_vertices;
			}
		}

		if (_skip)
			S.push_back(j);
		else
			S.push_back(j+1);		/// 因为j本来是index，比实际次数少1

		if (candidatesNLD.empty())
		{
			size_t num_collap = 0;
			while (!candidatesREM.empty() && stats.num_vertices != target_num_vertices)			// TODO：这里好像有问题，好像是跳过了==？
			{
				++num_collap;
				const detail::candidate_operation c = candidatesREM.pop();
				perform_collapse(c);
			}
			
			if (stats.num_vertices != target_num_vertices)
			{
				S.resize(1, -1);
				return { Mesh, S };				// dummy
			}

			if (S.size() % 2 == 0)
			{
				S.push_back(0);
				S.push_back(num_collap);
			}
			else
				S.push_back(num_collap);

			Mesh.vertices = obj.vertices;
			connectivity.on_triangles([&](const std::array<uint32_t, 3>& t) { Mesh.triangles.push_back(t); });			// 在reduce后更新了mesh的triangle？
			remove_standalone_vertices(Mesh, connectivity);
			obj.vertices.resize(ini_num_vertices);

			/// 临时debug加的
			mesh Mesh1;
			Mesh1.vertices = obj.vertices;
			connectivity.on_triangles([&](const std::array<uint32_t, 3>& t) { Mesh1.triangles.push_back(t); });			// 在reduce后更新了mesh的triangle？
			remove_standalone_vertices(Mesh1, connectivity);

			Mesh1.save("dinosaur2k_test.obj");

			return { Mesh, S };
		}
	}
	S.resize(1, -1);
	return { Mesh, S };
}