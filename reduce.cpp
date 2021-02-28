#include "reduce.h"
#include "mesh.h"
#include "ranges.h"
#include "geometry.h"
#include "half_edges.h"
#include "traversal.h"

#include <igl/opengl/glfw/Viewer.h>

#include <execution>
#include<math.h>
#include <iostream>
#include <time.h>
#include <string.h>

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

/// 遍历halfedge h所在的两个triangle的另外四个edge
/// --------------------------------------------
template<typename F>
void reduction::traverse_2_triangles(half_edge h, F f)						
{
	half_edge g = h.opposite();
	f(h = h.next());
	f(h = h.next());
	f(g = g.next());
	f(g = g.next());
}

void reduction::Liu_process_he(const uint32_t& hx)
{
	half_edge he = connectivity.handle(hx);
	bool he_delaunay_valid = metric.delaunay_valid(hx);

	if (!he_delaunay_valid && metric.flip_valid(hx))
	{
		connectivity.flip_edge(hx);
		candidatesNLD._delete(he.index);
		candidatesNLD._delete(he.opposite().index);
		// 看对其他边有无影响  （hx所在的两个三角形中的其他边）
		std::unordered_set<uint32_t> he2update;
		traverse_2_triangles(he, [&](half_edge& h) {
			he2update.insert(h.edge_index());
			});
		for (uint32_t h : he2update)		
			Liu_process_he(h);
	}
	else 
		add_split(he, he_delaunay_valid);			// 必须在else里，因为flip中he_delaunay_valid可能会改变，不能还用以前的值
		// 如果flip，该边最终肯定会在flip中的迭代里更新
}

/// Liu方法中的split函数
/// -------------------------------------------------------------------
void reduction::Liu_perform_split(const detail::candidate_operation& c)
{
	half_edge h = connectivity.handle(c.index);

	auto [nonDelaunay, vertexPosition] = metric.split_position(h);
	obj.vertices.push_back(vertexPosition);											/// 将新插入的点的坐标，推入 obj.vertices 这个 vector 的最后

	uint32_t vertex_index = obj.vertices.size() - 1;
	connectivity.split_edge(h.index, vertex_index);				
	stats.on_operation(Split);

	
	/// debug
	/// -----------------------------------------
		// 临时debug加的
	// ------------------------
	mesh Mesh1;
	Mesh1.vertices = obj.vertices;
	connectivity.on_triangles([&](const std::array<uint32_t, 3>& t) { Mesh1.triangles.push_back(t); });			
	remove_standalone_vertices(Mesh1, connectivity);
	Mesh1.save("delaunay_test_split" + std::to_string(stats.num_split) + "_NLD" + std::to_string(candidatesNLD.size()) + ".obj");
	// ------------------------
	/*
	// Plot the mesh
	mesh Mesh1;
	Mesh1.vertices = obj.vertices;
	connectivity.on_triangles([&](const std::array<uint32_t, 3>& t) { Mesh1.triangles.push_back(t); });		
	remove_standalone_vertices(Mesh1, connectivity);
	igl::opengl::glfw::Viewer viewer;
	viewer.data().set_mesh(Mesh1.matrix_vertices(), Mesh1.matrix_triangles());
	viewer.core().set_rotation_type(igl::opengl::ViewerCore::ROTATION_TYPE_TRACKBALL);
	viewer.launch();
	/// -----------------------------------------
	*/

	/// update neighborhoods
	std::unordered_set<uint32_t> he2update;
	const char* option = "edge";
	k_ring(connectivity, 1, option, vertex_index, [&](half_edge he) {
		he2update.insert(he.edge_index());
		});

	for (uint32_t h : he2update)					// 不能在traverse_k_ring_edge中直接更新，因为flip会改变结构，导致回不到start edge
	{
		Liu_process_he(h);
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

	metric.setup(obj, connectivity);

	candidatesNLD.clear();
	candidatesREM.clear();
	// Fill the candidate queue
	candidatesREM.reserve(connectivity.num_half_edges(), connectivity.num_half_edges());
	candidatesNLD.reserve(std::ceil(connectivity.num_half_edges()), std::ceil(connectivity.num_half_edges()));

	/// 构建NLD queue
	/// --------------------------------------------
	size_t _number = connectivity.num_half_edges();
	for (uint32_t hx = 0; hx < _number; hx++)
	{
		const half_edge he = connectivity.handle(hx);
		const uint32_t ho = he.opposite().index;
		if (!he.is_boundary() && hx > ho) continue;						/// he.index > ho 是为了防止重复遍历一条edge

		Liu_process_he(hx);
	}

	/// 进行split操作
	/// 直到NLD queue为空
	/// --------------------------
	while (!candidatesNLD.empty())
	{
		const detail::candidate_operation c = candidatesNLD.pop();
		Liu_perform_split(c);
	}

	nLiu = obj.num_vertices();
	sequence_d = 2 * (nLiu - ini_num_vertices);

	/// 临时debug加的
	mesh Mesh1;
	Mesh1.vertices = obj.vertices;
	Mesh1.triangles.clear();
	connectivity.on_triangles([&](const std::array<uint32_t, 3>& t) { Mesh1.triangles.push_back(t); });
	remove_standalone_vertices(Mesh1, connectivity);

	Mesh1.save("initialization_test.obj");

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

void reduction::add_split(const half_edge& he, const bool& he_delaunay_valid)
{
	if (he.is_boundary()) return;										
	if (he_delaunay_valid)													/// 如果本身是delaunay的，则不需要split
	{
		candidatesNLD._delete(he.index);									/// 该边之前可能是nonDelaunay的，现在变成了delaunay，所以要把NLD中的删除
		candidatesNLD._delete(he.opposite().index);
		return;
	}

	detail::candidate_operation c;
	c.index = he.edge_index();
	c.weight = metric.cost_split(he.index);
	candidatesNLD.push(c);
}

void reduction::perform_flip(const half_edge& he)							/// 更新flip后的边及其neighbour
{
	connectivity.flip_edge(he.index);

	candidatesNLD._delete(he.index);
	candidatesNLD._delete(he.opposite().index);

	stats.on_operation(Flip);								/// TODO: 考虑对临边delaunay情况的影响

	// he所在的两个三角形所有vertex的one ring都要Update一下（对于REM）

	// TODO: 修复bug
	// 如果在traverse的过程中flip了，则traverse会死循环

	std::unordered_set<uint32_t> REM2update;
	std::unordered_set<uint32_t> REM_NLD2update;
	const char* option = "edge_ring";
	k_ring(connectivity, option, REM2update, REM_NLD2update, he.index);

	for (uint32_t hx : REM_NLD2update)					// 不能在traverse_k_ring_edge中直接更新，因为flip会改变结构，导致回不到start edge
		process_he(connectivity.handle(hx));

	for (uint32_t hx : REM2update)
	{
		half_edge h = connectivity.handle(hx);
		if(h.is_valid())
			add_one_collapse(h);
	}
}

void reduction::process_he(const half_edge& he)
{
	if (he.is_valid())
	{
		bool he_delaunay_valid = metric.delaunay_valid(he.index);
		if (!he_delaunay_valid && metric.flip_valid(he.index))										/// 如果不是delaunay_valid，能flip先flip，不能进行collapse/split
			perform_flip(he);
		else
		{
			add_split(he, he_delaunay_valid);
		}
		add_collapse(he);
	}
}

void reduction::perform_collapse(const detail::candidate_operation& c)			/// 从pirority queue中移除掉删掉了的halfedges，并更新周围的neighbour
{
	half_edge he = connectivity.handle(c.index);								/// he: to_remove -> to_keep
	uint32_t to_keep = he.vertex();
	
	/// 将需要更新的halfedge的坐标存起来
	/// 需要更新的是以 to_remove 的 2 ring
	/// TO OPTIMIZE: REM需要更新2 ring，NLD只需要更新1 ring
	uint32_t to_remove = he.opposite().vertex();
	std::unordered_set<uint32_t> REM2update;
	std::unordered_set<uint32_t> REM_NLD2update;
	const char* option = "vertex_ring";
	k_ring(connectivity, option, REM2update, REM_NLD2update, to_remove);

	stats.on_operation(Collapse);

	/*
	* debug测试
	if (stats.num_total == 4968)
	{
		half_edge he_opp = he.opposite();
		half_edge he_opp_next = he.opposite().next();
		half_edge he_opp_next_opp = he.opposite().next().opposite();
		half_edge he_opp_next_opp_next = he.opposite().next().opposite().next();
		half_edge he_opp_next_opp_next_next = he.opposite().next().opposite().next().next();
		half_edge he_opp_next_opp_next_next_opposite = he.opposite().next().opposite().next().next().opposite();
		half_edge he_opp_next_opp_next_opp = he.opposite().next().opposite().next().opposite();
		half_edge he_opp_next_opp_next_opp_next = he.opposite().next().opposite().next().opposite().next();
		half_edge he_opp_next_opp_next_opp_next_next = he.opposite().next().opposite().next().opposite().next().next();
		half_edge he_opp_next_opp_next_opp_next_next_opp = he.opposite().next().opposite().next().opposite().next().next().opposite();
		half_edge he_opp_next_opp_next_opp_next_next_opp_next = he.opposite().next().opposite().next().opposite().next().next().opposite().next();
		half_edge he_opp_next_opp_next_opp_next_next_opp_next_next = he.opposite().next().opposite().next().opposite().next().next().opposite().next().next();
		int temp = 1;
	}	
	*/

	for (uint32_t h : connectivity.collapse_edge(c.index))				/// connectivity.collapse_edge()返回的是 index of removed half edges
	{
		candidatesREM._delete(h);										/// 从priority queue中移除这些删掉了的边
		candidatesNLD._delete(h);
	}
	
	/*
	if (stats.num_total >= 7650)
	{
		mesh Mesh1;
		Mesh1.vertices = obj.vertices;
		connectivity.on_triangles([&](const std::array<uint32_t, 3>& t) { Mesh1.triangles.push_back(t); });			
		remove_standalone_vertices(Mesh1, connectivity);
		Mesh1.save("delaunay_test_after_collapse" + std::to_string(stats.num_total) + ".obj");
		
		std::cout << "num_total = " << std::to_string(stats.num_total) << std::endl;
		system("pause");
	}		
	*/
	
	for (uint32_t h : REM_NLD2update)
	{
		half_edge he = connectivity.handle(h);
		if (he.is_valid())
		{
			add_collapse(he);
			add_split(he,metric.delaunay_valid(h));
		}
	}	

	for (uint32_t h : REM2update)
	{
		add_one_collapse(connectivity.handle(h));
	}
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

	/// update neighborhoods
	std::unordered_set<uint32_t> REM2update;
	std::unordered_set<uint32_t> REM_NLD2update;
	const char* option = "vertex_ring";
	k_ring(connectivity, option, REM2update, REM_NLD2update, vertex_index);

	for (uint32_t h : REM_NLD2update)					/// 不能在traverse_k_ring_edge中直接更新，因为flip会改变结构，导致回不到start edge
		process_he(connectivity.handle(h));

	for (uint32_t h : REM2update)
		add_one_collapse(connectivity.handle(h));
}

std::pair<mesh, std::vector<size_t>> reduction::reduce_stream(Eigen::ArrayXf X)
{
	stats.clear();
	stats.num_vertices = ini_num_vertices;

	/// delete connectivity;							
	obj.vertices.resize(ini_num_vertices);
	obj.vertices = ori_vertices;
	connectivity = obj.half_edges();

	metric.setup(obj, connectivity);

	candidatesNLD.clear();
	candidatesREM.clear();

	// Fill the candidate queue
	candidatesREM.reserve(connectivity.num_half_edges(), connectivity.num_half_edges());
	candidatesNLD.reserve(std::ceil(connectivity.num_half_edges() / 2.0), std::ceil(connectivity.num_half_edges() / 2.0));

	for (uint32_t hx = 0; hx < connectivity.num_half_edges(); hx++) /// TODO: parallel
	{
		half_edge he = connectivity.handle(hx);
		if (!he.is_boundary() && he.index > he.opposite().index) continue;

		process_he(he);																
	}

	/// Algorithm2
	/// 问题 S返回后好像没用到？？
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
			// 临时debug加的
			// for (; j < 6599; ++j)
			for (; j < X[i]; ++j)
			{
				if (candidatesNLD.empty())
				{
					_skip = true;
					break;
				}
				if (stats.num_split >= nLiu)
				{
					// set S(x) = Sdummy and return 空集
					S.resize(1);
					S[0] = -1;
					return { Mesh, S };			// dummy
				}

				const detail::candidate_operation c = candidatesNLD.pop();
				perform_split(c);
			}
		}

		if (i % 2 != 0)							/// even, E = Ec
		{
			// 临时debug改的
		    // for (; j < 500; ++j)
			for (; j < X[i]; ++j)
			{
				if (candidatesNLD.empty())
				{
					_skip = true;
					break;
				}
				if (stats.num_vertices == target_num_vertices)
				{
					_skip = true;
					break;
				}
				if (candidatesREM.empty())
				{
					_skip = true;
					// set S(x) = Sdummy and return 空集
					S.resize(1);
					S[0] = -1;
					return { Mesh, S };			// dummy
				}

				const detail::candidate_operation c = candidatesREM.pop();
				perform_collapse(c);
			}
		}

		if (_skip)
			S.push_back(j-1);
		else
			S.push_back(j);		/// 因为j本来是index，比实际次数少1

		if (candidatesNLD.empty())
		{
			size_t num_collap = 0;
			while (!candidatesREM.empty() && stats.num_vertices > target_num_vertices)			// TODO：这里好像有问题，好像是跳过了==？
			{
				++num_collap;
				const detail::candidate_operation c = candidatesREM.pop();
				perform_collapse(c);
			}
			
			if (stats.num_vertices != target_num_vertices)
			{
				S.resize(1);
				S[0] = -1;
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
			Mesh.triangles.clear();
			connectivity.on_triangles([&](const std::array<uint32_t, 3>& t) { Mesh.triangles.push_back(t); });			// 在reduce后更新了mesh的triangle
			remove_standalone_vertices(Mesh, connectivity);

			/// 临时debug加的
			mesh Mesh1;
			Mesh1.vertices = obj.vertices;
			Mesh1.triangles.clear();
			connectivity.on_triangles([&](const std::array<uint32_t, 3>& t) { Mesh1.triangles.push_back(t); });		
			remove_standalone_vertices(Mesh1, connectivity);

			Mesh1.save("test" + std::to_string(clock()) + ".obj");

			return { Mesh, S };
		}
	}
	S.resize(1);
	S[0] = -1;
	return { Mesh, S };
}