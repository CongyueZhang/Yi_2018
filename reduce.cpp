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

/// ����halfedge h���ڵ�����triangle�������ĸ�edge
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
		// ��������������Ӱ��  ��hx���ڵ������������е������ߣ�
		std::unordered_set<uint32_t> he2update;
		traverse_2_triangles(he, [&](half_edge& h) {
			he2update.insert(h.edge_index());
			});
		for (uint32_t h : he2update)		
			Liu_process_he(h);
	}
	else 
		add_split(he, he_delaunay_valid);			// ������else���Ϊflip��he_delaunay_valid���ܻ�ı䣬���ܻ�����ǰ��ֵ
		// ���flip���ñ����տ϶�����flip�еĵ��������
}

/// Liu�����е�split����
/// -------------------------------------------------------------------
void reduction::Liu_perform_split(const detail::candidate_operation& c)
{
	half_edge h = connectivity.handle(c.index);

	auto [nonDelaunay, vertexPosition] = metric.split_position(h);
	obj.vertices.push_back(vertexPosition);											/// ���²���ĵ�����꣬���� obj.vertices ��� vector �����

	uint32_t vertex_index = obj.vertices.size() - 1;
	connectivity.split_edge(h.index, vertex_index);				
	stats.on_operation(Split);

	
	/// debug
	/// -----------------------------------------
		// ��ʱdebug�ӵ�
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

	for (uint32_t h : he2update)					// ������traverse_k_ring_edge��ֱ�Ӹ��£���Ϊflip��ı�ṹ�����»ز���start edge
	{
		Liu_process_he(h);
	}
}

/// TODO: ����boundary�����
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

	/// ����NLD queue
	/// --------------------------------------------
	size_t _number = connectivity.num_half_edges();
	for (uint32_t hx = 0; hx < _number; hx++)
	{
		const half_edge he = connectivity.handle(hx);
		const uint32_t ho = he.opposite().index;
		if (!he.is_boundary() && hx > ho) continue;						/// he.index > ho ��Ϊ�˷�ֹ�ظ�����һ��edge

		Liu_process_he(hx);
	}

	/// ����split����
	/// ֱ��NLD queueΪ��
	/// --------------------------
	while (!candidatesNLD.empty())
	{
		const detail::candidate_operation c = candidatesNLD.pop();
		Liu_perform_split(c);
	}

	nLiu = obj.num_vertices();
	sequence_d = 2 * (nLiu - ini_num_vertices);

	/// ��ʱdebug�ӵ�
	mesh Mesh1;
	Mesh1.vertices = obj.vertices;
	Mesh1.triangles.clear();
	connectivity.on_triangles([&](const std::array<uint32_t, 3>& t) { Mesh1.triangles.push_back(t); });
	remove_standalone_vertices(Mesh1, connectivity);

	Mesh1.save("initialization_test.obj");

}

/// һ��addһ�� (��Ϊÿ��edge�仯�����Ӧ������halfedge��collapse�������仯)
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
		c.weight = metric.cost_collapse(he.index, to_keep, to_remove);				/// �������edge collapse��cost
		candidatesREM.push(c);
	}
}

void reduction::add_split(const half_edge& he, const bool& he_delaunay_valid)
{
	if (he.is_boundary()) return;										
	if (he_delaunay_valid)													/// ���������delaunay�ģ�����Ҫsplit
	{
		candidatesNLD._delete(he.index);									/// �ñ�֮ǰ������nonDelaunay�ģ����ڱ����delaunay������Ҫ��NLD�е�ɾ��
		candidatesNLD._delete(he.opposite().index);
		return;
	}

	detail::candidate_operation c;
	c.index = he.edge_index();
	c.weight = metric.cost_split(he.index);
	candidatesNLD.push(c);
}

void reduction::perform_flip(const half_edge& he)							/// ����flip��ı߼���neighbour
{
	connectivity.flip_edge(he.index);

	candidatesNLD._delete(he.index);
	candidatesNLD._delete(he.opposite().index);

	stats.on_operation(Flip);								/// TODO: ���Ƕ��ٱ�delaunay�����Ӱ��

	// he���ڵ���������������vertex��one ring��ҪUpdateһ�£�����REM��

	// TODO: �޸�bug
	// �����traverse�Ĺ�����flip�ˣ���traverse����ѭ��

	std::unordered_set<uint32_t> REM2update;
	std::unordered_set<uint32_t> REM_NLD2update;
	const char* option = "edge_ring";
	k_ring(connectivity, option, REM2update, REM_NLD2update, he.index);

	for (uint32_t hx : REM_NLD2update)					// ������traverse_k_ring_edge��ֱ�Ӹ��£���Ϊflip��ı�ṹ�����»ز���start edge
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
		if (!he_delaunay_valid && metric.flip_valid(he.index))										/// �������delaunay_valid����flip��flip�����ܽ���collapse/split
			perform_flip(he);
		else
		{
			add_split(he, he_delaunay_valid);
		}
		add_collapse(he);
	}
}

void reduction::perform_collapse(const detail::candidate_operation& c)			/// ��pirority queue���Ƴ���ɾ���˵�halfedges����������Χ��neighbour
{
	half_edge he = connectivity.handle(c.index);								/// he: to_remove -> to_keep
	uint32_t to_keep = he.vertex();
	
	/// ����Ҫ���µ�halfedge�����������
	/// ��Ҫ���µ����� to_remove �� 2 ring
	/// TO OPTIMIZE: REM��Ҫ����2 ring��NLDֻ��Ҫ����1 ring
	uint32_t to_remove = he.opposite().vertex();
	std::unordered_set<uint32_t> REM2update;
	std::unordered_set<uint32_t> REM_NLD2update;
	const char* option = "vertex_ring";
	k_ring(connectivity, option, REM2update, REM_NLD2update, to_remove);

	stats.on_operation(Collapse);

	/*
	* debug����
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

	for (uint32_t h : connectivity.collapse_edge(c.index))				/// connectivity.collapse_edge()���ص��� index of removed half edges
	{
		candidatesREM._delete(h);										/// ��priority queue���Ƴ���Щɾ���˵ı�
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

/// �ȴ�priority queue��ɾ�����е�Halfedges
/// �ٽ��¼ӵ�halfedge����priority queue
/// �ٸ���neighbourhood
void reduction::perform_split(const detail::candidate_operation& c)
{
	half_edge h = connectivity.handle(c.index);

	auto [nonDelaunay, vertexPosition] = metric.split_position(h);
	obj.vertices.push_back(vertexPosition);											/// ���²���ĵ�����꣬���� obj.vertices ��� vector �����
	
	// ������ͳ��nondelaunay���·�����Ӧ���ò�����
	// stats.num_nonDelaunay += nonDelaunay - 1;									/// -1����Ϊh������non-delaunay����

	uint32_t vertex_index = obj.vertices.size() - 1;
	connectivity.split_edge(h.index, vertex_index);									/// vector �� index ��0��ʼ

	stats.on_operation(Split);

	/// update neighborhoods
	std::unordered_set<uint32_t> REM2update;
	std::unordered_set<uint32_t> REM_NLD2update;
	const char* option = "vertex_ring";
	k_ring(connectivity, option, REM2update, REM_NLD2update, vertex_index);

	for (uint32_t h : REM_NLD2update)					/// ������traverse_k_ring_edge��ֱ�Ӹ��£���Ϊflip��ı�ṹ�����»ز���start edge
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
	/// ���� S���غ����û�õ�����
	/// ----------------------------------------
	std::vector<size_t> S;
	bool _skip = false;
	mesh Mesh;
	X = X.round();

	for (size_t i = 0; i <= X.size(); ++i)
	{
		_skip = false;
		size_t j = 0;
		if (i % 2 == 0)								/// odd, E = Es  (��ΪX�����0��ʼ������iΪevenʱ����Ӧ����odd)
		{
			// ��ʱdebug�ӵ�
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
					// set S(x) = Sdummy and return �ռ�
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
			// ��ʱdebug�ĵ�
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
					// set S(x) = Sdummy and return �ռ�
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
			S.push_back(j);		/// ��Ϊj������index����ʵ�ʴ�����1

		if (candidatesNLD.empty())
		{
			size_t num_collap = 0;
			while (!candidatesREM.empty() && stats.num_vertices > target_num_vertices)			// TODO��������������⣬������������==��
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
			connectivity.on_triangles([&](const std::array<uint32_t, 3>& t) { Mesh.triangles.push_back(t); });			// ��reduce�������mesh��triangle
			remove_standalone_vertices(Mesh, connectivity);

			/// ��ʱdebug�ӵ�
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