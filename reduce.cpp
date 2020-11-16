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

/// ��traversal.h��copy����
template<typename F>
void reduction::traverse_k_ring_edge(unsigned int k, uint32_t v, F f)
{
	visited_vertices.insert(v);
	traverse_1_ring(v, [&](const half_edge& h)				// ������Ǵ���tranverse_1_ring�ĺ���
		{
			if (visited_edges.count(h.index) == 0) f(h);			// ��ֹ�ظ�����ĳ��Halfedge
			visited_edges.insert(h.index);
			visited_edges.insert(h.opposite().index);

			const uint32_t vert = h.vertex();
			if (k > 1 && visited_vertices.count(vert) == 0) traverse_k_ring_edge(k - 1, vert, f);		// ͬʱ����v��ÿһ���Բ��
		});
}

template<typename F>
void reduction::traverse_1_ring(uint32_t v, F f)										// ����vertex�а���v�������ε�����halfedge������ÿ��halfedgeִ�� f ����
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

template<typename F>
void reduction::traverse_2_triangles(half_edge h, F f)						
{
	half_edge g;
	g = h.opposite();
	f(h);
	f(h = h.next());
	f(h = h.next());
	f(g = g.next());
	f(g = g.next());
}

void reduction::Liu_process_he(const uint32_t& hx)
{
	half_edge he = connectivity.handle(hx);
	if (!metric.delaunay_valid(hx) && metric.flip_valid(hx))
	{
		connectivity.flip_edge(hx);
		candidatesNLD._delete(std::min(hx, he.opposite().index));
		// ��������������Ӱ��  ��hx���ڵ������������е������ߣ�
		std::vector<uint32_t> he2update;
		traverse_2_triangles(he, [&](half_edge& h) {
			he2update.push_back(h.index);
			});
		for (uint32_t h : he2update)		
		{
			if (connectivity.handle(h).is_valid())
				Liu_process_he(h);
		}
	}
	add_split(he);
}

void reduction::Liu_perform_split(const detail::candidate_operation& c)
{
	half_edge h = connectivity.handle(c.index);

	auto [nonDelaunay, vertexPosition] = metric.split_position(h);
	obj.vertices.push_back(vertexPosition);											/// ���²���ĵ�����꣬���� obj.vertices ��� vector �����

	uint32_t vertex_index = obj.vertices.size() - 1;
	connectivity.split_edge(h.index, vertex_index);				
	stats.on_operation(Split);

	// ��ʱdebug�ӵ�
	// ------------------------
	//system("pause");
	/*
	mesh Mesh1;
	Mesh1.vertices = obj.vertices;
	connectivity.on_triangles([&](const std::array<uint32_t, 3>& t) { Mesh1.triangles.push_back(t); });			// ��reduce�������mesh��triangle��
	remove_standalone_vertices(Mesh1, connectivity);

	Mesh1.save("delaunay_test_split" + std::to_string(stats.num_split) + "_REM" + std::to_string(candidatesNLD.size()) + ".obj");	
	*/
	// ------------------------

	/// update neightbourhood
	visited_edges.clear();
	visited_vertices.clear();
	std::vector<uint32_t> he2update;
	traverse_k_ring_edge(1, vertex_index, [&](half_edge he) {
		he2update.push_back(h.index);
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
	candidatesNLD.reserve(std::ceil(connectivity.num_half_edges() / 2.0), std::ceil(connectivity.num_half_edges() / 2.0));

	/// �ڽ��в���ǰ�Ȱ�������Ҫflip�ı�flip��
	/// ͳ��non_delaunay_valid�ġ�edge��������
	/// Ѱ��l_max, l_min, angle_min
	size_t _number = connectivity.num_half_edges();

	// �Ȱ���filp��ȫflip��
	// ������NLD queue
	for (uint32_t hx = 0; hx < _number; hx++)
	{
		const half_edge he = connectivity.handle(hx);
		const uint32_t ho = he.opposite().index;
		if (!he.is_boundary() && hx > ho) continue;						/// he.index > ho ��Ϊ�˷�ֹ�ظ�����һ��edge

		Liu_process_he(hx);
	}

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

void reduction::add_split(const half_edge& he)
{
	if (he.is_boundary()) return;
	if (metric.delaunay_valid(he.index))										/// ���������delaunay�ģ�����Ҫsplit
	{
		candidatesNLD._delete(std::min(he.index, he.opposite().index));		/// �ñ�֮ǰ������nonDelaunay�ģ����ڱ����delaunay������Ҫ��NLD�е�ɾ��
		return;
	}

	detail::candidate_operation c;
	c.index = std::min(he.index,he.opposite().index);
	c.weight = metric.cost_split(he.index);
	candidatesNLD.push(c);
}

/// ����halfedge����priority queue
void reduction::add_operation(const half_edge& he)
{
	add_collapse(he);
	add_split(he);
}

void reduction::perform_flip(const half_edge& he)								/// ����flip��ı߼���neighbour
{
	connectivity.flip_edge(he.index);

	candidatesNLD._delete(std::min(he.index, he.opposite().index));

	stats.on_operation(Flip);								/// TODO: ���Ƕ��ٱ�delaunay�����Ӱ��

	// he���ڵ���������������vertex��one ring��ҪUpdateһ�£�����REM��
	visited_edges.clear();
	visited_vertices.clear();

	// TODO: �޸�bug
	// �����traverse�Ĺ�����flip�ˣ���traverse����ѭ��
	std::vector<uint32_t> he2update;
	traverse_k_ring_edge(1, he.vertex(),										
		[&](half_edge h) {
			he2update.push_back(h.index);
		});
	traverse_k_ring_edge(1, he.opposite().vertex(),	
		[&](half_edge h) {
			he2update.push_back(h.index);
		});
	traverse_k_ring_edge(1, he.next().vertex(),
		[&](half_edge h) {
			he2update.push_back(h.index);
		});
	traverse_k_ring_edge(1, he.opposite().next().vertex(),
		[&](half_edge h) {
			he2update.push_back(h.index);
		});

	for (uint32_t h : he2update)					// ������traverse_k_ring_edge��ֱ�Ӹ��£���Ϊflip��ı�ṹ�����»ز���start edge
	{
		if (connectivity.handle(h).is_valid())
			process_he(connectivity.handle(h));
	}
	// 
	/*
	* ������ЩӦ���Ǽ��NLD��
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
		if (!metric.delaunay_valid(he.index) && metric.flip_valid(he.index))										/// �������delaunay_valid����flip��flip�����ܽ���collapse/split
			perform_flip(he);
		else 
			add_operation(he);
	}
}

void reduction::perform_collapse(const detail::candidate_operation& c)			/// ��pirority queue���Ƴ���ɾ���˵�halfedges����������Χ��neighbour
{
	half_edge he = connectivity.handle(c.index);						/// he�Ǵ�Ҫ�Ƴ��ĵ�����
	uint32_t to_keep = he.vertex();

	// TODO: To optimize
	// ��he2updateֱ�Ӵ��unordered_set
	std::vector<uint32_t> he2update;
	
	/// ����Ҫ���µ�halfedge�����������
	/// ��Ҫ���µ����� to_remove �� 2 ring
	/// TO OPTIMIZE: REM��Ҫ����2 ring��NLDֻ��Ҫ����1 ring
	uint32_t to_remove = he.opposite().vertex();
	visited_edges.clear();
	visited_vertices.clear();
	traverse_k_ring_edge(2, to_remove,										/// collapse��split���ڶ������edge������Halfedge������Ҫ��ֹ�ظ�����ͬһ��edge
		[&](half_edge h) {
				he2update.push_back(h.index);										
		});	

	stats.on_operation(Collapse);

	for (uint32_t h : connectivity.collapse_edge(c.index))				/// connectivity.collapse_edge()���ص��� index of removed half edges
	{
		candidatesREM._delete(h);										/// ��priority queue���Ƴ���Щɾ���˵ı�
		candidatesNLD._delete(std::min(h, connectivity.handle(h).opposite().index));
	}

	for (uint32_t h : he2update)
	{
		if (connectivity.handle(h).is_valid())							/// ����δɾ����Halfedge
			process_he(connectivity.handle(h));
	}	

/// ��ʱdebug�ӵ�
/// ----------------------------------------
// Plot the mesh
	/*
	mesh Mesh;
	Mesh.vertices = obj.vertices;
	connectivity.on_triangles([&](const std::array<uint32_t, 3>& t) { Mesh.triangles.push_back(t); });			// ��reduce�������mesh��triangle��
	remove_standalone_vertices(Mesh, connectivity);

	Mesh.save("oneNLD.obj");		
	*/
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

	/// update neightbourhood
	/// ����: process_he������ĺ���������ֱ�Ӵ�������connectivity,stats�ȵ�ֻ����������
	/// ���õĽ������������غ������ó����ŵ���ߵ�����
	/// TO OPTIMIZE: REM��Ҫ����2 ring��NLDֻ��Ҫ����1 ring
	visited_edges.clear();
	visited_vertices.clear();
	std::vector<uint32_t> he2update;
	traverse_k_ring_edge(2, vertex_index, [&](half_edge he) {
		he2update.push_back(h.index);
		});

	for (uint32_t h : he2update)					// ������traverse_k_ring_edge��ֱ�Ӹ��£���Ϊflip��ı�ṹ�����»ز���start edge
	{
		if (connectivity.handle(h).is_valid())						
			process_he(connectivity.handle(h));
	}
}

std::pair<mesh, std::vector<size_t>> reduction::reduce_stream(Eigen::ArrayXf X)
{
	stats.clear();
	stats.num_vertices = ini_num_vertices;

	/// delete connectivity;			///TODO: ȷ��һ���Ƿ���Ҫ�����ͷſռ�
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
	/// ���� S����ʲô�õģ����غ����û�õ�����
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
			for (; j < X[i]; ++j)
			{
				// ���⣺ ����һ��operationǰ�Ƿ�Ҫ�жϣ���
				// ����������ĳ������ж��ٲ���

				if (candidatesNLD.empty())
				{
					_skip = true;
					break;
				}
				if (stats.num_split >= nLiu)
				{
					// set S(x) = Sdummy and return �ռ�
					S.resize(1, -1);
					return { Mesh, S };			// dummy
				}

				const detail::candidate_operation c = candidatesNLD.pop();
				perform_split(c);
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
					// set S(x) = Sdummy and return �ռ�
					S.resize(1, -1);
					return { Mesh, S };			// dummy
				}

				const detail::candidate_operation c = candidatesREM.pop();
				perform_collapse(c);
			}
		}

		if (_skip)
			S.push_back(j);
		else
			S.push_back(j+1);		/// ��Ϊj������index����ʵ�ʴ�����1

		if (candidatesNLD.empty())
		{
			size_t num_collap = 0;
			while (!candidatesREM.empty() && stats.num_vertices != target_num_vertices)			// TODO��������������⣬������������==��
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
			Mesh.triangles.clear();
			connectivity.on_triangles([&](const std::array<uint32_t, 3>& t) { Mesh.triangles.push_back(t); });			// ��reduce�������mesh��triangle��
			remove_standalone_vertices(Mesh, connectivity);

			/// ��ʱdebug�ӵ�
			mesh Mesh1;
			Mesh1.vertices = obj.vertices;
			Mesh1.triangles.clear();
			connectivity.on_triangles([&](const std::array<uint32_t, 3>& t) { Mesh1.triangles.push_back(t); });		
			remove_standalone_vertices(Mesh1, connectivity);

			Mesh1.save("test.obj");

			return { Mesh, S };
		}
	}
	S.resize(1, -1);
	return { Mesh, S };
}