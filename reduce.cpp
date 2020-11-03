#include "reduce.h"
#include "mesh.h"
#include "ranges.h"
#include "geometry.h"
#include "half_edges.h"

#include <execution>
#include<math.h>

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
void reduction::traverse_1_ring_edge(uint32_t v, F f)
{
	std::unordered_set<uint32_t> visited_edges;
	traverse_1_ring(v, [&](const half_edge& h)				// ������Ǵ���tranverse_1_ring�ĺ���
		{
			if (visited_edges.count(h.index) == 0) f(h);			// ��ֹ�ظ�����ĳ��Halfedge
			visited_edges.insert(h.index);
			visited_edges.insert(h.opposite().index);
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

/// TODO: ����boundary�����
void reduction::initialize()
{
	std::function<uint32_t(const detail::candidate_operation& c)> candidate_index = [](const detail::candidate_operation& c)->uint32_t {
		return c.index;
	};
	candidatesNLD.setup_indexof(candidate_index);
	candidatesREM.setup_indexof(candidate_index);

	double l_max = 0, l_min = std::numeric_limits<double>::infinity(), angle_max = 0;

	/// �ڽ��в���ǰ�Ȱ�������Ҫflip�ı�flip��
	/// ͳ��non_delaunay_valid�ġ�edge��������
	/// Ѱ��l_max, l_min, angle_max
	size_t _number = connectivity.num_half_edges();
	for (uint32_t hx = 0; hx < _number; hx++)	/// TODO: parallel
	{
		const half_edge he = connectivity.handle(hx);
		const uint32_t ho = he.opposite().index;
		if (!he.is_boundary() && he.index < ho) continue;						/// he.index < ho ��Ϊ�˷�ֹ�ظ�����һ��edge

		//if (!metric.delaunay_valid(hx))
		//{
			//if (metric.flip_valid(hx))
			//	connectivity->flip_edge(hx);
			//else
			//connectivity->half_edges[hx].delaunay_valid = false;
			//++stats.num_nonDelaunay;
		//}
		Eigen::Vector3d v1, v2;													/// he����������
		v1 = obj.vertices[he.vertex()];
		v2 = obj.vertices[he.opposite().vertex()];
		double l_v1v2 = (v1 - v2).norm();
		if (l_v1v2 < l_min)	l_min = l_v1v2;
		if (l_v1v2 > l_max) l_max = l_v1v2;

		auto [angle1, angle2] = opp_angle(obj, connectivity, hx);
		if (angle1 > angle_max)	angle_max = angle1;
		if (angle2 > angle_max) angle_max = angle2;
	}

	nLiu = ceil(ini_num_vertices * l_max / (l_min * sin(angle_max) * sin(angle_max)));
	sequence_d = 2 * (nLiu - ini_num_vertices);									/// operation sequence�ĳ���
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
	if (!he.is_boundary()) return;
	if (metric.delaunay_valid(he.index))										/// ���������delaunay�ģ�����Ҫsplit
	{
		candidatesNLD._delete(he.index);		/// �ñ�֮ǰ������nonDelaunay�ģ����ڱ����delaunay������Ҫ��NLD�е�ɾ��
		return;
	}

	detail::candidate_operation c;
	c.index = he.index;
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

	stats.on_operation(Flip);								/// TODO: ���Ƕ��ٱ�delaunay�����Ӱ��

	// ����he���ڵ������������е�����halfedges
	process_he(he);
	process_he(he.next());
	process_he(he.next().next());
	process_he(he.opposite().next());
	process_he(he.opposite().next().next());
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
	half_edge he = connectivity.handle(c.index);								/// he�Ǵ�Ҫ�Ƴ��ĵ�����

	std::vector<uint32_t> he2update;

	/// ����Ҫ���µ�halfedge�����������
	/// ��Ҫ���µ����� to_remove Ϊ��������������ε�����halfedge
	/// ��������non_delaunay�ıߣ�����Ҫ��NLD queue���Ƴ���������num_nonDelaunay ����Ϊcollapse����ЩhalfedgeҪô��ɾ����Ҫô���delaunay�ˣ�
	uint32_t to_remove = he.opposite().vertex();
	traverse_1_ring_edge(to_remove,										/// collapse��split���ڶ������edge������Halfedge������Ҫ��ֹ�ظ�����ͬһ��edge
		[&](half_edge h) {
				he2update.push_back(h.index);										
		});

	stats.on_operation(Collapse);

	for (uint32_t h : connectivity.collapse_edge(c.index))							/// connectivity.collapse_edge()���ص��� index of removed half edges
	{
		candidatesREM._delete(h);				/// ��priority queue���Ƴ���Щɾ���˵ı�
		candidatesNLD._delete(h);
	}
	for (uint32_t h : he2update)
	{
		if (connectivity.handle(h).is_valid())										/// ����δɾ����Halfedge
			process_he(connectivity.handle(h));
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
	connectivity.split_edge(h.index, vertex_index);								/// vector �� index ��0��ʼ

	stats.on_operation(Split);

	/// update neightbourhood
	/// ����: process_he������ĺ���������ֱ�Ӵ�������connectivity,stats�ȵ�ֻ����������
	/// ���õĽ������������غ������ó����ŵ���ߵ�����
	traverse_1_ring_edge(vertex_index, [&](half_edge he) {
		process_he(he);
		});
}

std::pair<mesh, std::vector<size_t>> reduction::reduce_stream(Eigen::ArrayXf X)
{
	stats.clear();
	stats.num_vertices = ini_num_vertices;

	/// delete connectivity;			///TODO: ȷ��һ���Ƿ���Ҫ�����ͷſռ�
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
		if (!he.is_boundary() && he.index > he.opposite().index) continue;								/// he.index > ho ��ֹ�ظ�����һ��edge

		process_he(he);																
	}

	/// Algorithm2
	/// ���� S����ʲô�õģ����غ����û�õ�����
	/// ----------------------------------------
	std::vector<size_t> S;
	bool _skip = false;
	mesh Mesh;
	X = X.round();
	for (size_t i = 1; i <= X.size(); ++i)
	{
		_skip = false;
		size_t j = 0;
		if (i % 2 == 1)								/// odd, E = Es
		{
			for (; j < X(i); ++j)
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
					return { NULL, S };
				}

				const detail::candidate_operation c = candidatesNLD.pop();
				perform_split(c);
				++stats.num_split;
				++stats.num_vertices;
			}
		}

		if (i % 2 == 0)		// even, E = Ec
		{
			for (; j < X(i); ++j)
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
					return { NULL, S };
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
			S.push_back(j+1);		/// ��Ϊj������index����ʵ�ʴ�����1

		if (candidatesNLD.empty())
		{
			size_t num_collap = 0;
			if (stats.num_vertices != target_num_vertices)
			{
				++num_collap;
				const detail::candidate_operation c = candidatesREM.pop();
				perform_collapse(c);
			}
			if (S.size() % 2 == 0)
			{
				S.push_back(0);
				S.push_back(num_collap);
			}
			else
				S.push_back(num_collap);

			Mesh.vertices = obj.vertices;
			connectivity.on_triangles([&](const std::array<uint32_t, 3>& t) { Mesh.triangles.push_back(t); });			// ��reduce�������mesh��triangle��
			remove_standalone_vertices(Mesh, connectivity);
			obj.vertices.resize(ini_num_vertices);
			return { Mesh, S };
		}
	}
	S.resize(1, -1);
	return { NULL, S };
}