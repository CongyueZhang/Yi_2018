#include "metric_Yi2018.h"
#include "geometry.h"
#include <corecrt_math_defines.h>

#include <CGAL/Circle_3.h>
#include <CGAL/Point_3.h>
#include <CGAL/Cartesian.h>

#include <set>
#include <algorithm>

typedef CGAL::Cartesian<double>  K;
typedef CGAL::Point_3< K >	Point;
typedef CGAL::Circle_3< K >	Circle;

//算cot值
//---------------------
double cot(double angle)
{
	return cos(angle) / sin(angle);
}

/// 求点和线段的交点
/// a, b, p1是圆的三点坐标；p1, p2是线段的两个顶点	
/// 其中p1是圆上的一点
Eigen::Vector3d intersec_l_c(const Eigen::Vector3d a, const Eigen::Vector3d b, const Eigen::Vector3d p1, const Eigen::Vector3d p2)
{
	Circle c = Circle(Point(a[0],a[1],a[2]), Point(b[0], b[1], b[2]), Point(p1[0], p1[1], p1[2]));
	Eigen::Vector3d center((c.center()).x(), (c.center()).y(), (c.center()).z());
	double r = c.squared_radius();

	Eigen::Vector3d l = p2 - p1;				// p1指向p2的向量
	l = l / l.norm();
	Eigen::Vector3d t = center - p1;			// p1指向center的向量

	return 2 * t.dot(l) * l + p1;				// return的是交点的【坐标】

}

/// 求点和线段的交点
/// p1, pe, p3是圆的三点坐标；a, b是线段的两个顶点	数轴以  a->b
/// 其中p3为 a或b
double intersec_l_c(const Eigen::Vector3d p1, const Eigen::Vector3d p2, const Eigen::Vector3d p3, const Eigen::Vector3d a, const Eigen::Vector3d b)
{
	Circle c = Circle(Point(p1[0], p1[1], p1[2]), Point(p2[0], p2[1], p2[2]), Point(p3[0], p3[1], p3[2]));
	Eigen::Vector3d center((c.center()).x(), (c.center()).y(), (c.center()).z());
	double r = c.squared_radius();

	Eigen::Vector3d dir = b - a;				
	dir = dir / dir.norm();						// dir是direction，即a->b的方向向量
	Eigen::Vector3d a2c = center - a;			// a指向center的向量

	return 2 * a2c.dot(dir);					// return的是在 a->b数轴 上的坐标
}

/*
/// sample交点几何
/// input: segment的两个端点a, b
/// output: sample这个segment的点集合
std::unordered_set<Eigen::Vector3d> segmentSet(const Eigen::Vector3d a, const Eigen::Vector3d b)
{
	Eigen::Vector3d segment = b - a;
	double t = segment.norm();
	Eigen::Vector3d d = segment / t;		//d是direction

	std::vector<Eigen::Vector3d> r(20);
	for (int i = 0; i < 20; i++)
	{
		r[i] = a + i*t/20*d;
	}

	std::unordered_set<Eigen::Vector3d> s(r.begin(), r.end());

	return s;
}
*/


void metric_Yi2018::setup(const mesh& object, half_edge_connectivity& connec)
{
	obj = &object;
	connectivity = &connec;
}

///判断halfedge是否delaunay
/// --------------------------------------------------------
bool metric_Yi2018::delaunay_valid(uint32_t h_index) const
{
	auto [angle1, angle2] = opp_angle(*obj, *connectivity, h_index);
	if (angle1 + angle2 > M_PI)
	{
		return false;
	}
	return true;
}

/// 设he是 a→b
/// 返回4个edge中【non】-delaunay的数量
std::pair<double, Eigen::Vector3d> metric_Yi2018::split_position(const half_edge he) const
{
	Eigen::Vector3d a, b, c, d, D1, D2, D3, D4;
	b = obj->vertices[he.vertex()];
	D2 = obj->vertices[he.opposite().next().next().opposite().next().vertex()];
	d = obj->vertices[he.next().vertex()];
	D3 = obj->vertices[he.next().opposite().next().vertex()];
	a = obj->vertices[he.next().next().vertex()];
	D4 = obj->vertices[he.next().next().opposite().next().vertex()];
	c = obj->vertices[he.opposite().next().vertex()];
	D1 = obj->vertices[he.opposite().next().opposite().next().vertex()];

	double pI0_1, pI0_2, pI1, pI2, pI3, pI4, pa, pb;

	// 以a->b建立数轴，a为原点
	pa = 0;
	pb = (b - a).norm();
	pI0_1 = intersec_l_c(c, d, a, a, b);				//靠近a点的
	pI0_2 = intersec_l_c(c, d, b, a, b);				//靠近b点的
	pI1 = intersec_l_c(c, D1, a, a, b);
	pI2 = intersec_l_c(c, D2, b, a, b);
	pI3 = intersec_l_c(d, D3, b, a, b);
	pI4 = intersec_l_c(d, D4, a, a, b);

	std::pair<double, double> I0(pI0_1, pI0_2);
	std::vector<std::pair<uint32_t, uint32_t>> I(4);
	I.push_back(std::pair(pI1, pb));
	I.push_back(std::pair(pa, pI2));
	I.push_back(std::pair(pa, pI3));
	I.push_back(std::pair(pI4, pb));

	for (std::vector<std::pair<uint32_t, uint32_t>>::iterator iter = I.begin(); iter != I.end(); )
	{
		if (iter->second < I0.first || iter->first > I0.second)	iter = I.erase(iter);
		else
		{
			if (iter->second > I0.second)	iter->second = I0.second;
			if (iter->first < I0.first)	iter->first = I0.first;
			iter++;
		}
	}

	std::pair<double, int> pmax((I0.first+I0.second)/2.0, 0), p;
	double delta = pb / 20;
	for (int i = 0; i < 20; ++i)
	{
		p.first = delta * i;
		p.second = 0;
		for (std::pair<uint32_t, uint32_t> s : I)
		{
			if (p.first > s.first && p.first < s.second)	++p.second;
		}
		if (p.second > pmax.second)	pmax = p;
		if (pmax.second == 4)	break;
	}

	return { 4 - pmax.second, a + (b - a) / (b - a).norm() * pmax.first };
}


// 计算对角的cot和
// NLD queue 的 key
//------------------------------------------------------
double metric_Yi2018::cost_split(uint32_t h_index) const
{
	auto [angle1, angle2] = opp_angle(*obj, *connectivity, h_index);
	return cot(angle1) + cot(angle2);
}

// 检测是否可以flip
//----------------------------------------------------
bool metric_Yi2018::flip_valid(uint32_t h_index) const
{
	Eigen::Vector3d n1, n2;
	const half_edge h = connectivity->handle(h_index);

	const auto [v1, v2] = connectivity->edge_vertices(h_index);		//h是 v1 -> v2
	Eigen::Vector3d h_vector = obj->vertices[v2] - obj->vertices[v1];

	n1 = h_vector.cross(obj->vertices[h.next().vertex()] - obj->vertices[h.vertex()]);
	n2 = h_vector.cross(obj->vertices[h.opposite().next().vertex()] - obj->vertices[h.opposite().vertex()]);

	double dihedralAngle = abs(n1.dot(n2) / (n1.norm() * n2.norm()));

	return dihedralAngle == 1.0;		//问题：精读会不会产生影响？？
}

// 检测是否removable
//------------------------------------------------
bool metric_Yi2018::remove_valid(uint32_t h_index)
{
	// 两种方法：
	// ①将这一部分mesh取出单独进行操作，判断后这一部分mesh丢掉
	// ②先collapse，然后再split复原

	// ②
	// 将即将被collapse的Halfedges的data保存下来
	// 新想法：除了delete掉的halfedge，射入to_remove的halfedge也变了

	// 在collapse前把原data存下来
	//------------------------------------------
	std::vector<half_edge_data> half_edges_data;						
	std::vector<uint32_t> half_edges_index;

	half_edge he = connectivity->handle(h_index);							// he: v -> v'   v是to_remove
	uint32_t start = uint32_t(-1);
	while (he.is_valid() && he.index != start)								// 遍历v的one ring
	{
		if (start == uint32_t(-1)) 
			start = he.index;

		half_edges_data.push_back(connectivity->half_edges[he.index]);		// 将所有入射v的halfedge的data存起来
		half_edges_index.push_back(he.index);

		he = he.next();

		half_edges_data.push_back(connectivity->half_edges[he.index]);		// 将所有入射v的halfedge的data存起来
		half_edges_index.push_back(he.index);

		he = he.next();

		half_edges_data.push_back(connectivity->half_edges[he.index]);		// 将所有入射v的halfedge的data存起来
		half_edges_index.push_back(he.index);

		he = he.opposite();
	}

	connectivity->collapse_edge_test(h_index);
	
	// 检测collapse后是否满足条件
	// ------------------------
	start = uint32_t(-1);
	bool validity = true;
	for (uint32_t index: half_edges_index)
	{
		if (!connectivity->half_edges[index].is_valid())	continue;		//说明该边已经被删除
		
		/// 检测这个he是否是delaunay的
		he = connectivity->handle(index);
		/// todo: 会有edge的重复遍历(因为每个edge有两个halfedge)   看下怎么优化减少重复   (可以通过visited_edges？)
		if (!(delaunay_valid(he.index)))	
		{
			validity = false;
			break;
		}
	}
	
	/// 用保存的data复原
	/// -----------------------------------------
	for (int i = 0; i<half_edges_data.size(); ++i)
	{
		connectivity->half_edges[half_edges_index[i]] = half_edges_data[i];
	}

	// TODO:  v' 的 vertex_to_half_edge 也要改回去
	/*
	* 将collapse_edge_test中所有可能会改变vertex_to_half_edge的代码都删除了
	he = connectivity->handle(h_index);
	connectivity->vertex_to_half_edge[he.opposite().vertex()] = he.index;

	connectivity->vertex_to_half_edge[he.vertex()] = he.opposite().index;
	*/


	return validity;
}

/// todo: 不知道这个qem是否能达到要求
double metric_Yi2018::cost_collapse(uint32_t half_edge, uint32_t to_keep, uint32_t to_remove) const
{
	/*
	if (lock_boundaries)
	{
		const std::pair<uint32_t, uint32_t> v = connectivity->edge_vertices(half_edge);
		if (connectivity->is_boundary_vertex(v.first) || connectivity->is_boundary_vertex(v.second))
			return std::numeric_limits<double>::infinity();
	}	
	*/

	const Eigen::Vector3d fallback = obj->vertices[to_keep];
	Eigen::Matrix4d q_toKeep, q_toRemove;
	const std::pair<Eigen::Vector3d, double> result = qem_merge(vertex_quadric((const mesh&)obj, (const half_edge_connectivity&)connectivity, to_keep), vertex_quadric((const mesh&)obj,(const half_edge_connectivity&)connectivity, to_remove), fallback);
	return result.second;
}

