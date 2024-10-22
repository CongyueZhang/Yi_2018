#include "metric_Yi2018.h"
#include "geometry.h"
#include "traversal.h"
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
double intersec_l_c(const Eigen::Vector3d& p1, const Eigen::Vector3d& p2, const Eigen::Vector3d& p3, const Eigen::Vector3d& a, const Eigen::Vector3d& b)
{
	Circle c = Circle(Point(p1[0], p1[1], p1[2]), Point(p2[0], p2[1], p2[2]), Point(p3[0], p3[1], p3[2]));
	Eigen::Vector3d center((c.center()).x(), (c.center()).y(), (c.center()).z());

	Eigen::Vector3d dir = (b - a).normalized();			// dir是direction，即a->b的方向向量
	Eigen::Vector3d a2c = center - a;					// a指向center的向量

	double I = 2 * a2c.dot(dir);

	return I;					// return的是I的长度
}

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
	// 检测sliver
	if (angle1 == 0 || angle2 == 0 || std::isnan(angle1) || std::isnan(angle2))
	{
		return false;
	}
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

	// 以a->b建立数轴，a为原点，坐标为0；b坐标为ab向量的长度
	pa = 0;
	pb = (b - a).norm();
	pI0_1 = intersec_l_c(c, d, a, a, b);					//靠近a点的
	pI0_2 = pb - intersec_l_c(c, d, b, b, a);				//靠近b点的
	pI1 = intersec_l_c(c, D1, a, a, b);
	pI2 = pb - intersec_l_c(c, D2, b, b, a);
	pI3 = pb - intersec_l_c(d, D3, b, b, a);
	pI4 = intersec_l_c(d, D4, a, a, b);

	std::pair<double, double> I0(pI0_2, pI0_1);
	std::vector<std::pair<double, double>> I(4);
	I[0] = std::make_pair(pI1, pb);
	I[1] = std::make_pair(pa, pI2);
	I[2] = std::make_pair(pa, pI3);
	I[3] = std::make_pair(pI4, pb);

	if (I0.second > pb)	I0.second = pb;
	if (I0.first < pa)	I0.first = pa;

	for (std::vector<std::pair<double, double>>::iterator iter = I.begin(); iter != I.end(); )
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

	// 在I0中采样
	size_t N = 51;				// 采样个数
	double delta = (I0.second - I0.first) / N;
	int j;
	for (int i = N/2+1; i < N; ++i)						// 从中间开始，更容易找到，找到的也更合适
	{
		// 向大的
		p.first = delta * i + I0.first;
		p.second = 0;
		for (std::pair<uint32_t, uint32_t> s : I)
		{
			if (p.first > s.first && p.first < s.second)	++p.second;
		}
		if (p.second > pmax.second)	pmax = p;
		if (pmax.second == 4)	break;

		// 向小的
		j = N - i;
		p.first = delta * j + I0.first;
		p.second = 0;
		for (std::pair<uint32_t, uint32_t> s : I)
		{
			if (p.first > s.first && p.first < s.second)	++p.second;
		}
		if (p.second > pmax.second)	pmax = p;
		if (pmax.second == 4)	break;
	}

	return { 4 - pmax.second, a + (b - a).normalized() * pmax.first };
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

	const auto [v1, v2] = connectivity->edge_vertices(h_index);		//h: v1 -> v2
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
	// 在collapse前把原data存下来
	//------------------------------------------
	std::unordered_map<uint32_t, half_edge_data> to_restore;

	half_edge he = connectivity->handle(h_index);							// he: v -> v'   v是to_remove
	uint32_t to_remove = he.opposite().vertex();
	std::unordered_set<uint32_t> visited_edges;
	uint32_t start = uint32_t(-1);

	std::vector<uint32_t> polygon_vertices;

	const char* option = "halfedge";
	k_ring(*connectivity, 1, option, to_remove, [&](const half_edge& h) {
			to_restore[h.index] = (*connectivity).half_edges[h.index];
			if(h.opposite().vertex() == to_remove)
				polygon_vertices.push_back(h.vertex());
		});

	// 检测内角和
	int n = polygon_vertices.size();
	for (int i = 0; i < n; ++i)
	{
		Eigen::Vector3d e0, e1, e2;
		double _angle1, _angle2, _angle;
		uint32_t index0 = i-1, index1 = i, index2 = i+1;
		if (i == 0)
			index0 = n - 1;
		else if (i == n - 1)
			index2 = 0;
		e0 = obj->vertices[polygon_vertices[index0]] - obj->vertices[polygon_vertices[index1]];
		e2 = obj->vertices[polygon_vertices[index2]] - obj->vertices[polygon_vertices[index1]];
		e1 = obj->vertices[to_remove] - obj->vertices[polygon_vertices[index1]];
		_angle = getAngle(e0, e2);
		_angle1 = getAngle(e0, e1);
		_angle2 = getAngle(e1, e2);
		if (_angle1 + _angle2 > M_PI)
			return false;
	}

	connectivity->collapse_edge_test(h_index);
	
	// 检测collapse后是否满足条件
	// ------------------------
	bool validity = true;

	// 检测是否是delaunay的
	if (validity)
	for (auto _data: to_restore)
	{
		if (!connectivity->handle(_data.first).is_valid())	continue;		//说明该边已经被删除

		if (visited_edges.count(_data.first) == 0)
		{
			/// 检测这个he是否是delaunay的
			if (!(delaunay_valid(_data.first)))
			{
				validity = false;
				break;
			}
		}
		visited_edges.insert(_data.first);
		visited_edges.insert(connectivity->handle(_data.first).opposite().index);
	}
	
	/// 用保存的data复原
	/// -----------------------------------------
	for (auto _data: to_restore)
	{
		connectivity->half_edges[_data.first] = _data.second;
	}

	return validity;
}

/// todo: 不知道这个qem是否能达到要求
double metric_Yi2018::cost_collapse(uint32_t half_edge, uint32_t to_keep, uint32_t to_remove) const
{
	if (lock_boundaries)
	{
		const std::pair<uint32_t, uint32_t> v = connectivity->edge_vertices(half_edge);
		if (connectivity->is_boundary_vertex(v.first) || connectivity->is_boundary_vertex(v.second))
			return std::numeric_limits<double>::infinity();
	}	


	const Eigen::Vector3d fallback = obj->vertices[to_keep];
	Eigen::Matrix4d q_toKeep, q_toRemove;
	const std::pair<Eigen::Vector3d, double> result = qem_merge(vertex_quadric((const mesh&)*obj, (const half_edge_connectivity&)*connectivity, to_keep), vertex_quadric((const mesh&)*obj,(const half_edge_connectivity&)*connectivity, to_remove), fallback);
	return result.second;
}

