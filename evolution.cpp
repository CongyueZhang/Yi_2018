#include<time.h>
#include <random>
#include <set>
#include<igl/hausdorff.h>
#include <Eigen/Sparse>
#include <Eigen/Core>

#include "reduce.h"
#include "evolution.h"
#include "mesh.h"
#include "metric.h"


mesh evolution_stream(mesh& obj, const size_t& target_num_vertices, const size_t& Np, const size_t& n_iter, const double& F, const double& Cr, const double& tao, const double& nc, reduction_metric& metric)
{
	printf("Evolution Stream Begin!\n");
	reduction _reduction(obj, metric, target_num_vertices);
	printf("Reduction Initialization Begin!\n");
	_reduction.initialize();
	printf("Reduction Initialization Done!\n");
	const size_t nLiu = _reduction.get_nLiu();
	const size_t Ymax = nLiu - obj.num_vertices();
	const size_t Xmax = nLiu - target_num_vertices;
	const size_t sequence_d = _reduction.get_sequence_d();

	mesh Mesh_min;
	double D_min = std::numeric_limits<double>::infinity();

	double f_k, f_k1;
	// std::vector<std::vector<size_t>> Xk(Np, std::vector<size_t>(sequence_d));
	// std::vector<std::vector<size_t>> Xk1(Np, std::vector<size_t>(sequence_d));
	Eigen::ArrayXXf Xk(Np, sequence_d);		//这些100原本都是Np
	Eigen::ArrayXXf Xk1(Np, sequence_d);
	// std::vector<double> Dk(Np);
	// std::vector<double> Dk1(Np);
	Eigen::ArrayXf Dk(Np);
	Eigen::ArrayXf Dk1(Np);

	/// initilization
	/// 这里的Xk是X0,Dk是D0
	/// ------------------------------------------------------------
	printf("Evolution Initialization Begin!\n");
	for (int j = 0; j<Np; ++j)
	{
		for (int i = 0; i < sequence_d; ++i)
		{
			std::random_device rd;
			std::mt19937 gen(rd());
			std::uniform_real_distribution<> dis(0, 1);
			if (i % 2 == 0)		/// i is odd	(因为数组原因，从0开始的，所以这里判断为even，对应的i是odd)
			{
				Xk1(j,i) = dis(gen) * double(Ymax);		/// 随机产生[0,1]均匀分布的浮点数  TODO: 验证一下
			}
			else                /// otherwise
			{
				Xk1(j,i) = dis(gen) * double(Xmax);

			}
		}

		auto [Mesh_i, S_i] = _reduction.reduce_stream(Xk1.row(j));

		/// Todo: 保证进行完一次reduce_stream后，obj还和原来一样
		if (S_i[1] == -1)
			Dk1[j] = std::numeric_limits<double>::infinity();
		else
			igl::hausdorff(Mesh_i.matrix_vertices(), Mesh_i.matrix_triangles(), obj.matrix_vertices(), obj.matrix_triangles(), Dk1[j]);
		if (D_min > Dk1[j])
		{
			D_min = Dk1[j];
			Mesh_min = Mesh_i;
		}
	}
	printf("Evolution Initialization done!\n");
	/// ---------------------------------------------------------------------
	size_t K = 0;
	
	f_k = D_min;

	double _nc = 0;

	while (true)
	{
		Xk = Xk1;
		Dk = Dk1;
		int j = 0;
		for (int j =0; j<Np; ++j)
		{
			/// Generate a mutative agent X'k,j by mutation
			/// -------------------------------------------
			/// TODO: 优化生成distinct random的部分
			Eigen::ArrayXf Xk_mutation(sequence_d);
			
			/// 生成随机数rand1,2,3和F
			/// rand1,2,3是互不相同且与j不相同的整数
			/// F是range(0,1)中的随机数
			std::set<int> rand3;
			rand3.insert(j);
			srand((unsigned)time(NULL));
			while (rand3.size() != 4)
			{
				rand3.insert(rand()%(Np));
			}
			rand3.erase(j);
			int _rand3[3];
			int i = 0;
			for (auto& n : rand3)
			{
				_rand3[i++] = n;
			}

			/// compute X'k,j
			Xk_mutation = Xk.row(_rand3[0])+ F*(Xk.row(_rand3[1]) - Xk.row(_rand3[2]));

			Eigen::ArrayXf temp_col_i(Np);
			for (int i = 0; i < sequence_d; ++i)
			{
				if (i % 2 != 0)		/// i is odd
				{
					Xk_mutation.row(i) = (Xk_mutation.row(i).max(0)).min(Ymax);
				}
				else                /// otherwise
				{
					Xk_mutation.row(i) = (Xk_mutation.row(i).max(0)).min(Xmax);
				}		
			}

			/// for odd i
			/// 因为index从0开始，所以even row对应着odd i
			Xk_mutation(Eigen::all, Eigen::seq(0, Eigen::last, 2)) = Xk_mutation(Eigen::all, Eigen::seq(0, Eigen::last, 2)).max(0).min(Ymax);

			/// for even i
			Xk_mutation(Eigen::all, Eigen::seq(1, Eigen::last, 2)) = Xk_mutation(Eigen::all, Eigen::seq(1, Eigen::last, 2)).max(0).min(Xmax);

			/// Generate a trial agent X''k,j by crossover
			/// ------------------------------------------
			Eigen::ArrayXf Xk_crossover = Eigen::ArrayXf::Random(sequence_d);
			Xk_crossover = Xk_crossover.abs();
			
			// Todo: To optimize
			for (int i = 0; i < sequence_d; ++i)
			{
				if (Xk_crossover(i) <= Cr)
				{
					Xk_crossover(i) = Xk_mutation(i);
				}
				else
					Xk_crossover(i) = Xk(j,i);
			}
			// i_rand is a random index in [1, d] to ensure that the trial agent 
			// get at least one element from the mutative agent.
			int i_rand = rand() % (sequence_d);
			Xk_crossover(i_rand) = Xk_mutation(i_rand);

			/// Compute the mapped sequence S(X''k,j) and the candidate mesh
			/// ------------------------------------------------------------
			auto[crossMesh, S] = _reduction.reduce_stream(Xk_crossover);
			double D_crossover;
			if (S[1] == -1)
				D_crossover = std::numeric_limits<double>::infinity();
			else
				igl::hausdorff(crossMesh.matrix_vertices(), crossMesh.matrix_triangles(), obj.matrix_vertices(), obj.matrix_triangles(), D_crossover);
			if (D_crossover < Dk(j))
			{
				Xk1(j) = Xk_crossover(j);
				Dk1(j) = D_crossover;
				if (D_crossover < D_min)
				{
					D_min = D_crossover;
					Mesh_min = crossMesh;
				}
			}
			else
			{
				Xk1(j) = Xk(j);
				Dk1(j) = Dk(j);
			}
		}
		// termination condition
		if ((Dk.minCoeff() - Dk1.minCoeff()) / Dk.minCoeff() && _nc++ >= nc || ++K > n_iter)
			break;
		else
			_nc = 0;
	}
	return Mesh_min;
}