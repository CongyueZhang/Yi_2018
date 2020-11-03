#include <charconv>
#include <memory>
#include <chrono>
#include <string>

#include "mesh.h"
#include "reduce.h"
#include "metric_Yi2018.h"
#include "geometry.h"
#include "evolution.h"

inline bool ci_equal(const char* a, const char* b)
{
	if (!a || !b) return false;
	while (*a && *b && std::tolower(*a) == std::tolower(*b))
		a++, b++;
	return !*a && !*b;
}


using metric_create_f = std::unique_ptr<reduction_metric>();

template<typename T>
std::unique_ptr<reduction_metric> make_metric() { return std::make_unique<T>(); }

constexpr const std::pair<const char*, metric_create_f*>  metric_register[] =
{
	{ "Yi_2018",	 make_metric<metric_Yi2018> }
};

std::unique_ptr<reduction_metric> make_metric(const char* name, const reduction_metric::serialized_parameters& params)
{
	for (const std::pair<const char*, metric_create_f*>& m : metric_register)
		if (ci_equal(name, m.first))
		{
			std::unique_ptr<reduction_metric> metric = (*m.second)();
			metric->load(params);
			return metric;
		}
	return nullptr;
}

//计算持续时间
//------------------------------
template<typename C, typename D>
double elapsed(const std::chrono::time_point<C, D>& start, const std::chrono::time_point<C, D>& end)
{
	return std::chrono::duration_cast<std::chrono::microseconds>(end - start).count() / 1000000.0;
}

mesh load_mesh(const char* filename)
{
	using clock = std::chrono::high_resolution_clock;

	const clock::time_point ts = clock::now();
	mesh m(filename);
	const clock::time_point te = clock::now();

	if (m.empty())
	{
		printf("Cannot read triangle mesh '%s'\n", filename);
		exit(1);
	}

	printf("Input mesh:\n"
		"\tvertices  : %zu\n"
		"\tedges     : %zu\n"
		"\ttriangles : %zu\n",
		m.num_vertices(), num_edges(m), m.num_triangles());

	printf("\tload time : %.6f s\n", elapsed(ts, te));
	return m;
}

struct command_line
{
	// Files
	const char* filename_in = nullptr;
	const char* filename_out = nullptr;

	// General options
	unsigned int num_vertices = 0;
	size_t Np = 100, n_iter = 100, Nc = 5;
	double tao = pow(10, -4);
	double F = 0.5, Cr = 0.9;

	command_line(unsigned int argc, char** argv)
	{
		if (argc == 4)
		{
			magic(argc, argv);
			return;
		}

		if (argc < 9)
		{
			printf("ERROR: not enough arguments\n");
			print_help();
			exit(0);
		}

		filename_in = argv[1];
		filename_out = argv[2];
		std::from_chars(argv[3], argv[3] + strlen(argv[3]), num_vertices);
		std::from_chars(argv[4], argv[4] + strlen(argv[4]), Np);
		std::from_chars(argv[5], argv[5] + strlen(argv[5]), n_iter);
		std::from_chars(argv[6], argv[6] + strlen(argv[6]), F);
		std::from_chars(argv[7], argv[7] + strlen(argv[7]), Cr);
		std::from_chars(argv[8], argv[8] + strlen(argv[8]), tao);
		std::from_chars(argv[9], argv[9] + strlen(argv[9]), Nc);
	}

	void magic(unsigned int, char** argv)
	{
		printf("Automagically filling missing parameters\n");
		filename_in = argv[1];
		filename_out = argv[2];
		std::from_chars(argv[3], argv[3] + strlen(argv[3]), num_vertices);
	}

	void print_help() const
	{
		printf("Usage:\n"
			"spectral-collapsing <file_in> <file_out> <num_vertices> <Np> <n_iter=100> <Cr> <tao=e-4> <nc = 5>[<param>=<value>]*\n"
			"\n"
			"\t<file_in>		: mesh file to load and process\n"
			"\t<file_out>		: mesh file to save\n"
			"\t<num_vertices>	: target number of vertices in the final mesh\n"
			"\t<Np>				: evolution population size\n"
			"\t<n_iter>			: evolution iteration number\n"
			"\t<F>				: mutation rate\n"
			"\t<Cr>				: evolution crossover rate\n"
			"\t<tao>			: evolution threshold\n"
			"\t<nc>				: consecutive nc iterations\n");
	}

	void print_parameters() const
	{
		printf("Parameters:\n"
			"\tFilename (in)	: %s\n"
			"\tFilename (out)	: %s\n"
			"\tNum vertices		: %u\n"
			"\tNp				: %u\n"
			"\tN_iter			: %u\n"
			"\tF				: %f\n"
			"\tCr				: %f\n"
			"\tTao				: %f\n"
			"\tNc				: %u\n",
			filename_in, filename_out, num_vertices, Np, n_iter, F, Cr, tao, Nc);
	}
};

double reduce(mesh& m, const command_line& cmd,  reduction_metric& metric)
{
	using clock = std::chrono::high_resolution_clock;
	const clock::time_point ts = clock::now();
	
	m = evolution_stream(m, cmd.num_vertices, cmd.Np, cmd.n_iter, cmd.F, cmd.Cr, cmd.tao, cmd.Nc, metric);

	const clock::time_point te = clock::now();
	const double duration = elapsed(ts, te);
	printf("Reduction time: %.6f s\n", duration);
	return duration;
}

int main(unsigned int argc, char** argv)
{
	const command_line cmd(argc, argv);
	cmd.print_parameters();

	const char* filename_in = cmd.filename_in;
	mesh m = load_mesh(filename_in);

	std::unique_ptr<reduction_metric> metric = make_metric("Yi_2018", {});
	reduce(m, cmd, *metric);

	m.save(cmd.filename_out);
}
