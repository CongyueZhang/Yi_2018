#pragma once

#include "mesh.h"
#include <charconv>

struct reduction_metric
{
	using serialized_parameters = std::vector<std::pair<std::string, std::string>>;

	virtual ~reduction_metric()
	{
	}

	virtual serialized_parameters save() const
	{
		return {};
	}

	virtual void load(const serialized_parameters&)
	{
	}

	//去掉了const
	virtual void setup(const mesh&, half_edge_connectivity&)
	{
	}

	virtual bool collapse_still_valid(uint32_t /* half_edge */) const
	{
		return true;
	}

	virtual bool flip_still_valid(uint32_t /* half_edge */) const
	{
		return true;
	}

	//改动了
	virtual double cost_collapse(uint32_t /* half_edge */, uint32_t /* to_keep */, uint32_t /* to_remove */) const
	{
		return 0;
	}

	// 操作结束后影响几个rings结构
	virtual unsigned int post_collapse(uint32_t /* half_edge */, uint32_t /* to_keep */, uint32_t /* to_remove */)
	{
		return 0;
	}

	virtual double cost_flip(uint32_t /* half_edge */) const
	{
		return std::numeric_limits<double>::quiet_NaN();
	}

	virtual unsigned int post_flip(uint32_t /* half_edge */)
	{
		return 0;
	}

	// 自己加的-----------------------------------------------------------------
	virtual double cost_split(uint32_t /* half_edge */) const
	{
		return std::numeric_limits<double>::quiet_NaN();
	}

	virtual bool remove_valid(uint32_t /* half_edge */) 
	{
		return true;
	}

	virtual bool delaunay_valid(uint32_t /* h_index */) const
	{
		return false;
	}

	virtual std::pair<double, Eigen::Vector3d> split_position(const half_edge he) const
	{
		return { 0, Eigen::Vector3d(0,0,0) };
	}

	virtual bool flip_valid(uint32_t /*half_edge*/) const
	{
		return false;
	}

	// ------------------------------------------------------------------------------
protected:
	template<typename T>
	static void serialize(serialized_parameters& p, const char* name, const T& value)
	{
		if constexpr (std::is_same_v<T, bool>)
			p.emplace_back(name, value ? "true" : "false");
		else
			p.emplace_back(name, std::to_string(value));
	}

	template<typename T>
	static void deserialize(const serialized_parameters& p, const char* name, T& value)
	{
		if constexpr (std::is_same_v<T, bool>)
		{
			for (const std::pair<std::string, std::string>& x : p)
				if (x.first == name)
					value = ci_equal(x.second.c_str(), "true")
					|| ci_equal(x.second.c_str(), "yes")
					|| ci_equal(x.second.c_str(), "1")
					|| ci_equal(x.second.c_str(), "on");
		}
		else
		{
			for (const std::pair<std::string, std::string>& x : p)
				if (x.first == name)
					std::from_chars(x.second.data(), x.second.data() + x.second.size(), value);
		}
	}
};
