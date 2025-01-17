#pragma once

#include <vector>

template<typename T, typename F1 = std::less<T>, typename F2 = std::equal_to<T>>
void sort_unique_inplace(std::vector<T>& v, F1 is_less = F1(), F2 is_same = F2())
{
	std::sort(v.begin(), v.end(), is_less);
	v.erase(std::unique(v.begin(), v.end(), is_same), v.end());
}

template<typename Ex, typename T, typename F1 = std::less<T>, typename F2 = std::equal_to<T>>
void sort_unique_inplace(Ex exec, std::vector<T>& v, F1 is_less = F1(), F2 is_same = F2())
{
	std::sort(exec, v.begin(), v.end(), is_less);
	v.erase(std::unique(v.begin(), v.end(), is_same), v.end());		
	// unique返回的是the new location of the element that followed the last element erased by the function call
	// 因此上行代码的功能是给v中的elements排序，并移除v中的重复elements
}

/// Arrays must be sorted
template<typename T>
std::vector<T> merge_unique(const std::vector<T>& a, const std::vector<T>& b)
{
	std::vector<T> c = a;
	c.insert(c.end(), b.begin(), b.end());
	sort_unique_inplace(c);
	return c;
}
