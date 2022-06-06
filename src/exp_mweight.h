#pragma once

#include <limits>
#include "./util.h"
#include "./error.h"
#include "./string.h"
#include "./mpi.h"
#include "./time.h"

template<typename T>
class NestLoopStack {
private:
	static const int max_nest_num = 10;
	int nest_num;
public:
	std::vector< std::vector< T > > data;
	int loop_size[max_nest_num];
	NestLoopStack() :nest_num(0) {}
	~NestLoopStack() { fini(); }
	void init(std::initializer_list<int> _loop_size) {
		fini();
		size_t total_loop_size = 1;
		for (auto const & __loop_size : _loop_size) {
			total_loop_size *= __loop_size;
			loop_size[nest_num++] = __loop_size;
		}
		assert(nest_num < max_nest_num);
		assert(total_loop_size < (std::numeric_limits<double>::max)());
		ERROR_CHECK(total_loop_size>(std::numeric_limits<double>::max)(), "large...loop");
		data.resize(total_loop_size);
	}
	void fini() {
		data.clear();
		nest_num = 0;
	}
	std::vector<T> &wptr(std::initializer_list<int> _index) {
		assert(_index.size() == nest_num);
		int index = 0;
		int i = 0;
		for (auto const &__index : _index) {
			assert(__index<loop_size[i]);
			index = index*loop_size[i] + __index;
			i++;
		}
		return data[index];
	}
	void push_back(std::initializer_list<int> _index, T t) {
		auto& _data = wptr(_index);
		_data.push_back(t);
	}
	size_t size() {
		size_t total_size = data.size();
		for (auto const& _data : data)
			total_size += _data.size();
		return total_size;
	}
};

// this matrix for exp_Mweight(nr_images*nr_classes*nr_directions*nr_psi*nr_over_rot*nr_translations*nr_over_trans)
// this is the matrix size for different oversampling and healpix-order:
// healpix_order,  oversampling,nr_images, nr_class,  nr_direction,           nr_rotation,    nr_over_rot      nr_trans, nr_over_trans,angle_precision,               size,          size(GB)
//             1,             1,        8,        8,            48(12*4)  ,            12,             8,            81,             4,            15,            95551488,      0.711914(GB)
//             2,             1,        8,        8,           192(12*4^2),            24,             8,            81,             4,           7.5,           764411904,       5.69531(GB)
//             3,             1,        8,        8,           768(12*4^3),            48,             8,            81,             4,          3.75,          6115295232,       45.5625(GB)
//             4,             1,        8,        8,          3072(12*4^4),            96,             8,            81,             4,         1.875,         48922361856,         364.5(GB)
//             5,             1,        8,        8,         12288(12*4^5),           192,             8,            81,             4,        0.9375,        391378894848,          2916(GB)
//             6,             1,        8,        8,         49152(12*4^6),           384,             8,            81,             4,       0.46875,       3131031158784,         23328(GB)
//             7,             1,        8,        8,        196608(12*4^7),           768,             8,            81,             4,      0.234375,      25048249270272,        186624(GB)
//             8,             1,        8,        8,        786432(12*4^8),          1536,             8,            81,             4,      0.117188,     200385994162176,   1.49299e+06(GB)
//             9,             1,        8,        8,       3145728(12*4^9),          3072,             8,            81,             4,     0.0585938,    1603087953297408,   1.19439e+07(GB)
//            10,             1,        8,        8,      12582912(12*4^10),         6144,             8,            81,             4,     0.0292969,   12824703626379264,   9.55515e+07(GB)
// NOTE : TILE size for direction can be 1,4,16,16*16,....for 'NEST healpix scheme'
// also need consider the opposite direction will share the same memory access for the 3D volume
class Exp_Mweight_old {
	class Sparse_Data;
	Sparse_Data * const sparse_data;

	int				_nr_images;
	int				_nr_classes;
	int				_nr_dir;
	int				_nr_psi;
	int				_capacity;
	std::string		current_fn;
	std::ofstream	os_file;
	double			max_nr_significant_ihidden;
	double			max_density;
	double			max_size;

public:
	int const & nr_images;
	int const & nr_classes;
	int const & nr_dir;
	int const & nr_psi;
	int const & capacity;

	Exp_Mweight_old();
	~Exp_Mweight_old();
	void init(int _nr_images, int _nr_classes, int _nr_dir, int _nr_psi, int _capacity);
	void fini();
	void clear_image(int iimage);
	void reset();

	std::vector< std::pair<int, double> > &wptr_sparse(int iimage, int iclass, int idir, int ipsi);

public:

	void analysis(std::string output_fn, std::string note);
		// do some analysis for exp_Mweight....(the compare matrix)

#ifndef DONT_INCLUDE_SAMPLING
	void printSamplingInfo();
#endif
};
