#include <vector>
#include <algorithm>
#include "./resmap_error.h"
#include "mkl.h"

void findBestPointsUnitTestUnitTest();

class FindBestPoints {
public:
	FindBestPoints();
	virtual ~FindBestPoints();

	void setWeightCount(int weight_count, double target_weight, double sum_weight);

	int get_weight_count() const;

	double* get_sorted_weights_ptr_if_avail();

	void insert(double w);

	void prepareToPop();

	double pop();

	struct Impl;
private:
	Impl* _impl;
};


