#include "learningHPBF.h"

int main()
{

	learningHPBF learner;

	poly<double> p[5];

	vector<int> index;


	index.resize(1);
	index[0]=1;
	p[0].addTerm(index,-1);
	index[0]=2;
	p[0].addTerm(index,1);
	index[0]=3;
	p[0].addTerm(index,1);

	index.resize(1);
	index[0]=1;
	p[1].addTerm(index,1);
	index[0]=2;
	p[1].addTerm(index,2);
	index[0]=3;
	p[1].addTerm(index,2);

	index.resize(2);
	index[0] = 1; index[1] = 2;
	p[2].addTerm(index,1);
	index[0] = 1; index[1] = 3;
	p[2].addTerm(index,1);
	index[0] = 2; index[1] = 3;
	p[2].addTerm(index,1);

	index.resize(2);
	index[0] = 1; index[1] = 2;
	p[3].addTerm(index,1);
	index[0] = 1; index[1] = 3;
	p[3].addTerm(index,2);
	index[0] = 2; index[1] = 3;
	p[3].addTerm(index,2);

	index.resize(1);
	index[0]=1;
	p[4].addTerm(index,10);
	index[0]=2;
	p[4].addTerm(index,10);
	index[0]=3;
	p[4].addTerm(index,10);

	learner.addComponent(&p[0]);
	learner.addComponent(&p[1]);
	learner.addComponent(&p[2]);
	learner.addComponent(&p[3]);
	learner.addComponent(&p[4]);

	vector<bool> y(3);
	y[0]=0;
	y[1]=1;
	y[2]=1;

	vector<double> w(3);
	vector<double>  slacks(3);

	learner.learn(y,w, 10, slacks, 3);
}