#include "learningHPBF.h"


int main()
{
	learningHPBF learner;

	poly<double> p;

	std::vector<int> v1(3,3);

	p.addTerm(v1,1);

	learner.addComponent(&p);

	printf("last var =  %d %d\n",p.lastVar(),learner.lastVar());

	int total=learner.fillTermID();

	poly<int>::TERMS::iterator it = learner.getTermID()->firstTerm();
	poly<int>::TERMS::iterator it2 = learner.getTermID()->lastTerm();
		
	for(;it!=it2;it++)
	{
		vector<int> tv = it->first;
		printf("< ");
		for(int i=0;i<tv.size();i++)
		{
			printf("%d ",tv[i]);
		}
		printf("> ");

		printf("%d\n",it->second);
	}

	printf("total=%d\n",total);

	total=learner.fillPosiID(2);

	it = learner.getPosiID()->firstTerm();
	it2 = learner.getPosiID()->lastTerm();
		
	for(;it!=it2;it++)
	{
		vector<int> tv = it->first;
		printf("< ");
		for(int i=0;i<tv.size();i++)
		{
			printf("%d ",tv[i]);
		}
		printf("> ");

		printf("%d\n",it->second);
	}

	printf("total=%d\n",total);

	printf("=============================================\n",total);

	poly<double> temp;
	vector<int> posi(6);

	posi[0] = 1;
	posi[1] = -2;
	posi[2] = -3;
	posi[3] = -4;
	posi[4] = -5;
	posi[5] = -6;

	learner.posi2poly( posi, temp,  1);

	poly<double>::TERMS::iterator it3 = temp.firstTerm();
	poly<double>::TERMS::iterator it4 = temp.lastTerm();
		
	for(;it3!=it4;it3++)
	{
		vector<int> tv = it3->first;
		printf("< ");
		for(int i=0;i<tv.size();i++)
		{
			printf("%d ",tv[i]);
		}
		printf("> ");

		printf("%f\n",it3->second);
	}


	return 0;
}