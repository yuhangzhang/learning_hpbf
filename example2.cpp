#include "learningHPBF.h"

int main(int argc, char* argv[])
{
	//8 variable 00001111

	learningHPBF learner;

	//6 basis functions
	poly<double> p[8];
	poly<double> mg;

	vector<int> index;

	srand(17);

	index.resize(1);
	index[0]=1;mg.addTerm(index,-1);
	index[0]=2;mg.addTerm(index,-1);
	index[0]=3;mg.addTerm(index,-1);
	index[0]=4;mg.addTerm(index,-1);
	index[0]=5;mg.addTerm(index,1);
	index[0]=6;mg.addTerm(index,1);
	index[0]=7;mg.addTerm(index,1);
	index[0]=8;mg.addTerm(index,1);

	index.resize(1);
	for(int i=1;i<=8;i++)
	{
		index[0]=i;
		p[0].addTerm(index,double(rand()%200)/100-1);
	}


	index.resize(1);
	for(int i=1;i<=8;i++)
	{
		index[0]=i;
		p[1].addTerm(index,double(rand()%200)/100-1);
	}


	index.resize(2);
	for(int i=1;i<=8;i++)
	{
		for(int j=i-2;j<=i+2;j++)
		{
			if(j<1)
			{
				index[0]=i;
				index[1]=j+8;
			}
			else if(j>8)
			{
				index[0]=j-8;
				index[1]=i;
			}
			else if(j==i)
			{
				continue;
			}
			else if(i<j)
			{
				index[0]=i;
				index[1]=j;
			}
			else
			{
				index[0]=j;
				index[1]=i;
			}

			p[2].addTerm(index,double(rand()%200)/100-1);
		}		
	}


	index.resize(2);
	for(int i=1;i<=8;i++)
	{
		for(int j=i-2;j<=i+2;j++)
		{
			if(j<1)
			{
				index[0]=i;
				index[1]=j+8;
			}
			else if(j>8)
			{
				index[0]=j-8;
				index[1]=i;
			}
			else if(j==i)
			{
				continue;
			}
			else if(i<j)
			{
				index[0]=i;
				index[1]=j;
			}
			else
			{
				index[0]=j;
				index[1]=i;
			}

			p[3].addTerm(index,double(rand()%200)/100-1);
		}		
	}

	index.resize(2);
	for(int i=1;i<=8;i++)
	{
		for(int j=i-2;j<=i+2;j++)
		{
			if(j<1)
			{
				index[0]=i;
				index[1]=j+8;
			}
			else if(j>8)
			{
				index[0]=j-8;
				index[1]=i;
			}
			else if(j==i)
			{
				continue;
			}
			else if(i<j)
			{
				index[0]=i;
				index[1]=j;
			}
			else
			{
				index[0]=j;
				index[1]=i;
			}

			p[4].addTerm(index,double(rand()%200)/100-1);
		}		
	}

	index.resize(2);
	for(int i=1;i<=8;i++)
	{
		for(int j=i-2;j<=i+2;j++)
		{
			if(j<1)
			{
				index[0]=i;
				index[1]=j+8;
			}
			else if(j>8)
			{
				index[0]=j-8;
				index[1]=i;
			}
			else if(j==i)
			{
				continue;
			}
			else if(i<j)
			{
				index[0]=i;
				index[1]=j;
			}
			else
			{
				index[0]=j;
				index[1]=i;
			}

			p[5].addTerm(index,double(rand()%200)/100-1);
		}		
	}

	index.resize(3);
	index[0]=1;index[1]=2;index[2]=3;p[6].addTerm(index,double(rand()%200)/100-1);
	index[0]=2;index[1]=3;index[2]=4;p[6].addTerm(index,double(rand()%200)/100-1);
	index[0]=3;index[1]=4;index[2]=5;p[6].addTerm(index,double(rand()%200)/100-1);
	index[0]=4;index[1]=5;index[2]=6;p[6].addTerm(index,double(rand()%200)/100-1);
	index[0]=5;index[1]=6;index[2]=7;p[6].addTerm(index,double(rand()%200)/100-1);
	index[0]=6;index[1]=7;index[2]=8;p[6].addTerm(index,double(rand()%200)/100-1);
	index[0]=1;index[1]=7;index[2]=8;p[6].addTerm(index,double(rand()%200)/100-1);
	index[0]=1;index[1]=2;index[2]=8;p[6].addTerm(index,double(rand()%200)/100-1);

	index.resize(3);
	index[0]=1;index[1]=2;index[2]=3;p[7].addTerm(index,double(rand()%200)/100-1);
	index[0]=2;index[1]=3;index[2]=4;p[7].addTerm(index,double(rand()%200)/100-1);
	index[0]=3;index[1]=4;index[2]=5;p[7].addTerm(index,double(rand()%200)/100-1);
	index[0]=4;index[1]=5;index[2]=6;p[7].addTerm(index,double(rand()%200)/100-1);
	index[0]=5;index[1]=6;index[2]=7;p[7].addTerm(index,double(rand()%200)/100-1);
	index[0]=6;index[1]=7;index[2]=8;p[7].addTerm(index,double(rand()%200)/100-1);
	index[0]=1;index[1]=7;index[2]=8;p[7].addTerm(index,double(rand()%200)/100-1);
	index[0]=1;index[1]=2;index[2]=8;p[7].addTerm(index,double(rand()%200)/100-1);

	//printf("%f\n",p[5].getTerm(index));getchar();

	learner.addComponent(&mg);
	learner.addComponent(&p[0]);
	learner.addComponent(&p[1]);
	learner.addComponent(&p[2]);
	learner.addComponent(&p[3]);
	learner.addComponent(&p[4]);
	learner.addComponent(&p[5]);
	learner.addComponent(&p[6]);
	learner.addComponent(&p[7]);

	vector<double> w(9);
	vector<double>  slacks(8);

	vector<bool> y(8);
	y[0]=false;
	y[1]=false;
	y[2]=false;
	y[3]=false;
	y[4]=true;
	y[5]=true;
	y[6]=true;
	y[7]=true;

	double optval=learner.learn(y,w, 10, slacks, 4);

	FILE *fpw=fopen(argv[1],"wt");

	fprintf(fpw,"-1 %f \n",optval);

	for(int i=0;i<256;i++)
	{
		y[0] = (i%2>0)? 1:0;
		y[1] = (i%4>1)? 1:0;
		y[2] = (i%8>3)? 1:0;
		y[3] = (i%16>7)? 1:0;
		y[4] = (i%32>15)? 1:0;
		y[5] = (i%64>31)? 1:0;
		y[6] = (i%128>63)? 1:0;
		y[7] = (i%256>127)? 1:0;

		double value = learner.evaluate2(y);

		for(int j=0;j<8;j++)
		{
			fprintf(fpw,"%d ",int(y[j]));
			//printf("%d ",int(y[j]));
		}
		fprintf(fpw,"%d %f \n",i,value);
		//printf("%d %f \n",i,value);getchar();
	}

	fclose(fpw);

	return 0;

}


