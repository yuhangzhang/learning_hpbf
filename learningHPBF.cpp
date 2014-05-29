#include "learningHPBF.h"


learningHPBF::learningHPBF()
{
	_lastvar = 0;
	_component = NULL;
	_para.resize(0);

	return;
}

int learningHPBF::addComponent(poly<double>* A_i)
{
	component *tcomp = new component;

	tcomp->A_i = A_i;
	tcomp->cid = int(_para.size())+1;
	tcomp->next = _component;
	_component = tcomp;

	_para.resize(_para.size()+1,0);

	if(_lastvar<A_i->lastVar()) _lastvar = A_i->lastVar();

	return tcomp->cid;
}


int learningHPBF::lastVar()
{
	return _lastvar;
}

int learningHPBF::lastComp()
{
	if(_component == NULL) return -1;
	else return _component->cid;
}


void learningHPBF::setPara(const vector<double>& p)
{
	_para = p;
	return;
}

void learningHPBF::getPara(vector<double>& p)
{
	p = _para;
	return;
}


int learningHPBF::fillTermID()
{	
	int counter=1;
	vector<int> index;
	index.resize(0);
	
	return fill(1,_lastvar,index,counter)-1;
}

int learningHPBF::fill(int start, int end, vector<int>& index, int& counter)
{
	int size = int(index.size());
		
	for(int i=start;i<=end;i++)
	{
		index.resize(size+1);
		index[size] = i;
		_TermID.addTerm(index,counter++);
		fill(i+1,_lastvar,index,counter);
	}



	return counter;
}


int learningHPBF::fillPosiID(int maxorder)
{	
	int counter=1;
	vector<int> index;
	index.resize(0);
	
	return dualfill(1,_lastvar,index,counter,maxorder)-1;
}

int learningHPBF::dualfill(int start, int end, vector<int>& index, int& counter, int maxorder)
{
	int size = int(index.size());

	if(size<maxorder)
	{		
		for(int i=start;i<=end;i++)
		{
			index.resize(size+1);
			index[size] = i;
			_PosiID.addTerm(index,counter++);
			dualfill(i+1,_lastvar,index,counter,maxorder);

			index.resize(size+1);
			index[size] = -i;
			_PosiID.addTerm(index,counter++);
			dualfill(i+1,_lastvar,index,counter,maxorder);
		}
	}
	return counter;
}



double learningHPBF::learn(const vector<bool>& y,vector<double> &w, const double& lambda, vector<double> & slacks, int maxorder)
{
	_TermID.destroy();
	_PosiID.destroy();
	fillTermID();
	fillPosiID(maxorder);

	std::vector< vector<double> > tripletList;
	tripletList.reserve(_PosiID.numTerm()*10);//a rough estimation of nonzero entries in the linear constraints
	
	printf("Construct Matrix A ...\n");

	for(component *k=_component;k!=NULL;k=k->next)
	{
		if(k->A_i==NULL) continue;
		
		for(poly<double>::TERMS::iterator it=k->A_i->firstTerm();it!=k->A_i->lastTerm();it++)
		{
			vector<double> index(3);
				
			index[0] = _TermID.getTerm(it->first);//row
			index[1] = k->cid;//column
			index[2] = it->second;


			tripletList.push_back(index);//remember we count from 1
		}
	}

	for(poly<int>::TERMS::iterator it=_PosiID.firstTerm();it!=_PosiID.lastTerm();it++)
	{
		vector<double> index(3);
		index[1] = lastComp()+it->second;//column

		poly<double> cols;
		posi2poly(it->first, cols, 1);
		
		for(poly<double>::TERMS::iterator it2=cols.firstTerm();it2!=cols.lastTerm();it2++)
		{
			index[0] = _TermID.getTerm(it2->first);
			index[2] = -it2->second;
			tripletList.push_back(index);
		}

		cols.destroy();

		if(ifZero(it->first,y)==1)
		{
			index[0] = _TermID.numTerm()+1;
			index[2] = 1;
			tripletList.push_back(index);
		}
	}

	for(int i=1;i<=_lastvar;i++)
	{
		vector<double> index(3);
		vector<int> index2(1,i);
		index[0] = _TermID.getTerm(index2);
		index[1] = lastComp()+_PosiID.numTerm()+i;
		index[2] = y[i-1]*-2+1;
		tripletList.push_back(index);
	}


	printf("passing to matlab\n");

	Engine *eg = engOpen(NULL);
	engEvalString(eg,"clear all");
	engEvalString(eg,"cd C:\\projects\\cvx;");
	engEvalString(eg,"cvx_setup;");		



	mxArray *scalar = mxCreateDoubleMatrix(1,1,mxREAL);

	*((double *) mxGetPr(scalar))=tripletList.size();
	engPutVariable(eg,"nzentry",scalar);
	*((double *) mxGetPr(scalar))=_TermID.numTerm();
	engPutVariable(eg,"numpoly",scalar);
	*((double *) mxGetPr(scalar))=_PosiID.numTerm();
	engPutVariable(eg,"numposi",scalar);
	*((double *) mxGetPr(scalar))=_para.size();
	engPutVariable(eg,"numcom",scalar);
	*((double *) mxGetPr(scalar))=_lastvar;
	engPutVariable(eg,"numvar",scalar);
	*((double *) mxGetPr(scalar))=lambda;
	engPutVariable(eg,"lambda",scalar);
	

	mxArray *index_i = mxCreateDoubleMatrix(tripletList.size(),1,mxREAL);
	mxArray *index_j = mxCreateDoubleMatrix(tripletList.size(),1,mxREAL);
	mxArray *value_s = mxCreateDoubleMatrix(tripletList.size(),1,mxREAL);

	for(unsigned int i=0;i<tripletList.size();i++)
	{
		mxGetPr(index_i)[i] = tripletList[i][0];
		mxGetPr(index_j)[i] = tripletList[i][1];
		mxGetPr(value_s)[i] = tripletList[i][2];
		tripletList[i].clear();
	}

	tripletList.clear();

	engPutVariable(eg,"index_i",index_i);
	engPutVariable(eg,"index_j",index_j);
	engPutVariable(eg,"value_s",value_s);

	mxDestroyArray(index_i);
	mxDestroyArray(index_j);
	mxDestroyArray(value_s);

	engEvalString(eg,"A = sparse(index_i,index_j,value_s,numpoly+1,numcom+numposi+numvar,nzentry);");
printf("cvx\n");//getchar();

	engEvalString(eg,"cvx_begin");

	engEvalString(eg,"variable d(numvar);");
	engEvalString(eg,"variable w(numcom);");
	engEvalString(eg,"variable p(numposi);");



	engEvalString(eg,"minimize(norm(w,1)+sum(d));");
	engEvalString(eg,"subject to");
	engEvalString(eg,"A*[w;p;d]==zeros(numpoly+1,1);");
	engEvalString(eg,"w(1)==1;");
	engEvalString(eg,"d>=zeros(numvar,1);");
	engEvalString(eg,"p>=zeros(numposi,1);");
	engEvalString(eg,"cvx_end");

printf("cvx end\n");


	mxArray *slacksmx = engGetVariable(eg,"d");
	mxArray *weightsmx = engGetVariable(eg,"w");
	mxArray *optvalmx = engGetVariable(eg,"cvx_optval");
	for(int i=0;i<_para.size();i++)
	{
		_para[i] = mxGetPr(weightsmx)[i];
	}

	for(int i=0;i<_lastvar;i++)
	{
		//printf("i=%d\n",i);
		slacks[i] = mxGetPr(slacksmx)[i];
	}

	
	double optval = mxGetPr(optvalmx)[0];

	return optval;



	return 0;

}

poly<int> * learningHPBF::getTermID()
{
	return &_TermID;
}

poly<int> * learningHPBF::getPosiID()
{
	return &_PosiID;
}

void learningHPBF::posi2poly(vector<int> posi,poly<double>& index, double coeff)
{
	for(int i=0;i<posi.size();i++)
	{
		if(posi[i]<0)
		{
			posi[i] = -posi[i];
			posi2poly(posi,index,-coeff);
			posi.erase(posi.begin()+i);
			posi2poly(posi,index,coeff);
			break;
		}
		else if(i == posi.size()-1) 
		{
			index.addTerm(posi,coeff);
		}
	}

	return;
}

bool learningHPBF::ifZero(const vector<int>& posi,const vector<bool>& y)
{
	bool value=1;

	for ( vector<int>::const_iterator it = posi.begin() ; it != posi.end(); it++)
	{
		if(*it<0&&y[-*it]==1) 
		{
			value = 0;
			break;
		}
		else if(*it>0&&y[*it]==0)
		{
			value = 0;
			break;
		}
	}

	return value;
}

