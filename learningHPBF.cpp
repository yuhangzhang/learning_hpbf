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
	tcomp->cid = _para.size()+1;
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
	int size = index.size();
		
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
	int size = index.size();

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



double learningHPBF::learn(const vector<bool>& y,vector<double> &w, const vector<double>& lambda, const vector<double> & slacks, int maxorder)
{
	_TermID.destroy();
	_PosiID.destroy();
	fillTermID();
	fillPosiID(maxorder);



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
