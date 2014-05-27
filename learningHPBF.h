#include "../poly/poly.h"
#include <engine.h>



class learningHPBF
{
public:

	typedef struct component
	{
		poly<double>* A_i;
		int cid;
		component* next;
	}component;

	learningHPBF();
	int addComponent(poly<double>* A_i); //return the id of the added component

	//void remove_cQPBF(int cid);//remove component-cid
	
	double learn(const vector<bool>& y,vector<double> &w, const vector<double>& lambda, const vector<double> & slacks, int maxorder);
	
	int lastVar();
	int lastComp();
	double optimize(vector<bool>& y);
	void setPara(const vector<double>& p);
	void getPara(vector<double>& p);

	poly<int> * getTermID();
	int fillTermID();

	poly<int> * getPosiID();
	int fillPosiID(int maxorder);

private:
	vector<double> _para;
	component* _component;
	int _lastvar;

	poly<int> _TermID;
	poly<int> _PosiID;

	int fill(int start, int end, vector<int>& index, int& counter);
	int dualfill(int start, int end, vector<int>& index, int& counter, int maxorder);
};



