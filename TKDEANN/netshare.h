#ifndef __NET_SHARED
#define __NET_SHARED

#include <queue>

#define MAX_DIST 9999999

typedef FastList<int> IDset;

int NodeNum;
int NNnum=1;
float* weights;
bool isWeightUse=false;	// only for d_{max}
bool isSumDist=true;

struct QueryPoint {
	int Ni,Nj;
	float xc,yc,ndist,vP;	// fields xc,yc,ndist for genquery only
};

int gQrySize;
QueryPoint* gQryPts;

void ReadQueryFile(const char* queryf) {
	gQrySize=0;
	FILE* cquery=fopen(queryf,"r");
	CheckFile(cquery,queryf);
	fscanf(cquery,"%d\n",&gQrySize);
	gQryPts=new QueryPoint[gQrySize];
	
	//printf("s: %d\n",gQrySize);
	for (int s=0;s<gQrySize;s++) {
		fscanf(cquery,"%d %d %f\n",&(gQryPts[s].Ni),&(gQryPts[s].Nj),&(gQryPts[s].vP));
		//printf("%d %d %f\n",gQryPts[s].Ni,gQryPts[s].Nj,gQryPts[s].vP);
	}
	fclose(cquery);
}

struct edge {
	int FirstRow;
	int Ni,Nj;
	float dist;
	FastArray<float> pts;
	float dLB,dUB;		// for gendata use only
};

typedef map<int,edge*> EdgeMapType;

// build AdjList on fly
FastArray<int>* AdjList;
EdgeMapType EdgeMap;	// key: i0*NodeNum+j0

inline int getKey(int i,int j) {
	int i0=min(i,j),j0=max(i,j);	// be careful of the order in other places
	return (i0*NodeNum+j0);
}

inline void breakKey(int key,int& Ni,int& Nj) {
	Ni=key/NodeNum;	Nj=key%NodeNum;
}

void printEdgeKeys() {
	int Ni,Nj;
	printf("EdgeKeys:\n");
	EdgeMapType::iterator p=EdgeMap.begin();
	while (p!=EdgeMap.end()) {
		edge* e=p->second;		
		breakKey(p->first,Ni,Nj);
		printf("%d %d %d\n",Ni,Nj,(e==NULL));		
		p++;
	}
}


//-------
// Step
//-------
struct StepEvent {
	float dist;	
	int node,ClusID;	// ClusID for multiple expansion
	float accDist;		// posDist for gendata.cc only
	
	float gdist,hdist;
	float xc,yc;
	bool isSortAcc;
};

struct StepComparison {
	bool operator () (StepEvent left,StepEvent right) const
	{ return left.dist > right.dist; }	// to turn off this weighted feature for other exp. !!! 
//    {	// only for weighted cases !!!
//    	bool isSortAcc=left.isSortAcc;
//    	if (isWeightUse&&isSortAcc) {	// weighted expansion !
//    		return 
//    		(left.dist*weights[left.ClusID] > 
//    		right.dist*weights[right.ClusID]);
//    	} else 
//    		return (left.dist > right.dist);
//    }
    
};

typedef	priority_queue<StepEvent,vector<StepEvent>,StepComparison> StepQueue;

void printStepEvent(StepEvent& event)  {
	printf("(%d,%0.3f)\n",event.node,event.dist);
}


#include <math.h>
typedef map<int,float> DISTMAP;
DISTMAP elem_map;
struct DistElem {
	int id;
	float dist;
};
typedef vector<DistElem> DistQueue;


//typedef less<float> DescFloatComp;
//typedef priority_queue<float,vector<float>,DescFloatComp> DistQueue;

inline void InitDQ(DistQueue& dQ) {
	dQ.empty();
	elem_map.clear();

//	while (!dQ.empty()) dQ.pop();	// clear kNN-dist heap
}

inline void processDQ(DistQueue& dQ,float& topdist,float newdist,int id) {
	if (elem_map.count(id)>0) {
		if (elem_map[id]<=newdist) return;	// topdist unchanged
		
		// remove itself
		elem_map[id]=newdist;
		for (int i=0;i<dQ.size();i++)
			if (dQ[i].id==id) {
				dQ.erase(dQ.begin()+i);
				break;
			}
		
		//if (dQ.size()>=NNnum) topdist=dQ[NNnum-1].dist;	// still quite cheap !
	}
	
	DistElem new_elem;
	new_elem.id=id;		new_elem.dist=newdist;	
	if (dQ.size()<NNnum) {		// bestdist still infinity
		elem_map[id]=newdist;	// unsorted here !
		dQ.push_back(new_elem);		
		
		// must be sorted in order for the case below to be valid !
		int pos=dQ.size()-1;
		while (pos>0&&dQ[pos-1].dist>=newdist) pos--;
		for (int i=dQ.size()-1;i>pos;i--) dQ[i]=dQ[i-1];
		dQ[pos]=new_elem;
	} else {
		if (newdist<dQ[NNnum-1].dist) {
			int old_id=dQ[NNnum-1].id;
			elem_map.erase(old_id);
			elem_map[id]=newdist;
			
			int pos=NNnum-1;
			while (pos>0&&dQ[pos-1].dist>=newdist) pos--;
			for (int i=NNnum-1;i>pos;i--) dQ[i]=dQ[i-1];
			dQ[pos]=new_elem;			
		}
	}
	if (dQ.size()>=NNnum) topdist=dQ[NNnum-1].dist;	// still quite cheap !
	
//	if (dQ.size()<NNnum) 
//		dQ.push(newdist);	// bestdist still infinity
//	else {
//		if (newdist<dQ.top()) {
//			dQ.push(newdist);
//			dQ.pop();
//		}
//	}
//	if (dQ.size()>=NNnum) topdist=dQ.top();	// still quite cheap !
}

#endif // __NET_SHARED


