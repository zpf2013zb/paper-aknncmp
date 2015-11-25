#include <stdlib.h>
#include "rtree.h"

using namespace std;
#include <vector>

typedef vector<DATA*> DataStore;

struct BranchEvent {
	float mindist;	
	RTNode* pnode;	// null means current is a real data
	int son;
};

struct BranchComp {
	bool operator () (BranchEvent left,BranchEvent right) const
    { return left.mindist > right.mindist; }
};

typedef	priority_queue<BranchEvent,vector<BranchEvent>,BranchComp> BranchQueue;



float getMinDist(float *bounces1,float *bounces2) {
    float sum = 0.0;
    for(int i = 0; i < DIMENSION; i++) {
    	float diff=0;
		if ( bounces1[2*i+1] < bounces2[2*i])
		    diff = bounces2[2*i]-bounces1[2*i+1];
		else if ( bounces2[2*i+1] < bounces1[2*i])
		    diff = bounces1[2*i]-bounces2[2*i+1];
		sum += pow(diff,2);
    }
    return sqrtf(sum);
}

// for rect DATA !
float getMinDist(DATA* p,float *bounces) {
	return getMinDist(p->data,bounces);
}

// for rect DATA !
float getDataDist(DATA* p1,DATA* p2) {
	return getMinDist(p1->data,p2->data);
}

float getSumDataDist(DATA* p,DataStore& queries) {
	float sumVal=0;
	for (int s=0;s<queries.size();s++) {
		float curVal=getDataDist(p,queries[s]);
		if (isWeightUse) curVal=curVal*weights[s];
		sumVal+=curVal;
	}
	return sumVal;
}

float getMaxDataDist(DATA* p,DataStore& queries) {
	float maxVal=0;
	for (int s=0;s<queries.size();s++) {
		float curVal=getDataDist(p,queries[s]);
		if (isWeightUse) curVal=curVal*weights[s];
		maxVal=max(maxVal,curVal);
	}
	return maxVal;
}

float getNextNN(RTree* rt,BranchQueue& bQ,int& totna,DATA* point,DATA* NNpt) {
	while (!bQ.empty()) {
		BranchEvent event=bQ.top();
		bQ.pop();
		//printf("dist: %f\n",event.mindist);
		
		RTNode* node;
		if (event.pnode==NULL) 
			node=rt->root_ptr;	// assume root already loaded
		else {
			if (event.pnode->is_data_node()) {
				(*NNpt)=((RTDataNode*)event.pnode)->data[event.son];
				return event.mindist;
				continue;	// play safe, otherwise recursive looping
			} else
				node=((RTDirNode*)event.pnode)->entries[event.son].get_son();
		}
		totna++;
		if (node->is_data_node()) {
			RTDataNode* datanode=(RTDataNode *)node;
		    for (int i=0; i<datanode->num_entries; i++) {
		    	BranchEvent newevent;
		        newevent.pnode=node;	newevent.son=i;
		        newevent.mindist=getDataDist(point,&(datanode->data[i]));
		        bQ.push(newevent);
		    }
		} else {
			RTDirNode* dirnode=(RTDirNode *)node;
			for (int i=0; i<dirnode->num_entries; i++) {
				BranchEvent newevent;
				newevent.pnode=node;	newevent.son=i;
				newevent.mindist=getMinDist(point,dirnode->entries[i].bounces);
				bQ.push(newevent);
		    }
		}
	}
	return MAXREAL;	// no solution
}


float getNextAggNN(RTree* rt,BranchQueue& bQ,int& totna,DataStore& qrypts,DATA* NNpt) {
	float finaldist;
	RTNode* node;
	
	while (!bQ.empty()) {
		BranchEvent event=bQ.top();
		bQ.pop();		
		//printf("dist: %f\n",event.mindist);
		
		if (event.pnode==NULL) 
			node=rt->root_ptr;	// assume root already loaded
		else {
			if (event.pnode->is_data_node()) {
				(*NNpt)=((RTDataNode*)event.pnode)->data[event.son];
				// correct, already agg. dist
				return event.mindist;	
				continue;	// play safe, otherwise recursive looping
			} else 
				node=((RTDirNode*)event.pnode)->entries[event.son].get_son();
		}
		totna++;
		if (node->is_data_node()) {
			RTDataNode* datanode=(RTDataNode *)node;
		    for (int i=0; i<datanode->num_entries; i++) {
		    	BranchEvent newevent;
		        newevent.pnode=node;	newevent.son=i;
		        float finaldist=0;	// agg dist.
		    	for (int s=0;s<qrypts.size();s++) {
		    		float dist=getDataDist(qrypts[s],&(datanode->data[i]));
		    		if (isWeightUse) dist=dist*weights[s];	// wgt
		    		if (isSumDist) 
		    			finaldist+=dist;
		    		else 
		    			finaldist=max(dist,finaldist);
				}
		        newevent.mindist=finaldist;
		        bQ.push(newevent);
		    }
		} else {
			RTDirNode* dirnode=(RTDirNode *)node;
			for (int i=0; i<dirnode->num_entries; i++) {
				BranchEvent newevent;
				newevent.pnode=node;	newevent.son=i;
				finaldist=0;	// agg dist.
		    	for (int s=0;s<qrypts.size();s++) {
		    		float dist=getMinDist(qrypts[s],dirnode->entries[i].bounces);
		    		if (isWeightUse) dist=dist*weights[s];	// wgt
		    		if (isSumDist) 
		    			finaldist+=dist;
		    		else 
		    			finaldist=max(dist,finaldist);
		        }
				newevent.mindist=finaldist;
				bQ.push(newevent);
		    }
		}
	}
	return MAXREAL;	// no solution
}

/*void INN(RTree* rt,DATA* point,int& totna) {
	// need to be changed for rect rt !
	exit(0);
	BranchQueue bQ;
	
	{BranchEvent event;
	event.pnode=NULL;
	event.son=-1;
	event.mindist=0;	// always search from root node
	bQ.push(event);
	}
	
	DATA* NNpt=new DATA;
	printf("query: (%0.3f,%0.3f)\n",point->data[0],point->data[1]);
	for (int k=0;k<10;k++) {
		float dist=getNextNN(rt,bQ,totna,point,NNpt);
		printf("(%0.3f,%0.3f): %f\n",NNpt->data[0],NNpt->data[1],dist);
	}
	delete NNpt;
	
	printf("size: %d\n",bQ.size());
}*/

/*
float* getMBR(DataStore& queries) {
    float* mbr=new float[2*DIMENSION];
    
    // assume 2D
    mbr[0]=mbr[2]=MAXREAL;
    mbr[1]=mbr[3]=-MAXREAL;   
    for (int j=0; j < queries.size(); j++) {
		DATA* rect=queries[j];
		for (int i = 0; i < DIMENSION; i ++) {
			mbr[2*i]   = min(mbr[2*i],   rect->data[2*i]);
			mbr[2*i+1] = max(mbr[2*i+1], rect->data[2*i+1]);
       	}
    }
    return mbr;
}*/
