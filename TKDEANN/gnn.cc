#include "utility.h"
#include "netshare.h"
#include "diskbased.h"
#include "rtree.h"

typedef map<int,float> DISTMAP;
DISTMAP* DistMaps;	// only use cheap storage for range search
BitStore* astar_Visits;		// tmp. bitmap for A* search

typedef vector<int> 	IntVec;
typedef vector<float> 	FloatVec;

FloatVec M_xcrd,M_ycrd,compdist,bestcompdist;

BitStore onSameEdge;
float** partdists;
float bestdist;	// best value found so far

bool isNNquery=true;
bool isEuclidSpace=true;
char ptgrprtfn[255];

#include "rtacc.h"
DistQueue dQ_graph;

float getDist(float x1,float y1,float x2,float y2) {
	return float(sqrt((x2-x1)*(x2-x1)+(y2-y1)*(y2-y1)));
}

inline float getNetworkDist(float vJ,BitStore& onSameEdge,float eKdist) {
	float dist=0;
	for (int sp=0;sp<gQrySize;sp++) {
		float distL=partdists[sp][0]+vJ;
		float distR=partdists[sp][1]+eKdist-vJ;
		float distT=min(distL,distR);
		if (onSameEdge[sp]) {
			float distQ=ValAbs(gQryPts[sp].vP-vJ);
			if (distQ<distT) distT=distQ;
		}
		
		compdist[sp]=distT;	// for log (log the original one) !
		if (isWeightUse) distT=distT*weights[sp];
		
		if (isSumDist)
			dist+=distT;
		else 
			dist=max(distT,dist);
	}
	return dist;
}

inline void UpdateOnePt(float vJ,bool& hasChanged,BitStore& onSameEdge,float eKdist,int ptid) {
	float dist=getNetworkDist(vJ,onSameEdge,eKdist);
	if (dist<bestdist) {	// for S-dist
		hasChanged=true;
		processDQ(dQ_graph,bestdist,dist,ptid);
		bestcompdist.assign(compdist.begin(),compdist.end());
	}
}

inline float getOptPos(int NodeID,int NewNodeID, float eKdist,float& bestPos) {
	float vJ,bestVal,tmpPos,tmpVal;
	int sid=0,eid=1;
	if (NewNodeID<NodeID) { sid=1; eid=0;}
	
	for (int sp=0;sp<gQrySize;sp++) {
		partdists[sp][0]=partdists[sp][1]=MAX_DIST;
		if (DistMaps[sp].count(NodeID)>0) 
			partdists[sp][sid]=DistMaps[sp][NodeID];
		if (DistMaps[sp].count(NewNodeID)>0) 
			partdists[sp][eid]=DistMaps[sp][NewNodeID];
		int qNi=gQryPts[sp].Ni,qNj=gQryPts[sp].Nj;
		onSameEdge[sp]=(NodeID==qNi&&NewNodeID==qNj)||(NodeID==qNj&&NewNodeID==qNi);
	}
	
	// for isSumDist only
	bestPos=0;
	bestVal=getNetworkDist(bestPos,onSameEdge,eKdist);

	tmpPos=eKdist;
	tmpVal=getNetworkDist(tmpPos,onSameEdge,eKdist);
	if (tmpVal<bestVal) { bestPos=tmpPos; bestVal=tmpVal; }
	
	int NUM_STEPS=100;
	if (!isWeightUse) {
		if (isSumDist) {
			for (int sp=0;sp<gQrySize;sp++) {
				if (onSameEdge[sp]) {
					tmpPos=gQryPts[sp].vP;
					tmpVal=getNetworkDist(tmpPos,onSameEdge,eKdist);
					if (tmpVal<bestVal) { bestPos=tmpPos; bestVal=tmpVal; }
				}
			}
		} else {
			// can find a loose lb for ...
			//bestVal-=eKdist;
			for (int j=0;j<NUM_STEPS;j++) {	// scan
				tmpPos=j*eKdist/(NUM_STEPS-1);
				tmpVal=getNetworkDist(tmpPos,onSameEdge,eKdist);
				if (tmpVal<bestVal) { bestPos=tmpPos; bestVal=tmpVal; }
			}
		}
	} else {
		for (int j=0;j<NUM_STEPS;j++) {	// scan
			tmpPos=j*eKdist/(NUM_STEPS-1);
			tmpVal=getNetworkDist(tmpPos,onSameEdge,eKdist);
			if (tmpVal<bestVal) { bestPos=tmpPos; bestVal=tmpVal; }
		}
	}
	return bestVal;
}

bool UpdatePoints(int NodeID,int NewNodeID,int PtGrpKey, float eKdist) {
	float vJ;
	bool hasChanged=false;
	
	int sid=0,eid=1;
	if (NewNodeID<NodeID) { sid=1; eid=0;}
	
	for (int sp=0;sp<gQrySize;sp++) {
		partdists[sp][0]=partdists[sp][1]=MAX_DIST;
		if (DistMaps[sp].count(NodeID)>0) 
			partdists[sp][sid]=DistMaps[sp][NodeID];
		if (DistMaps[sp].count(NewNodeID)>0) 
			partdists[sp][eid]=DistMaps[sp][NewNodeID];
		int qNi=gQryPts[sp].Ni,qNj=gQryPts[sp].Nj;
		onSameEdge[sp]=(NodeID==qNi&&NewNodeID==qNj)||(NodeID==qNj&&NewNodeID==qNi);
	}
	
	if (isNNquery) {
		int PtGrpAddr,PtGrpSize,dummy;
		PtGrpAddr=pointQuery(PtTree,PtGrpKey,dummy);
		getFixedF(SIZE_P,Ref(PtGrpSize),PtGrpAddr);
	 	for (int j=0;j<PtGrpSize;j++) {	// scan
			getVarE(PT_P,Ref(vJ),PtGrpAddr,j);
			UpdateOnePt(vJ,hasChanged,onSameEdge,eKdist,PtGrpKey+j);
		}
	} else {
		vJ=0;
		getOptPos(NodeID,NewNodeID,eKdist,vJ);
		UpdateOnePt(vJ,hasChanged,onSameEdge,eKdist,-1);
	}
	//if (bestdist<MAX_DIST) {PrintElapsed();	exit(0);}
	/*if (hasChanged) 
		printf("c: %d %d %f\n",NodeID,NewNodeID,eKdist);*/
	return hasChanged;
}

/* // need to be changed for weighted case !!!
inline float getMinOptPos(int NodeID,int NewNodeID, float eKdist) {
	float dist=0,mincpdist,ldist,rdist;
	for (int sp=0;sp<gQrySize;sp++) {
		mincpdist=MAX_DIST;
		ldist=rdist=MAX_DIST;
		if (DistMaps[sp].count(NodeID)>0) ldist=DistMaps[sp][NodeID];
		if (DistMaps[sp].count(NewNodeID)>0) rdist=DistMaps[sp][NewNodeID];
		if (ldist<MAX_DIST) {
			if (rdist<MAX_DIST) 
				mincpdist=min(ldist,rdist);
			else 
				mincpdist=max(0,ldist-eKdist);
		} else {
			if (rdist<MAX_DIST) 
				mincpdist=max(0,rdist-eKdist);
			else 
				mincpdist=0;
		}
		if (onSameEdge[sp]) mincpdist=0;
		mincpdist=max(0,mincpdist);
		
		if (isSumDist)
			dist+=mincpdist;
		else 
			dist=max(mincpdist,dist);
	}
	return dist;
}*/

struct RgnType {
	float partial_dist;
	int crash;
};

typedef map<int,RgnType> REGIONMAP;

/*void ConcurrentExpansion(StepQueue& sQ) {
	int MaxHeapSize=0,AggPopSize=0;
	REGIONMAP epsNbrList;
	BitStore epsBits;
	int AdjListSize,NewNodeID,AdjGrpAddr,PtGrpKey;
	float eKdist,xVal,yVal;
	
	// initialize variables
	if (!isEuclidSpace) return;
	epsBits.assign(NodeNum,false);
	int maxNbrSize=0;
	while (!sQ.empty()) {
		if (MaxHeapSize<sQ.size()) MaxHeapSize=sQ.size();
		StepEvent event=sQ.top();
		sQ.pop();	AggPopSize++;
		
		int NodeID=event.node;
		int id=event.ClusID;
		
		if (DistMaps[id].count(NodeID)>0) continue;	// already found
		DistMaps[id][NodeID]=event.dist;
		if (event.dist>=bestdist) break;	// for both S-dist and M-dist: bound almost correct 
		
		AdjGrpAddr=getAdjListGrpAddr(NodeID);
		getFixedF(SIZE_A,Ref(AdjListSize),AdjGrpAddr);	// read # entries
		getFixedF(XCRD_A,Ref(xVal)	,AdjGrpAddr);
		getFixedF(YCRD_A,Ref(yVal)	,AdjGrpAddr);
		
		maxNbrSize=max(maxNbrSize,epsNbrList.size());
		float epsilon=(isSumDist)?(bestdist/gQrySize):(bestdist);
		if (epsilon>event.dist) { 	// any bug ?
			if (epsBits[NodeID]==false) {
				epsBits[NodeID]=true;
				epsNbrList[NodeID].crash=0;
				epsNbrList[NodeID].partial_dist=0;
			}
		}
		
		if (epsNbrList.count(NodeID)>0) {
			RgnType& rgt=epsNbrList[NodeID];
			rgt.crash+=1;	// update crash
			bool willErase=(rgt.crash>=gQrySize);
			
			if (isSumDist) {	// update cur dist and check if need prune
				epsNbrList[NodeID].partial_dist+=event.dist;
				if (!willErase) willErase=((gQrySize-rgt.crash)*event.dist+rgt.partial_dist>=bestdist);
			} else {	// update cur dist and check if need prune
				epsNbrList[NodeID].partial_dist=max(event.dist,epsNbrList[NodeID].partial_dist);
				if (!willErase) willErase=(rgt.partial_dist>=bestdist);	// current lb
			}
			
			// not useful at all
//			if (!willErase) { // correct if maxedgedist used ??
//				float curXdist=rgt.partial_dist;
//				for (int s=0;s<gQrySize;s++) 
//					if (rgt.bset[s]==false) {
//						float tmpdist=getDist(xVal,yVal,M_xcrd[s],M_ycrd[s]);
//						tmpdist=max(tmpdist,event.dist);
//						if (isSumDist) 
//							curXdist+=tmpdist;
//						else 
//							curXdist=max(tmpdist,curXdist);
//					}
//				if (curXdist>=bestdist) willErase=false;
//			}
			if (willErase) epsNbrList.erase(NodeID);
		}
		
		if (bestdist<MAX_DIST) {	// at least one found
			if (epsNbrList.size()==0) break;	// all pruned => exit
		}
		
		bool needPruneEps=false;
		for (int z=0;z<AdjListSize;z++) {
			getVarE(ADJNODE_A	,Ref(NewNodeID)	,AdjGrpAddr,z);
			getVarE(DIST_A		,Ref(eKdist)	,AdjGrpAddr,z);
			
			StepEvent newevent=event;	// copy ...
			newevent.node=NewNodeID;
			newevent.dist=event.dist+eKdist;
			if (DistMaps[id].count(NewNodeID)==0) sQ.push(newevent);	// propagation
			
			// if error, add cond. to (... && not on same edge) [sp]
			bool needUpdate=true;
			for (int sp=0;sp<gQrySize;sp++) 
				if (DistMaps[sp].count(NodeID)==0&&DistMaps[sp].count(NewNodeID)==0)
					needUpdate=false;
			if (needUpdate&&isNNquery) {
				float tmpposition;
				float bestVal=getOptPos(NodeID,NewNodeID,eKdist,tmpposition);
				if (bestVal>=bestdist) needUpdate=false;
			}
			if (!needUpdate) continue;
			
			getVarE(PTKEY_A	,Ref(PtGrpKey)	,AdjGrpAddr,z);
			if (PtGrpKey>=0||!isNNquery) {	// valid PtGrpKey  (-1 for invalid key)
				bool hasChanged=UpdatePoints(NodeID,NewNodeID,PtGrpKey,eKdist);
			}
		}
	}
	printf("Heap: max_size %d, pop_size %d\n",MaxHeapSize,AggPopSize);
	printf("bestdist: %f\n",bestdist);
	if (!sQ.empty()) printf("topdist: %f\n",sQ.top().dist);
	int totalmapsize=0;
	for (int sp=0;sp<gQrySize;sp++) {
		//printf("%d) td_sz: %d\n",sp,DistMaps[sp].size());
		totalmapsize+=DistMaps[sp].size();
	}
	printf("totalmapsize: %d, maxNbrSz: %d\n",totalmapsize,maxNbrSize);
}*/

// memory requirement: r bitmap and r tempdiststore
void ConcurrentExpansion(StepQueue& sQ) {
	int MaxHeapSize=0,AggPopSize=0;
	REGIONMAP epsNbrList;
	BitStore epsBits;
	int AdjListSize,NewNodeID,AdjGrpAddr,PtGrpKey;
	float eKdist,xVal,yVal;
	
	// initialize variables
	if (!isEuclidSpace) return;
	epsBits.assign(NodeNum,false);
	int maxNbrSize=0;
	while (!sQ.empty()) {
		if (MaxHeapSize<sQ.size()) MaxHeapSize=sQ.size();
		StepEvent event=sQ.top();
		sQ.pop();	AggPopSize++;
		
		int NodeID=event.node;
		int id=event.ClusID;
		
		if (DistMaps[id].count(NodeID)>0) continue;	// already found
		DistMaps[id][NodeID]=event.dist;
		
		float topbound=event.dist;	// change here !
		if (isWeightUse) topbound=topbound*weights[id];
		if (topbound>=bestdist) continue;	// for both S-dist and M-dist: bound almost correct 
		
		AdjGrpAddr=getAdjListGrpAddr(NodeID);
		getFixedF(SIZE_A,Ref(AdjListSize),AdjGrpAddr);	// read # entries
		getFixedF(XCRD_A,Ref(xVal)	,AdjGrpAddr);
		getFixedF(YCRD_A,Ref(yVal)	,AdjGrpAddr);
		
		maxNbrSize=max(maxNbrSize,epsNbrList.size());
		float epsilon=(isSumDist)?(bestdist/gQrySize):(bestdist);
		if (epsilon>topbound) { 	// any bug ?
			if (epsBits[NodeID]==false) {
				epsBits[NodeID]=true;
				epsNbrList[NodeID].crash=0;
				epsNbrList[NodeID].partial_dist=0;
			}
		}
		
		if (epsNbrList.count(NodeID)>0) {
			RgnType& rgt=epsNbrList[NodeID];
			rgt.crash+=1;	// update crash
			bool willErase=(rgt.crash>=gQrySize);
			
			float curcpdist=event.dist;
			if (isWeightUse) curcpdist=curcpdist*weights[id];
			if (isSumDist) {	// update cur dist and check if need prune
				rgt.partial_dist+=curcpdist;
//				if (!willErase&&!isWeightUse) 
//					willErase=((gQrySize-rgt.crash)*event.dist+rgt.partial_dist>=bestdist);
			} else {	// update cur dist and check if need prune
				rgt.partial_dist=max(curcpdist,rgt.partial_dist);
//				if (!willErase&&!isWeightUse) 
//					willErase=(rgt.partial_dist>=bestdist);	// current lb
			}
			
			// not useful at all
			if (!willErase) { // correct if maxedgedist used ??
				float curXdist=rgt.partial_dist;
				for (int s=0;s<gQrySize;s++) 
					if (DistMaps[s].count(NodeID)==0) {
						float tmpdist=getDist(xVal,yVal,M_xcrd[s],M_ycrd[s]);
						tmpdist=max(tmpdist,event.dist);
						//tmpdist=event.dist;	// i.e., Euclidean bounds not used
						
						if (isWeightUse) tmpdist=tmpdist*weights[s];
						if (isSumDist) 
							curXdist+=tmpdist;
						else 
							curXdist=max(tmpdist,curXdist);
					}
				if (curXdist>=bestdist) willErase=true;
			}
			if (willErase) {
				// print info.
				/*for (int s=0;s<gQrySize;s++) 
					if (DistMaps[s].count(NodeID)==0)
						if (getDist(xVal,yVal,M_xcrd[s],M_ycrd[s])>event.dist)
							printf("%d %f %f\n",NodeID,event.dist,
								getDist(xVal,yVal,M_xcrd[s],M_ycrd[s]));*/
				
				epsNbrList.erase(NodeID);
			}
		}
		
		if (bestdist<MAX_DIST) {	// at least one found
			if (epsNbrList.size()==0) break;	// all pruned => exit
		}
		
		for (int z=0;z<AdjListSize;z++) {
			getVarE(ADJNODE_A	,Ref(NewNodeID)	,AdjGrpAddr,z);
			getVarE(DIST_A		,Ref(eKdist)	,AdjGrpAddr,z);
			
			StepEvent newevent=event;	// copy ...
			newevent.node=NewNodeID;
			newevent.dist=event.dist+eKdist;
			if (DistMaps[id].count(NewNodeID)==0) sQ.push(newevent);	// propagation
			
			// if error, add cond. to (... && not on same edge) [sp]
			bool needUpdate=true;
			for (int sp=0;sp<gQrySize;sp++) 
				if (DistMaps[sp].count(NodeID)==0&&DistMaps[sp].count(NewNodeID)==0)
					needUpdate=false;
			if (needUpdate&&isNNquery) {
				float tmpposition;
				float bestVal=getOptPos(NodeID,NewNodeID,eKdist,tmpposition);
				if (bestVal>=bestdist) needUpdate=false;
			}
			if (!needUpdate) continue;
			
			getVarE(PTKEY_A	,Ref(PtGrpKey)	,AdjGrpAddr,z);
			if (PtGrpKey>=0||!isNNquery) {	// valid PtGrpKey  (-1 for invalid key)
				bool hasChanged=UpdatePoints(NodeID,NewNodeID,PtGrpKey,eKdist);
			}
		}
	}
	printf("Heap: max_size %d, pop_size %d\n",MaxHeapSize,AggPopSize);
	printf("bestdist: %f\n",bestdist);
	if (!sQ.empty()) printf("topdist: %f\n",sQ.top().dist);
	int totalmapsize=0;
	for (int sp=0;sp<gQrySize;sp++) {
		//printf("%d) td_sz: %d\n",sp,DistMaps[sp].size());
		totalmapsize+=DistMaps[sp].size();
	}
	printf("totalmapsize: %d, maxNbrSz: %d\n",totalmapsize,maxNbrSize);
}

float A_star(StepQueue& aQ,BitStore& isVisited,DISTMAP& curdistmap,
			 int dest,int& PopSize,int& VisitedSize) {
	// gdist: from src to cur ; hdist: from cur to dist
	float x_dest,y_dest,x_cur,y_cur;
	int TmpAdjGrpAddr=getAdjListGrpAddr(dest);
	getFixedF(XCRD_A,Ref(x_dest),TmpAdjGrpAddr);
	getFixedF(YCRD_A,Ref(y_dest),TmpAdjGrpAddr);
	
	while (!aQ.empty()) {
		StepEvent event=aQ.top();
		aQ.pop();	PopSize++;
		
		int NodeID=event.node;
		if (isVisited[NodeID]) continue;
		isVisited[NodeID]=true;		VisitedSize++;
		curdistmap[NodeID]=event.gdist;	// only gdist means the real dist. from src. !!!
		
		float eKdist;
		int AdjListSize,NewNodeID,AdjGrpAddr;
		AdjGrpAddr=getAdjListGrpAddr(NodeID);
		getFixedF(SIZE_A,Ref(AdjListSize),AdjGrpAddr);	// read # entries
		for (int z=0;z<AdjListSize;z++) {
			getVarE(ADJNODE_A	,Ref(NewNodeID)	,AdjGrpAddr,z);
			getVarE(DIST_A		,Ref(eKdist)	,AdjGrpAddr,z);
			
			StepEvent newevent=event;	// copy ...
			newevent.node=NewNodeID;
			newevent.gdist+=eKdist;
			newevent.hdist=0;
			
			TmpAdjGrpAddr=getAdjListGrpAddr(NewNodeID);
			getFixedF(XCRD_A,Ref(x_cur)	,TmpAdjGrpAddr);
			getFixedF(YCRD_A,Ref(y_cur)	,TmpAdjGrpAddr);
			newevent.xc=x_cur;		newevent.yc=y_cur;
						
			if (isEuclidSpace) newevent.hdist=getDist(x_cur,y_cur,x_dest,y_dest);
			newevent.dist=newevent.gdist+newevent.hdist;
			
			// pathmax equation for non-monotonic heuristic ?
			if (newevent.dist<event.dist) newevent.dist=event.dist;
			if (isVisited[NewNodeID]==false) aQ.push(newevent);	// propagation		
		}
		if (NodeID==dest) return event.dist;
	}
	return MAX_DIST;
}


struct ValuePair {
	float value;
	int id;
};

FastArray<ValuePair> covernodeset;

void Repair_astar(StepQueue& aQ,BitStore& isVisited,int dest) {
	int tsize=0;
	StepEvent* tmpAry=new StepEvent[aQ.size()];
	
	float x_dest,y_dest,x_cur,y_cur;
	int TmpAdjGrpAddr=getAdjListGrpAddr(dest);
	getFixedF(XCRD_A,Ref(x_dest),TmpAdjGrpAddr);
	getFixedF(YCRD_A,Ref(y_dest),TmpAdjGrpAddr);
	
	while (!aQ.empty()) {
		StepEvent event=aQ.top();
		aQ.pop();
		if (isVisited[event.node]) continue;
		
		// assume xc yc fields initialized by the caller !
		event.hdist=0;	// for Dijkstra
		if (isEuclidSpace) event.hdist=getDist(x_dest,y_dest,event.xc,event.yc);		
		event.dist=event.gdist+event.hdist;
		tmpAry[tsize]=event;
		tsize++;
	}
	for (int i=0;i<tsize;i++) aQ.push(tmpAry[i]);
	delete[] tmpAry;
}

void IER(StepQueue* mtQ) {
	int rt_totna = 0; // total na
	int PopSize=0,VisitedSize=0;
	float eKdist;
	int NodeID,NewNodeID,PtGrpKey,PtGrpAddr,dummy,AdjGrpAddr,AdjListSize;
	DataStore result,queries;
	
	RTree *rt = new RTree(ptgrprtfn,128);
    rt->file->isReadOnly=true;	// don't really write blocks to disk
    rt->load_root();
    
    // for rect RTree
    for (int s=0;s<gQrySize;s++) {
		DATA* newrect=new DATA;
		newrect->data[0]=newrect->data[1]=M_xcrd[s];		
		newrect->data[2]=newrect->data[3]=M_ycrd[s];
		queries.push_back(newrect);
	}
	DATA AggNNpt;
    BranchQueue AggNN_queue;
	{BranchEvent event;
	event.pnode=NULL;
	event.son=-1;
	event.mindist=0;	// always search from root node
	AggNN_queue.push(event);
	}
	
//	int nna_counter=0,zna_cnt=0;
//	BitStore* test_Visits=new BitStore[gQrySize];
//	for (int i=0;i<gQrySize;i++)
//		test_Visits[i].assign(NodeNum,false);	
	
	int test_cnt=0;
	while (true) {
		float mindist=getNextAggNN(rt,AggNN_queue,rt_totna,queries,&AggNNpt);
		if (mindist>=bestdist) break;	// using mbr of edge => correct
		test_cnt++;
		//printf("%d %d: %f %f %d\n",AggNNpt.id,AggNNpt.id2,mindist,bestdist,PopSize);
		
		NodeID=AggNNpt.id;
		AdjGrpAddr=getAdjListGrpAddr(NodeID);
		getFixedF(SIZE_A,Ref(AdjListSize),AdjGrpAddr);	// read # entries
		for (int z=0;z<AdjListSize;z++) {
			getVarE(ADJNODE_A	,Ref(NewNodeID)	,AdjGrpAddr,z);
			getVarE(DIST_A		,Ref(eKdist)	,AdjGrpAddr,z);
			getVarE(PTKEY_A		,Ref(PtGrpKey)	,AdjGrpAddr,z);
			if (NewNodeID==AggNNpt.id2) break;
		}
		if (PtGrpKey<0&&isNNquery) continue; // check if PtGrpKey valid  (-1 for invalid key)
		
		//nna_counter++;
		for (int s=0;s<gQrySize;s++) {
			DISTMAP& curDistMap=DistMaps[s];
			
			//if (test_Visits[s][NodeID]==false) {test_Visits[s][NodeID]=true;zna_cnt++;}
			if (curDistMap.count(NodeID)==0) {
				Repair_astar(mtQ[s],astar_Visits[s],NodeID);
				A_star(mtQ[s],astar_Visits[s],curDistMap,NodeID,PopSize,VisitedSize);
			}
			
			//if (test_Visits[s][NodeID]==false) {test_Visits[s][NodeID]=true;zna_cnt++;}
			if (curDistMap.count(NewNodeID)==0) {
				Repair_astar(mtQ[s],astar_Visits[s],NewNodeID);
				A_star(mtQ[s],astar_Visits[s],curDistMap,NewNodeID,PopSize,VisitedSize);
			}
		}
		bool needUpdate=true;
		if (isNNquery) {
			float tmpposition;
			float bestVal=getOptPos(NodeID,NewNodeID,eKdist,tmpposition);
			if (bestVal>=bestdist) needUpdate=false;
		}
		if (needUpdate) UpdatePoints(NodeID,NewNodeID,PtGrpKey,eKdist);
	}
	//printf("bestdist: %f, net-na: %d,%d\n",bestdist,nna_counter,zna_cnt);
	printf("bestdist: %f\n",bestdist);
	printf("pop size: %d, VisitedSize: %d, test_cnt: %d\n",PopSize,VisitedSize,test_cnt);
	printf("page faults: %d, totna: %d\n",rt->file->page_faults,rt_totna);
}

inline void printSetting() {
	if (isNNquery) printf(",point"); else printf(",center");
	if (isSumDist) printf(",sum"); else printf(",max");
}

enum ALGTYPE {Alg_CE,Alg_IER};

// Only find the necessary dist, not necessarily all pairs dist
void FindSoln(ALGTYPE curAlg) {	// "node" reused as the next nodeID instead !!!
	StepQueue sQ;
	StepEvent stepL,stepR;
	StepQueue mtQ[gQrySize];
	
	// initialization
	bestdist=MAX_DIST;
	InitDQ(dQ_graph);	// clear kNN-dist heap
	
	onSameEdge.assign(gQrySize,false);
	partdists=new float*[gQrySize];
	for (int i=0;i<gQrySize;i++) partdists[i]=new float[2];
	compdist.assign(gQrySize,MAX_DIST);
	bestcompdist.assign(gQrySize,MAX_DIST);
	
	float vP,eDist,eKdist,x1,y1,x2,y2,x,y;
	int AdjListSize,Ni,Nj,NewNodeID,TmpAdjGrpAddr;
	DATA* centroid=new DATA;
	
	centroid->data[0]=centroid->data[1]=0;
	while (!sQ.empty()) sQ.pop();
	for (int s=0;s<gQrySize;s++) {
		vP=gQryPts[s].vP;
		Ni=gQryPts[s].Ni;
		Nj=gQryPts[s].Nj;
		
		// for greedy search and rapid elimination
		TmpAdjGrpAddr=getAdjListGrpAddr(Ni);
		getFixedF(XCRD_A,Ref(x1),TmpAdjGrpAddr);
		getFixedF(YCRD_A,Ref(y1),TmpAdjGrpAddr);
		
		// lookup eDist from adj. list. of Ni;
		getFixedF(SIZE_A,Ref(AdjListSize),TmpAdjGrpAddr);	// read # entries
		for (int z=0;z<AdjListSize;z++) {
			getVarE(ADJNODE_A	,Ref(NewNodeID)	,TmpAdjGrpAddr,z);
			getVarE(DIST_A		,Ref(eKdist)	,TmpAdjGrpAddr,z);
			if (NewNodeID==Nj) eDist=eKdist;
		}
		
		TmpAdjGrpAddr=getAdjListGrpAddr(Nj);
		getFixedF(XCRD_A,Ref(x2),TmpAdjGrpAddr);
		getFixedF(YCRD_A,Ref(y2),TmpAdjGrpAddr);
		x=x1+(vP/eDist)*(x2-x1);	y=y1+(vP/eDist)*(y2-y1);
		M_xcrd.push_back(x);		M_ycrd.push_back(y);
		
		stepL.ClusID=stepR.ClusID=s;
		stepL.dist=vP;			stepL.node=Ni;
		stepR.dist=eDist-vP;	stepR.node=Nj;
		
		stepL.isSortAcc=true;	stepR.isSortAcc=true;	// ??? 4/3/2004 19:00
		sQ.push(stepL);			sQ.push(stepR);
		
		stepL.xc=x1;	stepL.yc=y1;
		stepR.xc=x2;	stepR.yc=y2;
		
		if (isWeightUse) {
			centroid->data[0]+=x*weights[s];	centroid->data[1]+=y*weights[s];
		} else {
			centroid->data[0]+=x;	centroid->data[1]+=y;
		}
		stepL.gdist=vP;			stepL.hdist=0;
		stepR.gdist=eDist-vP;	stepR.hdist=0;
		
		stepL.isSortAcc=false;	stepR.isSortAcc=false;	// ??? 4/3/2004 19:00
		mtQ[s].push(stepL);		mtQ[s].push(stepR);
	}
	centroid->data[0]/=gQrySize;	centroid->data[1]/=gQrySize;
	
	for (int sp=0;sp<gQrySize;sp++) 
		DistMaps[sp].clear(); // clear first (safe before use)
	
	if (curAlg==Alg_CE) {
		printf("\n[CE");
		printSetting();
		printf("]\n");
		ConcurrentExpansion(sQ);
		//RR_Expand(mtQ);
	}
		
	// later, change the ... for cache
	if (curAlg==Alg_IER) {
		printf("\n[IER");
		printSetting();
		printf("]\n");
		astar_Visits=new BitStore[gQrySize];
		for (int i=0;i<gQrySize;i++)
			astar_Visits[i].assign(NodeNum,false);
		IER(mtQ);
	}
}

// error for 10000 times: 2.804688
int main(int argc,char** argv) {
	ConfigType cr;
	AddConfigFromFile(cr,"config.prop");
	AddConfigFromCmdLine(cr,argc,argv);
	ListConfig(cr);
	
	string filename=getConfigStr("DBASE",cr);
	filename=filename+getConfigStr("data",cr);
	NNnum=getConfigInt("NNnum",cr);
	const char* fileprefix=filename.c_str();
  	const char* edgert=getConfigStr("edgert",cr,false,"");
  	if (strcmp(edgert,"")==0)
  		sprintf(ptgrprtfn,"%s.p_rt",fileprefix);
  	else 
  		sprintf(ptgrprtfn,"data/%s.ed_rt",edgert);
  	
  	char queryf[255];
  	sprintf(queryf,"query/%s.qry",getConfigStr("query",cr));
  	ReadQueryFile(queryf);
  	
	InitClock();	// side effect: initialize the seeds for 2 random generators
	//OpenDiskComm(fileprefix,100);	// 100,200,400,800
	OpenDiskComm(fileprefix,DEFAULT_CACHESIZE);
	InitClock();	// don't count file reading time
	
	// if seg. fault, check binary file type
	//isEuclidSpace=false;
	//isNNquery=false;	// center query
	//isSumDist=false;	// d_{max}
	if (!isNNquery) NNnum=1;	// force to find 1 center only !
	DistMaps=new DISTMAP[gQrySize];
	
	//isWeightUse=true;	// default: false
	if (isWeightUse) {
		weights=new float[gQrySize];
		float theta=getConfigFloat("qryskew",cr);
		float avgSkew=0,maxSkew=0;
	  	for(int s=1;s<=gQrySize;s++)  {
	    	weights[s-1]=1.0/pow( (double)s, (double)theta);
	    	avgSkew+=weights[s-1];
	    	maxSkew=max(maxSkew,weights[s-1]);
	    }
		avgSkew=avgSkew/gQrySize;
	    // normalize: sum => sum(wgts)=|Q| ; max => max(wgts)=1
	    printf("wgts=(");
	    for (int s=0;s<gQrySize;s++) {
	    	weights[s]/=(isSumDist)?(avgSkew):(maxSkew);
	    	printf("%0.3f ",weights[s]);
	    }
	    printf(")\n");
    }
	srand(0); srand48(0);	// reduce random effect
	
//	for (int t=0;t<10;t++) { 
		if (isWeightUse) {
			for (int x=0;x<100;x++) {  // shuffle weights
				int y=rand()%gQrySize,z=rand()%gQrySize;
				float tmpw=weights[y];
				weights[y]=weights[z];
				weights[z]=tmpw;
			}
		}
	
		RefreshCache();
		FindSoln(Alg_CE);
		//FindSoln(Alg_IER);		
//	}
	
	CloseDiskComm();
	PrintElapsed();
	
	if (isWeightUse) {
		for (int s=0;s<gQrySize;s++) {
			printf("%d: %f\t",s,bestcompdist[s]);
			if ((s+1)%4==0) printf("\n");
		}
	}
	printf("\n");
	return 0;
}

/*
// memory requirement: r bitmap and r tempdiststore
void ConcurrentExpansion(StepQueue& sQ) {
	int MaxHeapSize=0,AggPopSize=0;
	
	// initialize variables
	if (!isEuclidSpace) return;
	
	// find crd. of r query points
	REGIONMAP epsNbrList;
	BitStore epsBits;
	epsBits.assign(NodeNum,false);
	int maxNbrSize=0;
	float global_topdist=0;
	while (!sQ.empty()) {
		if (MaxHeapSize<sQ.size()) MaxHeapSize=sQ.size();
		StepEvent event=sQ.top();
		sQ.pop();	AggPopSize++;
		
		int NodeID=event.node;
		int id=event.ClusID;
		
		if (DistMaps[id].count(NodeID)>0) continue;	// already found
		DistMaps[id][NodeID]=event.dist;
		
		int AdjListSize,NewNodeID,AdjGrpAddr,PtGrpKey;
		float eKdist,xVal,yVal;
		AdjGrpAddr=getAdjListGrpAddr(NodeID);
		getFixedF(SIZE_A,Ref(AdjListSize),AdjGrpAddr);	// read # entries
		getFixedF(XCRD_A,Ref(xVal)	,AdjGrpAddr);
		getFixedF(YCRD_A,Ref(yVal)	,AdjGrpAddr);
		
		global_topdist=max(global_topdist,event.dist);
		maxNbrSize=max(maxNbrSize,epsNbrList.size());
		bool needPushNeighbor=true;
		
		if (epsBits[NodeID]==false) {
			epsBits[NodeID]=true;
			epsNbrList[NodeID].crash=0;
			epsNbrList[NodeID].partial_dist=0;
		}
		RgnType& rgt=epsNbrList[NodeID];
		float tmp_lbagdist=0;
		rgt.crash+=1;	// update crash
		if (isSumDist) {	// update cur dist and check if need prune
			epsNbrList[NodeID].partial_dist+=event.dist;
			if ((gQrySize-rgt.crash)*event.dist+rgt.partial_dist>=bestdist) 
				needPushNeighbor=false;
			tmp_lbagdist=(gQrySize-rgt.crash)*event.dist+rgt.partial_dist;
		} else { 	// update cur dist and check if need prune
			epsNbrList[NodeID].partial_dist=max(event.dist,epsNbrList[NodeID].partial_dist);
			if (rgt.partial_dist>=bestdist)
				needPushNeighbor=false;
			tmp_lbagdist=rgt.partial_dist;
		}
		
		for (int z=0;z<AdjListSize;z++) {
			getVarE(ADJNODE_A	,Ref(NewNodeID)	,AdjGrpAddr,z);
			getVarE(DIST_A		,Ref(eKdist)	,AdjGrpAddr,z);
			
			StepEvent newevent=event;	// copy ...
			newevent.node=NewNodeID;
			newevent.dist=event.dist+eKdist;
			
			if (DistMaps[id].count(NewNodeID)==0) {// unvisited
				if (needPushNeighbor) sQ.push(newevent);
				else printf("%d, %f\n",id,event.dist);
//				if (isSumDist) {
//					if (tmp_lbagdist<bestdist+gQrySize*eKdist) sQ.push(newevent);
//				} else {
//					if (tmp_lbagdist<bestdist+eKdist) sQ.push(newevent);
//				}
			}
			
			// if error, add cond. to (... && not on same edge) [sp]
			bool needUpdate=true;
			for (int sp=0;sp<gQrySize;sp++) 
				if (DistMaps[sp].count(NodeID)==0&&DistMaps[sp].count(NewNodeID)==0)
					needUpdate=false;
			if (needUpdate&&isNNquery) {
				float bestVal=getOptPos(NodeID,NewNodeID,eKdist);
				if (bestVal>=bestdist) needUpdate=false;
			}
			if (!needUpdate) continue;
			getVarE(PTKEY_A	,Ref(PtGrpKey)	,AdjGrpAddr,z);
			if (PtGrpKey>=0||!isNNquery) 	// valid PtGrpKey  (-1 for invalid key)
				bool hasChanged=UpdatePoints(NodeID,NewNodeID,PtGrpKey,eKdist);
		}
		
	}
	printf("Heap: max_size %d, pop_size %d\n",MaxHeapSize,AggPopSize);
	printf("bestdist: %f\n",bestdist);
	if (!sQ.empty()) printf("topdist: %f\n",sQ.top().dist);
	printf("gtd: %f\n",global_topdist);
	int totalmapsize=0;
	for (int sp=0;sp<gQrySize;sp++) {
		//printf("%d) td_sz: %d\n",sp,DistMaps[sp].size());
		totalmapsize+=DistMaps[sp].size();
	}
	printf("totalmapsize: %d, maxNbrSz: %d\n",totalmapsize,maxNbrSize);
}

*/


