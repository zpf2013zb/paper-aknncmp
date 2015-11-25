#include "utility.h"
#include "netshare.h"
#include "diskbased.h"

typedef map<int,float> DISTMAP;
DISTMAP* DistMaps;	// only use cheap storage for range search
BitStore *raVisited,*saVisited;		// tmp. bitmap for A* search

typedef vector<float> 	FloatVec;

FloatVec M_xcrd,M_ycrd,compdist,bestcompdist;

BitStore onSameEdge;
float** partdists;
float bestdist;	// best value found so far

bool isNNquery=true;
bool isEuclidSpace=true;

DistQueue dQ_graph;
int PrintLimit;

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
			// bestVal-=eKdist;	// this line is different from that in gnn.cc !!!
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
	/*if (hasChanged) 
		printf("c: %d %d %f\n",NodeID,NewNodeID,eKdist);*/
	return hasChanged;
}

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
		
		if (isWeightUse) mincpdist=mincpdist*weights[sp];
		if (isSumDist)
			dist+=mincpdist;
		else 
			dist=max(mincpdist,dist);
	}
	return dist;
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
		// currently just degen. to Dijkstra
		event.hdist=0;
		event.hdist=getDist(x_dest,y_dest,event.xc,event.yc);		
		event.dist=event.gdist+event.hdist;
		tmpAry[tsize]=event;
		tsize++;
	}
	for (int i=0;i<tsize;i++) aQ.push(tmpAry[i]);
	delete[] tmpAry;
}

void getNextNode(float& nextdist,int& NextNodeID,StepQueue& saQueue,BitStore& isVisited,float& xc,float& yc) {
	float eKdist,xVal,yVal;
	int AdjListSize,NewNodeID,AdjGrpAddr;
	while (!saQueue.empty()) {
		StepEvent event=saQueue.top();
		saQueue.pop();
		int NodeID=event.node;
		if (isVisited[NodeID]) continue;
		isVisited[NodeID]=true;	// !!!
		
		AdjGrpAddr=getAdjListGrpAddr(NodeID);
		getFixedF(XCRD_A,Ref(xVal)	,AdjGrpAddr);
		getFixedF(YCRD_A,Ref(yVal)	,AdjGrpAddr);
		getFixedF(SIZE_A,Ref(AdjListSize),AdjGrpAddr);	// read # entries
		for (int z=0;z<AdjListSize;z++) {
			getVarE(ADJNODE_A	,Ref(NewNodeID)	,AdjGrpAddr,z);
			getVarE(DIST_A		,Ref(eKdist)	,AdjGrpAddr,z);
			StepEvent newevent=event;	// copy ...
			newevent.node=NewNodeID;
			newevent.dist=event.dist+eKdist;
			
			//printf("%f %f %f\n",event.dist,eKdist,newevent.dist);
			
			if (isVisited[NewNodeID]==false) saQueue.push(newevent);	// propagation		
		}
		
		xc=xVal;	yc=yVal;
		NextNodeID=event.node;
		nextdist=event.dist;
		return;
	}
	NextNodeID=0;
	nextdist=MAX_DIST;	// meaning empty queue
}

struct RgnType {
	float partial_dist;
	int crash;
	BitStore bset;
};

typedef map<int,RgnType> REGIONMAP;

/*void NRA(StepQueue* saQ) {
	int PopSize=0,VisitedSize=0;
	float eKdist;
	int AdjListSize,NewNodeID,AdjGrpAddr,PtGrpKey;
	
	// dijkstra's expansion in round robin
	int curround,totround=0;	// # of sorted access
	int NodeID;
	float nextdist;
	float cur_x,cur_y;
	
	float lb_dist[gQrySize];
	fill(lb_dist,lb_dist+gQrySize,0);
	
	BitStore epsBits;
	epsBits.assign(NodeNum,false);
	REGIONMAP epsNbrList;
	
	while (true) {
		// get next of curround 
		curround=totround%gQrySize;
		getNextNode(nextdist,NodeID,saQ[curround],saVisited[curround],cur_x,cur_y);
		if (nextdist==MAX_DIST) {
			totround++;	continue;
		}
		
		DistMaps[curround][NodeID]=nextdist;	// assume valid result
		lb_dist[curround]=nextdist;
		
		if (bestdist<MAX_DIST) {
			if (epsNbrList.size()==0) break;
		}
		
		float mindist=0;
		for (int s=0;s<gQrySize;s++) {
			if (isSumDist)
				mindist+=lb_dist[s];
			else
				mindist=max(mindist,lb_dist[s]);
			//if (totround%PrintLimit==PrintLimit-1) printf("%f ",lb_dist[s]);
		}
		if (totround%PrintLimit==PrintLimit-1) {
			//printf("\n");
			printf("%d: %f %f %d\n",totround,mindist,bestdist,epsNbrList.size());
			PrintElapsed();
		}
		//if (mindist>=bestdist) break;
		
		if (epsBits[NodeID]==false) {
			epsBits[NodeID]=true;
			epsNbrList[NodeID].crash=1;
			epsNbrList[NodeID].partial_dist=nextdist;
			epsNbrList[NodeID].bset.assign(gQrySize,false);
			epsNbrList[NodeID].bset[curround]=true;
		} else if (epsNbrList.count(NodeID)>0) {
			// compute lb and prune ?
			epsNbrList[NodeID].crash+=1;	// update crash
			RgnType& rgt=epsNbrList[NodeID];
			epsNbrList[NodeID].bset[curround]=true;
			
			if (isSumDist) {	// update cur dist and check if need prune
				rgt.partial_dist+=nextdist;
			} else {
				rgt.partial_dist=max(nextdist,rgt.partial_dist);
			}
		}
		
		// pruning only after soln. found !
		if (totround%PrintLimit==PrintLimit-1&&bestdist<MAX_DIST) {
			//printf("b: %d\n",epsNbrList.size());
			REGIONMAP::iterator curIter=epsNbrList.begin();
			while (curIter!=epsNbrList.end()) {
				RgnType rgt=curIter->second;
				bool willErase=(rgt.crash>=gQrySize);
				float agdist_s=rgt.partial_dist;
				for (int s=0;s<gQrySize;s++)
					if (rgt.bset[s]==false) {
						float tmpdist=getDist(cur_x,cur_y,M_xcrd[s],M_xcrd[s]);
						tmpdist=0;
						if (lb_dist[s]>tmpdist) tmpdist=lb_dist[s];
						if (isSumDist) 
							agdist_s+=tmpdist;
						else
							agdist_s=max(agdist_s,tmpdist);
					}
				// check if need prune
				if (!willErase) willErase=(agdist_s>=bestdist);	// current lb
				if (willErase) epsNbrList.erase(curIter->first);
				curIter++;
			}
			//printf("a: %d\n",epsNbrList.size());
		}
		
		AdjGrpAddr=getAdjListGrpAddr(NodeID);
		getFixedF(SIZE_A,Ref(AdjListSize),AdjGrpAddr);	// read # entries
		for (int z=0;z<AdjListSize;z++) {
			getVarE(ADJNODE_A	,Ref(NewNodeID)	,AdjGrpAddr,z);
			getVarE(DIST_A		,Ref(eKdist)	,AdjGrpAddr,z);
			getVarE(PTKEY_A		,Ref(PtGrpKey)	,AdjGrpAddr,z);
			if (PtGrpKey>=0||!isNNquery) {	// valid PtGrpKey  (-1 for invalid key)
				bool needUpdate=true;
				if (needUpdate&&isNNquery) {
					float bestVal=getOptPos(NodeID,NewNodeID,eKdist);
					if (bestVal>=bestdist) needUpdate=false;
				}
				if (needUpdate) UpdatePoints(NodeID,NewNodeID,PtGrpKey,eKdist);
			}
		}
		totround++;
	}
	printf("bestdist: %f\n",bestdist);
	printf("totround: %d, pop size: %d, VisitedSize: %d\n",totround,PopSize,VisitedSize);
}*/

// round robin sorted access (Dijkstra's alg.)
/*void TA(StepQueue* raQ,StepQueue* saQ) {
	float eKdist;
	int AdjListSize,NewNodeID,AdjGrpAddr,PtGrpKey;
	int PopSize=0,VisitedSize=0,NodeID;
	float cur_x,cur_y,nextdist;
	int curround,totround=0;	// # of sorted access
	
	float lb_dist[gQrySize];
	fill(lb_dist,lb_dist+gQrySize,0);
	while (true) {
		int sumsz=0;
		for (int s=0;s<gQrySize;s++) sumsz+=saQ[s].size();
		if (sumsz==0) break;
		
		// get next of curround
		curround=totround%gQrySize;
		while (saQ[curround].size()==0) {
			totround++;
			curround=totround%gQrySize;
		}
		getNextNode(nextdist,NodeID,saQ[curround],saVisited[curround],cur_x,cur_y);
		if (nextdist==MAX_DIST) {
			totround++;	continue;
		}
		
		DistMaps[curround][NodeID]=nextdist;	// assume valid result
		lb_dist[curround]=nextdist;
		
		float mindist=0,curnetdist=0;
		for (int s=0;s<gQrySize;s++) {
			float cpdist=getDist(cur_x,cur_y,M_xcrd[s],M_xcrd[s]);
			if (DistMaps[s].count(NodeID)>0) cpdist=DistMaps[s][NodeID];
			if (isSumDist) {
				mindist+=lb_dist[s];
				curnetdist+=cpdist;
			} else {
				mindist=max(mindist,lb_dist[s]);
				curnetdist=max(curnetdist,cpdist);
			}
		}
		if (totround%PrintLimit==PrintLimit-1) {
			printf("%d: %f %f\n",totround,mindist,bestdist);
			PrintElapsed();
		}
		totround++;
		if (mindist>=bestdist) break;
		
		// incorrect bound, just use to speed up !!! 
		if (isSumDist) {
			if (curnetdist>=bestdist+gQrySize*eKdist) continue;		// cannot have better distance
		} else {
			if (curnetdist>=bestdist+gQrySize) continue;			// cannot have better distance
		}
		
		
//		bool mustPrune=true;
//		for (int z=0;z<AdjListSize;z++) {
//			getVarE(ADJNODE_A	,Ref(NewNodeID)	,AdjGrpAddr,z);
//			getVarE(DIST_A		,Ref(eKdist)	,AdjGrpAddr,z);
//			getVarE(PTKEY_A		,Ref(PtGrpKey)	,AdjGrpAddr,z);
//			// L0^- filter
//			if (isSumDist) {
//				if (curnetdist<bestdist+gQrySize*eKdist) mustPrune=false;
//			} else {
//				if (curnetdist<bestdist+eKdist) mustPrune=false;
//			}
//		}
//		if (mustPrune) continue;
		
		AdjGrpAddr=getAdjListGrpAddr(NodeID);
		getFixedF(SIZE_A,Ref(AdjListSize),AdjGrpAddr);	// read # entries
		for (int s=0;s<gQrySize;s++) {
			DISTMAP& curDistMap=DistMaps[s];
			if (curDistMap.count(NodeID)==0) {
				Repair_astar(raQ[s],raVisited[s],NodeID);
				A_star(raQ[s],raVisited[s],curDistMap,NodeID,PopSize,VisitedSize);
			}
		}
		for (int z=0;z<AdjListSize;z++) {
			getVarE(ADJNODE_A	,Ref(NewNodeID)	,AdjGrpAddr,z);
			getVarE(DIST_A		,Ref(eKdist)	,AdjGrpAddr,z);
			getVarE(PTKEY_A		,Ref(PtGrpKey)	,AdjGrpAddr,z);
			
			// L0 filter
			if (isSumDist) {
				if (curnetdist>=bestdist+gQrySize*eKdist) continue;		// cannot have better distance
			} else {
				if (curnetdist>=bestdist+eKdist) continue;		// cannot have better distance
			}
			// L1 filter
			if (getMinOptPos(NodeID,NewNodeID,eKdist)>=bestdist) continue;
			if (PtGrpKey>=0||!isNNquery) {	// valid PtGrpKey  (-1 for invalid key)
				for (int s=0;s<gQrySize;s++) {
					DISTMAP& curDistMap=DistMaps[s];
					if (curDistMap.count(NewNodeID)==0) {
						Repair_astar(raQ[s],raVisited[s],NewNodeID);
						A_star(raQ[s],raVisited[s],curDistMap,NewNodeID,PopSize,VisitedSize);
					}
				}
				bool needUpdate=true;
				if (needUpdate&&isNNquery) {
					float bestVal=getOptPos(NodeID,NewNodeID,eKdist);
					if (bestVal>=bestdist) needUpdate=false;
				}
				if (needUpdate) UpdatePoints(NodeID,NewNodeID,PtGrpKey,eKdist);
			}
		}
	}
	printf("bestdist: %f\n",bestdist);
	printf("totround: %d, pop size: %d, VisitedSize: %d\n",totround,PopSize,VisitedSize);
}*/

// round robin sorted access (Dijkstra's alg.)
void TA_EW(StepQueue* raQ,StepQueue& sQ) {
	float eKdist;
	int AdjListSize,NewNodeID,AdjGrpAddr,PtGrpKey;
	int PopSize=0,VisitedSize=0,NodeID;
	float cur_x,cur_y,nextdist;
	int curround,totround=0;	// # of sorted access
	
	float lb_dist=0;
	while (!sQ.empty()) {
		StepEvent event=sQ.top();
		sQ.pop();
		NodeID=event.node;
		curround=event.ClusID;
		if (saVisited[curround][NodeID]) continue;
		saVisited[curround][NodeID]=true;	// !!!
		
		AdjGrpAddr=getAdjListGrpAddr(NodeID);
		getFixedF(SIZE_A,Ref(AdjListSize),AdjGrpAddr);	// read # entries
		for (int z=0;z<AdjListSize;z++) {
			getVarE(ADJNODE_A	,Ref(NewNodeID)	,AdjGrpAddr,z);
			getVarE(DIST_A		,Ref(eKdist)	,AdjGrpAddr,z);
			StepEvent newevent=event;	// copy ...
			newevent.node=NewNodeID;
			newevent.dist=event.dist+eKdist;
			if (saVisited[curround][NewNodeID]==false) sQ.push(newevent);	// propagation		
		}
		nextdist=event.dist;
		
		// get next of curround 
		DistMaps[curround][NodeID]=nextdist;	// assume valid result
		lb_dist=nextdist;
		
		float mindist=0,curnetdist=0;
		for (int s=0;s<gQrySize;s++) {
			float cpdist=getDist(cur_x,cur_y,M_xcrd[s],M_xcrd[s]);
			if (DistMaps[s].count(NodeID)>0) cpdist=DistMaps[s][NodeID];
			
			if (isWeightUse) cpdist=cpdist*weights[s];
			if (isSumDist) {
				mindist+=lb_dist;
				curnetdist+=cpdist;
			} else {
				mindist=max(mindist,lb_dist);
				curnetdist=max(curnetdist,cpdist);
			}
		}
		if (totround%PrintLimit==PrintLimit-1) {
			printf("%d: %f %f\n",totround,mindist,bestdist);
			PrintElapsed();
		}
		totround++;
		if (mindist>=bestdist) break;
		
		bool mustPrune=true;
		for (int z=0;z<AdjListSize;z++) {
			getVarE(ADJNODE_A	,Ref(NewNodeID)	,AdjGrpAddr,z);
			getVarE(DIST_A		,Ref(eKdist)	,AdjGrpAddr,z);
			getVarE(PTKEY_A		,Ref(PtGrpKey)	,AdjGrpAddr,z);
			// L0^- filter
			if (isSumDist) {
				if (curnetdist<bestdist+gQrySize*eKdist) mustPrune=false;
			} else {
				if (curnetdist<bestdist+eKdist) mustPrune=false;
			}
		}
		if (mustPrune) continue;
				
		AdjGrpAddr=getAdjListGrpAddr(NodeID);
		getFixedF(SIZE_A,Ref(AdjListSize),AdjGrpAddr);	// read # entries
		for (int s=0;s<gQrySize;s++) {
			DISTMAP& curDistMap=DistMaps[s];
			if (curDistMap.count(NodeID)==0) {
				Repair_astar(raQ[s],raVisited[s],NodeID);
				A_star(raQ[s],raVisited[s],curDistMap,NodeID,PopSize,VisitedSize);
			}
		}
		for (int z=0;z<AdjListSize;z++) {
			getVarE(ADJNODE_A	,Ref(NewNodeID)	,AdjGrpAddr,z);
			getVarE(DIST_A		,Ref(eKdist)	,AdjGrpAddr,z);
			getVarE(PTKEY_A		,Ref(PtGrpKey)	,AdjGrpAddr,z);
			
			// L0 filter
			if (isSumDist) {
				if (curnetdist>=bestdist+gQrySize*eKdist) continue;		// cannot have better distance
			} else {
				if (curnetdist>=bestdist+eKdist) continue;		// cannot have better distance
			}
			
			// L1 filter
			if (getMinOptPos(NodeID,NewNodeID,eKdist)>=bestdist) continue;
			if (PtGrpKey>=0||!isNNquery) {	// valid PtGrpKey  (-1 for invalid key)
				for (int s=0;s<gQrySize;s++) {
					DISTMAP& curDistMap=DistMaps[s];
					if (curDistMap.count(NewNodeID)==0) {
						Repair_astar(raQ[s],raVisited[s],NewNodeID);
						A_star(raQ[s],raVisited[s],curDistMap,NewNodeID,PopSize,VisitedSize);
					}
				}
				bool needUpdate=true;
				if (needUpdate&&isNNquery) {
					float tmppossition;
					float bestVal=getOptPos(NodeID,NewNodeID,eKdist,tmppossition);
					if (bestVal>=bestdist) needUpdate=false;
				}
				if (needUpdate) UpdatePoints(NodeID,NewNodeID,PtGrpKey,eKdist);
			}
		}
	}
	printf("bestdist: %f\n",bestdist);
	printf("totround: %d, pop size: %d, VisitedSize: %d\n",totround,PopSize,VisitedSize);
}

inline void printSetting() {
	if (isNNquery) printf(",point"); else printf(",center");
	if (isSumDist) printf(",sum"); else printf(",max");
}

enum ALGTYPE {Alg_TA,Alg_TA_EW,Alg_NRA};

// Only find the necessary dist, not necessarily all pairs dist
void FindSoln(ALGTYPE curAlg) {	// "node" reused as the next nodeID instead !!!
	StepQueue sQ;
	StepEvent stepL,stepR;
	StepQueue raQ[gQrySize],saQ[gQrySize];
	
	// initialization
	bestdist=MAX_DIST;
	InitDQ(dQ_graph);	// clear kNN-dist heap
	
	onSameEdge.assign(gQrySize,false);
	partdists=new (float*)[gQrySize];
	for (int i=0;i<gQrySize;i++) partdists[i]=new float[2];
	compdist.assign(gQrySize,MAX_DIST);
	bestcompdist.assign(gQrySize,MAX_DIST);
	
	float vP,eDist,eKdist,x1,y1,x2,y2,x,y;
	int AdjListSize,Ni,Nj,NewNodeID,TmpAdjGrpAddr;
	
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
		sQ.push(stepL);			sQ.push(stepR);
		
		stepL.xc=x1;	stepL.yc=y1;
		stepR.xc=x2;	stepR.yc=y2;
		
		stepL.gdist=vP;			stepL.hdist=0;
		stepR.gdist=eDist-vP;	stepR.hdist=0;
		
		stepL.isSortAcc=false;	stepR.isSortAcc=false;	// ??? 4/3/2004 18:55
		raQ[s].push(stepL);		raQ[s].push(stepR);
		
		stepL.isSortAcc=true;	stepR.isSortAcc=true;	// ??? 4/3/2004 18:55
		saQ[s].push(stepL);		saQ[s].push(stepR);
	}
	
	for (int sp=0;sp<gQrySize;sp++) 
		DistMaps[sp].clear(); // clear first (safe before use)
	
	saVisited=new BitStore[gQrySize];
	raVisited=new BitStore[gQrySize];
	for (int i=0;i<gQrySize;i++) {
		saVisited[i].assign(NodeNum,false);
		raVisited[i].assign(NodeNum,false);
	}
	
	if (curAlg==Alg_TA) {
		printf("\n[TA");
		printSetting();
		printf("]\n");
		//TA(raQ,saQ);
	}
	
	if (curAlg==Alg_TA_EW) {
		printf("\n[TA_EW");
		printSetting();
		printf("]\n");
		TA_EW(raQ,sQ);
	}
	
	if (curAlg==Alg_NRA) {
		printf("\n[NRA");
		printSetting();
		printf("]\n");
		//NRA(saQ);
	}
}

// note: this alg. causes thrashing (blowing the cache) when ...
int main(int argc,char** argv) {
	ConfigType cr;
	AddConfigFromFile(cr,"config.prop");
	AddConfigFromCmdLine(cr,argc,argv);
	ListConfig(cr);
	
	string filename=getConfigStr("DBASE",cr);
	filename=filename+getConfigStr("data",cr);
	NNnum=getConfigInt("NNnum",cr);
	
	const char* fileprefix=filename.c_str();
	
  	char queryf[255];
  	sprintf(queryf,"query/%s.qry",getConfigStr("query",cr));
  	ReadQueryFile(queryf);
  	  	
	InitClock();	// side effect: initialize the seeds for 2 random generators
	//OpenDiskComm(fileprefix,100);	// 100,200,400,800
	OpenDiskComm(fileprefix,DEFAULT_CACHESIZE);
	InitClock();	// don't count file reading time
	
	//isEuclidSpace=false;	// replace A* by Dijkstra !
	isSumDist=false;	// d_{max}
	if (!isNNquery) NNnum=1;	// force to find 1 center only !
	
	DistMaps=new DISTMAP[gQrySize];
	PrintLimit=100000;
	
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
	if (isWeightUse) {
		for (int x=0;x<100;x++) {  // shuffle weights
			int y=rand()%gQrySize,z=rand()%gQrySize;
			float tmpw=weights[y];
			weights[y]=weights[z];
			weights[z]=tmpw;
		}
	}
	
	RefreshCache();
	FindSoln(Alg_TA_EW);
	
	// don't consider !!!
	//FindSoln(Alg_TA);
	//FindSoln(Alg_NRA);
	
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
