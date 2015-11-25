#include "utility.h"
#include "netshare.h"
#include "diskbased.h"
#include "rtree.h"

typedef map<int,float> DISTMAP;
DISTMAP* DistMaps;	// only use cheap storage for range search

typedef vector<float> 	FloatVec;

FloatVec M_vP,M_xcrd,M_ycrd,compdist,bestcompdist;

int maxCrdSq;
FastArray<int> Partition;
FastArray<int> bnodecnt;
FastArray<int>* bnodeset;
BitStore isBoundary;

BitStore onSameEdge;
float** partdists;
float bestdist;	// best value found so far

float* weights;
bool isWeightUse=false;	// only for d_{max}
bool isNNquery=true;
bool isSumDist=true;
char nodertrfn[255];

#include "rtacc.h"
DistQueue dQ_graph;


bool isEuclidSpace=true;
bool isPrintVis=false;
StepEvent** initSteps;

inline void printNodeVis(int node) {
	if (isPrintVis) {
		float x_cur,y_cur;
		int TmpAdjGrpAddr=getAdjListGrpAddr(node);
		getFixedF(XCRD_A,Ref(x_cur)	,TmpAdjGrpAddr);
		getFixedF(YCRD_A,Ref(y_cur)	,TmpAdjGrpAddr);
		printf("%f %f\n",x_cur,y_cur);
	}
}

inline void printQueueVis(StepQueue& aQ,BitStore& isVisited) {
	if (isPrintVis) {
		printf("----------------\n");
		while (!aQ.empty()) {
			StepEvent event=aQ.top();
			aQ.pop();
			int NodeID=event.node;
			if (isVisited[NodeID]) continue;
			isVisited[NodeID]=true;
			printNodeVis(NodeID);
		}
	}
}

float getDist(float x1,float y1,float x2,float y2) {
	return float(sqrt((x2-x1)*(x2-x1)+(y2-y1)*(y2-y1)));
}

void ReadHiTiDist(const char* matf) {
	int dummy_NodeNum;
	float value;
	
	FILE* cmat=fopen(matf,"r");
	CheckFile(cmat,matf);
	
	fread(&dummy_NodeNum,1,sizeof(int),cmat);
	if (dummy_NodeNum!=NodeNum) exit(0);
	
	fread(&maxCrdSq,1,sizeof(int),cmat);
	Partition.assign(NodeNum,-1);
	fread(&(Partition[0]),NodeNum,sizeof(int),cmat);
	
	// adjlist, edgemap for L1 graph now !    31/1/2004
	bnodecnt.assign(maxCrdSq,0);
	bnodeset=new FastArray<int>[maxCrdSq];
	isBoundary.assign(NodeNum,false);
	AdjList=new FastArray<int>[NodeNum];
	EdgeMap.clear();
	for (int part=0;part<maxCrdSq;part++) {
		int bn;
		fread(&bn,1,sizeof(int),cmat);		// # of border nodes
		bnodeset[part].assign(bn,-1);
		fread(&(bnodeset[part][0]),bn,sizeof(int),cmat);	// border node IDs
		//printf("[%d: %d]\n",part,bn);
		
		for (int j=0;j<bn;j++)
			isBoundary[bnodeset[part][j]]=true;
		
		for (int cur_pos=1;cur_pos<bn;cur_pos++) {
			// write within links
			for (int t=0;t<cur_pos;t++) {
				int Ni=bnodeset[part][t];
				int Nj=bnodeset[part][cur_pos];
				
				fread(&value,1,sizeof(float),cmat);
				edge* e=new edge;
				e->dist=value;
				e->Ni=Ni;	e->Nj=Nj;	// enforce the constriant Ni<Nj
				AdjList[Ni].push_back(Nj);	AdjList[Nj].push_back(Ni);
			   	EdgeMap[getKey(Ni,Nj)]=e;
				//printf("%f ",value);
			}
			//printf("\n");
		}
	}
	printf("|E1| %d, filesize: %d\n",EdgeMap.size(),ftell(cmat));
	fclose(cmat);
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

inline void UpdateOnePt(float vJ,bool& hasChanged,BitStore& onSameEdge,float eKdist) {
	float dist=getNetworkDist(vJ,onSameEdge,eKdist);
	if (dist<bestdist) {	// for S-dist
		hasChanged=true;
		processDQ(dQ_graph,bestdist,dist);
		bestcompdist.assign(compdist.begin(),compdist.end());
	}
}

inline float getOptPos(int NodeID,int NewNodeID, float eKdist) {
	float vJ,bestPos,bestVal,tmpPos,tmpVal;
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
		int NUM_STEPS=100;
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
			UpdateOnePt(vJ,hasChanged,onSameEdge,eKdist);
		}
	} else {
		if (isSumDist) {
			UpdateOnePt(0,hasChanged,onSameEdge,eKdist);
			UpdateOnePt(eKdist,hasChanged,onSameEdge,eKdist);
			for (int sp=0;sp<gQrySize;sp++)
				if (onSameEdge[sp])
					UpdateOnePt(gQryPts[sp].vP,hasChanged,onSameEdge,eKdist);
		} else {
			int NUM_STEPS=100;
			for (int j=0;j<NUM_STEPS;j++) {	// scan
				vJ=j*eKdist/(NUM_STEPS-1);
				UpdateOnePt(vJ,hasChanged,onSameEdge,eKdist);
			}
		}
	}
	return hasChanged;
}


// SPAH
void A_star(DISTMAP& oldDistMap,int seed,int dest,int& PopSize,int& VisitedSize) {
//	BitStore isVisited;
//	isVisited.assign(NodeNum,false);
	
	DISTMAP nodeVisited;
	
	StepQueue aQ;
	aQ.push(initSteps[seed][0]);
	aQ.push(initSteps[seed][1]);
	
	// gdist: from src to cur ; hdist: from cur to dist
	float x_dest,y_dest,x_cur,y_cur;
	int TmpAdjGrpAddr=getAdjListGrpAddr(dest);
	getFixedF(XCRD_A,Ref(x_dest),TmpAdjGrpAddr);
	getFixedF(YCRD_A,Ref(y_dest),TmpAdjGrpAddr);
	
	while (!aQ.empty()) {
		StepEvent event=aQ.top();
		aQ.pop();	PopSize++;
		int NodeID=event.node;

		if (nodeVisited.count(NodeID)>0) continue;
		nodeVisited[NodeID]=event.gdist;
		
		// those visited cannot be cached, only those on shortest path
		
//		if (oldDistMap.count(NodeID)>0&&event.gdist<oldDistMap[NodeID])
//			printf("%f %f\n",event.gdist,oldDistMap[NodeID]);
			
		if (oldDistMap.count(NodeID)==0||event.gdist<oldDistMap[NodeID])
			oldDistMap[NodeID]=event.gdist;
		
		
//		if (isVisited[NodeID]) continue;
//		isVisited[NodeID]=true;		
		
		
		VisitedSize++;
		
		if (isBoundary[NodeID]) {
			// adj. list for level 1 HiTi-graph
			FastArray<int>& CurAdjList=AdjList[NodeID];
			for (int z=0;z<CurAdjList.size();z++) {
				int NewNodeID=CurAdjList[z];
				edge* e=EdgeMap[getKey(NodeID,NewNodeID)];
				StepEvent newevent=event;	// copy ...
				newevent.node=NewNodeID;
				newevent.gdist+=e->dist;
				newevent.hdist=0;
				
				if (nodeVisited.count(NewNodeID)>0) continue;
				//if (isVisited[NewNodeID]) continue;
				
				TmpAdjGrpAddr=getAdjListGrpAddr(NewNodeID);
				getFixedF(XCRD_A,Ref(x_cur)	,TmpAdjGrpAddr);
				getFixedF(YCRD_A,Ref(y_cur)	,TmpAdjGrpAddr);
				if (isEuclidSpace) newevent.hdist=getDist(x_cur,y_cur,x_dest,y_dest);
				newevent.dist=newevent.gdist+newevent.hdist;
				
				// pathmax equation for non-monotonic heuristic ?
				if (newevent.dist<event.dist) newevent.dist=event.dist;
				aQ.push(newevent);	// propagation	
			}
		}
		
		float eKdist;
		int AdjListSize,NewNodeID,AdjGrpAddr;
		AdjGrpAddr=getAdjListGrpAddr(NodeID);
		getFixedF(SIZE_A,Ref(AdjListSize),AdjGrpAddr);	// read # entries
		for (int z=0;z<AdjListSize;z++) {
			getVarE(ADJNODE_A	,Ref(NewNodeID)	,AdjGrpAddr,z);
			getVarE(DIST_A		,Ref(eKdist)	,AdjGrpAddr,z);
			
			if (nodeVisited.count(NewNodeID)>0) continue;
			//if (isVisited[NewNodeID]) continue;
			
			if (isBoundary[NodeID]&&!isBoundary[NewNodeID]) {
				if (Partition[NewNodeID]!=Partition[dest]) continue;
			}
			
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
			aQ.push(newevent);	// propagation		
		}
		if (NodeID==dest) {
			return;	// or gdist ?
		}
	}
	//return MAX_DIST;
}

void IER_HiTi() {
	int rt_totna = 0; // total na
	int PopSize=0,VisitedSize=0;
	float eKdist;
	int AdjListSize,NewNodeID,AdjGrpAddr,PtGrpKey;
	DataStore result,queries;
	
	RTree *rt = new RTree(nodertrfn,128);
    rt->file->isReadOnly=true;	// don't really write blocks to disk
    rt->load_root();
    
    for (int s=0;s<gQrySize;s++) {
		DATA* newpt=new DATA;
		newpt->data[0]=M_xcrd[s];		newpt->data[1]=M_ycrd[s];
		queries.push_back(newpt);
	}
	
    DATA AggNNpt;
    BranchQueue AggNN_queue;
	{BranchEvent event;
	event.pnode=NULL;
	event.son=-1;
	event.mindist=0;	// always search from root node
	AggNN_queue.push(event);
	}
	
	while (true) {
		float mindist=getNextAggNN(rt,AggNN_queue,rt_totna,queries,&AggNNpt);
		if (mindist>=bestdist) break;	// ??? using mbr(node and adj. edges) in order to be correct
		int NodeID=AggNNpt.id;
		
		//printf("%d: %f %f %d %d\n",AggNNpt.id,mindist,bestdist,PopSize,VisitedSize);
		AdjGrpAddr=getAdjListGrpAddr(NodeID);
		getFixedF(SIZE_A,Ref(AdjListSize),AdjGrpAddr);	// read # entries
		
		for (int s=0;s<gQrySize;s++) {
			DISTMAP& curDistMap=DistMaps[s];
			if (curDistMap.count(NodeID)==0) {
				A_star(curDistMap,s,NodeID,PopSize,VisitedSize);
			}
		}
		
		for (int z=0;z<AdjListSize;z++) {
			getVarE(ADJNODE_A	,Ref(NewNodeID)	,AdjGrpAddr,z);
			getVarE(DIST_A		,Ref(eKdist)	,AdjGrpAddr,z);
			getVarE(PTKEY_A		,Ref(PtGrpKey)	,AdjGrpAddr,z);
			if (PtGrpKey>=0||!isNNquery) {	// valid PtGrpKey  (-1 for invalid key)
				for (int s=0;s<gQrySize;s++) {
					DISTMAP& curDistMap=DistMaps[s];
					if (curDistMap.count(NewNodeID)==0) {
						A_star(curDistMap,s,NewNodeID,PopSize,VisitedSize);
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
	
	int tcnt=0;
	for (int s=0;s<gQrySize;s++) tcnt+=DistMaps[s].size();
	
	printf("bestdist: %f, storesize: %d\n",bestdist,tcnt);
	printf("pop size: %d, VisitedSize: %d\n",PopSize,VisitedSize);
	printf("page faults: %d, totna: %d\n",rt->file->page_faults,rt_totna);
}

inline void printSetting() {
	if (isNNquery) printf(",point"); else printf(",center");
	if (isSumDist) printf(",sum"); else printf(",max");
}

// Only find the necessary dist, not necessarily all pairs dist
void FindSoln() {	// "node" reused as the next nodeID instead !!!
	StepEvent stepL,stepR;
	
	// initialization
	bestdist=MAX_DIST;
	while (!dQ_graph.empty()) dQ_graph.pop();	// clear kNN-dist heap
	
	onSameEdge.assign(gQrySize,false);
	partdists=new (float*)[gQrySize];
	initSteps=new (StepEvent*)[gQrySize];
	for (int i=0;i<gQrySize;i++) {
		partdists[i]=new float[2];
		initSteps[i]=new StepEvent[2];
	}
	compdist.assign(gQrySize,MAX_DIST);	
	bestcompdist.assign(gQrySize,MAX_DIST);
	
	float vP,eDist,eKdist,x1,y1,x2,y2,x,y;
	int AdjListSize,Ni,Nj,NewNodeID,TmpAdjGrpAddr;
	for (int s=0;s<gQrySize;s++) {
		vP=gQryPts[s].vP;	Ni=gQryPts[s].Ni;	Nj=gQryPts[s].Nj;
		
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
		
		stepL.xc=x1;	stepL.yc=y1;
		stepR.xc=x2;	stepR.yc=y2;
		
		stepL.gdist=vP;			stepL.hdist=0;
		stepR.gdist=eDist-vP;	stepR.hdist=0;
		
		initSteps[s][0]=stepL;
		initSteps[s][1]=stepR;
	}
	
	for (int sp=0;sp<gQrySize;sp++) 
		DistMaps[sp].clear(); // clear first (safe before use)
	
	printf("\n[IER-HiTi");
	printSetting();
	printf("]\n");
	IER_HiTi();
}

int main(int argc,char** argv) {
	ConfigType cr;
	AddConfigFromFile(cr,"config.prop");
	AddConfigFromCmdLine(cr,argc,argv);
	ListConfig(cr);
	
	string filename=getConfigStr("DBASE",cr);
	filename=filename+getConfigStr("data",cr);
	NNnum=getConfigInt("NNnum",cr);
	const char* fileprefix=filename.c_str();
  	sprintf(nodertrfn,getConfigStr("nodertr",cr));
  	
  	char queryf[255];
  	sprintf(queryf,"query/%s.qry",getConfigStr("query",cr));
  	ReadQueryFile(queryf);
  	  	
	InitClock();	// side effect: initialize the seeds for 2 random generators
	OpenDiskComm(fileprefix,DEFAULT_CACHESIZE);
	ReadHiTiDist(getConfigStr("rhiti",cr));		// don't count time
	InitClock();	// don't count file reading time
	
	// if seg. fault, check binary file type
	//isNNquery=false;	// center query
	//isSumDist=false;	// d_{max}
	DistMaps=new DISTMAP[gQrySize];
	
	isWeightUse=false;
	RefreshCache();
	FindSoln();
	CloseDiskComm();
	PrintElapsed();
	printf("\n");
	return 0;
}

/*void HiTiSearch(int src,int dest) {
	int PopSize=0,VisitedSize=0;
	BitStore isVisited;
	isVisited.assign(NodeNum,false);
	
	StepQueue aQ;
	// gdist: from src to cur ; hdist: from cur to dist
	
	StepEvent initEvt;
	initEvt.node=src;
	initEvt.gdist=0;	initEvt.hdist=0;
	
	float x_dest,y_dest,x_cur,y_cur;
	int TmpAdjGrpAddr=getAdjListGrpAddr(dest);
	getFixedF(XCRD_A,Ref(x_dest),TmpAdjGrpAddr);
	getFixedF(YCRD_A,Ref(y_dest),TmpAdjGrpAddr);
	
	TmpAdjGrpAddr=getAdjListGrpAddr(src);
	getFixedF(XCRD_A,Ref(x_cur)	,TmpAdjGrpAddr);
	getFixedF(YCRD_A,Ref(y_cur)	,TmpAdjGrpAddr);
	
	if (isEuclidSpace) initEvt.hdist=getDist(x_cur,y_cur,x_dest,y_dest);
	initEvt.dist=initEvt.gdist+initEvt.hdist;
	printf("dist: %f\n",initEvt.dist);
	aQ.push(initEvt);
	
	while (!aQ.empty()) {
		StepEvent event=aQ.top();
		aQ.pop();
		PopSize++;
		
		int NodeID=event.node;
		if (isVisited[NodeID]) continue;
		isVisited[NodeID]=true;		VisitedSize++;
		//printNodeVis(NodeID);
		// only gdist corresponds to the real dis. from src. to current node !!!
		
		//printf("%d\t%f\t%f\t%f\n",NodeID,event.dist,event.gdist,event.hdist);
		if (isBoundary[NodeID]) {
			// adj. list for level 1 HiTi-graph
			FastArray<int>& CurAdjList=AdjList[NodeID];
			for (int z=0;z<CurAdjList.size();z++) {
				int NewNodeID=CurAdjList[z];
				if (isVisited[NewNodeID]) continue;
				
				edge* e=EdgeMap[getKey(NodeID,NewNodeID)];
				StepEvent newevent=event;	// copy ...
				newevent.node=NewNodeID;
				newevent.gdist+=e->dist;
				newevent.hdist=0;
				
				TmpAdjGrpAddr=getAdjListGrpAddr(NewNodeID);
				getFixedF(XCRD_A,Ref(x_cur)	,TmpAdjGrpAddr);
				getFixedF(YCRD_A,Ref(y_cur)	,TmpAdjGrpAddr);
				if (isEuclidSpace) newevent.hdist=getDist(x_cur,y_cur,x_dest,y_dest);
				newevent.dist=newevent.gdist+newevent.hdist;
				
				// pathmax equation for non-monotonic heuristic ?
				if (newevent.dist<event.dist) newevent.dist=event.dist;
				aQ.push(newevent);	// propagation	
			}
		}
		
		int AdjListSize,NewNodeID,AdjGrpAddr;
		float eKdist;
		AdjGrpAddr=getAdjListGrpAddr(NodeID);
		getFixedF(SIZE_A,Ref(AdjListSize),AdjGrpAddr);	// read # entries
		for (int z=0;z<AdjListSize;z++) {
			getVarE(ADJNODE_A	,Ref(NewNodeID)	,AdjGrpAddr,z);
			getVarE(DIST_A		,Ref(eKdist)	,AdjGrpAddr,z);
			if (isVisited[NewNodeID]) continue;
			
			if (isBoundary[NodeID]&&!isBoundary[NewNodeID]) {
				if (Partition[NewNodeID]!=Partition[dest]) continue;
			}
			
			StepEvent newevent=event;	// copy ...
			newevent.node=NewNodeID;
			newevent.gdist+=eKdist;
			newevent.hdist=0;
			
			TmpAdjGrpAddr=getAdjListGrpAddr(NewNodeID);
			getFixedF(XCRD_A,Ref(x_cur)	,TmpAdjGrpAddr);
			getFixedF(YCRD_A,Ref(y_cur)	,TmpAdjGrpAddr);
			
			if (isEuclidSpace) newevent.hdist=getDist(x_cur,y_cur,x_dest,y_dest);
			newevent.dist=newevent.gdist+newevent.hdist;
			
			// pathmax equation for non-monotonic heuristic ?
			if (newevent.dist<event.dist) newevent.dist=event.dist;
			aQ.push(newevent);	// propagation
		}
		
		if (NodeID==dest) {
			printf("dist: %0.3f (%0.3f + %0.3f) \n",event.dist,event.gdist,event.hdist);
			// placed at the end of while loop so as to collect back all the entries
			break;
		}
	}
	printf("pop size: %d, VisitedSize: %d\n",PopSize,VisitedSize);
	//printQueueVis(aQ,isVisited);
}

int HiTimain(int argc,char** argv) {
	ConfigType cr;
	AddConfigFromFile(cr,"config.prop");
	AddConfigFromCmdLine(cr,argc,argv);
	ListConfig(cr);
	
	string filename=getConfigStr("DBASE",cr);
	filename=filename+getConfigStr("data",cr);
	const char* fileprefix=filename.c_str();
  	
	InitClock();	// side effect: initialize the seeds for 2 random generators
	OpenDiskComm(fileprefix,DEFAULT_CACHESIZE);
	
	ReadHiTiDist(getConfigStr("rhiti",cr));
	
	InitClock();	// don't count file reading time
	
	int src=0;
	int dest=2445;	// for C (OL_dnw)
	
	isEuclidSpace=true;		// using R^2 space or not
	//isEuclidSpace=false;
	
	//isPrintVis=true;
	
	HiTiSearch(src,dest);
	
	CloseDiskComm();
	
	PrintElapsed();
	return 0;
}
*/

