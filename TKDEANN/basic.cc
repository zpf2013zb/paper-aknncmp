#include "utility.h"
#include "netshare.h"
#include "diskbased.h"

typedef map<int,float> DISTMAP;
DISTMAP* DistMaps;	// only use cheap storage for range search

BitStore onSameEdge;
float** partdists;
FastArray<float> compdist,bestcompdist;
float bestdist;		// best value found so far

bool isEuclidSpace=true;
bool isNNquery=true;

inline void UpdateOnePt(float vJ,bool& hasChanged,BitStore& onSameEdge,float eKdist) {
	float dist=0,tmpnmdist=0;
	for (int sp=0;sp<gQrySize;sp++) {
		float distL=partdists[sp][0]+vJ;
		float distR=partdists[sp][1]+eKdist-vJ;
		float distT=min(distL,distR);
		if (onSameEdge[sp]) {
			float distQ=ValAbs(gQryPts[sp].vP-vJ);
			if (distQ<distT) distT=distQ;
		}
		dist+=distT;
		compdist[sp]=distT;	// for log
		tmpnmdist=max(distT,tmpnmdist);
	}
	
	if (isSumDist) {
		if (dist<bestdist) {	// for S-dist
			hasChanged=true;
			bestdist=dist;
			bestcompdist.assign(compdist.begin(),compdist.end());
		}
	} else {
		if (tmpnmdist<bestdist) {	// for NM-dist
			hasChanged=true;
			bestdist=tmpnmdist;
			bestcompdist.assign(compdist.begin(),compdist.end());
		}
	}
}

inline void LinearUpdate(bool& hasChanged,BitStore& onSameEdge,float eKdist) {
	UpdateOnePt(0,hasChanged,onSameEdge,eKdist);
	UpdateOnePt(eKdist,hasChanged,onSameEdge,eKdist);
	for (int sp=0;sp<gQrySize;sp++) 
		if (onSameEdge[sp]) 
			UpdateOnePt(gQryPts[sp].vP,hasChanged,onSameEdge,eKdist);
}

inline void OptUpdate(bool& hasChanged,BitStore& onSameEdge,float eKdist) {
	const int NUM_STEP=1000;
	if (isSumDist)
		LinearUpdate(hasChanged,onSameEdge,eKdist);
	else {
		for (int j=0;j<NUM_STEP;j++) {	// scan
			float vJ=j*eKdist/(NUM_STEP-1);
			UpdateOnePt(vJ,hasChanged,onSameEdge,eKdist);
		}
	}
	// later derive soln. for max. update
}

bool UpdatePoints(int NodeID,int NewNodeID,int PtGrpKey, float eKdist) {
	float vJ;
	int PtGrpAddr,PtGrpSize,dummy;
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
		PtGrpAddr=pointQuery(PtTree,PtGrpKey,dummy);
		getFixedF(SIZE_P,Ref(PtGrpSize),PtGrpAddr);
	 	for (int j=0;j<PtGrpSize;j++) {	// scan
			getVarE(PT_P,Ref(vJ),PtGrpAddr,j);
			UpdateOnePt(vJ,hasChanged,onSameEdge,eKdist);
			//printf("%f %f\n",vJ,eKdist);
		}
	} else 
		OptUpdate(hasChanged,onSameEdge,eKdist);
	return hasChanged;
}

// debug later !!!
void BruteExpansion(StepQueue& sQ) {
	int AdjListSize,NewNodeID,AdjGrpAddr;
	float eKdist;
	int MaxHeapSize=0,AggPopSize=0;	
	
	FastArray<float> NodeDists[gQrySize];
	for (int id=0;id<gQrySize;id++) 
		NodeDists[id].assign(NodeNum,MAX_DIST);
	
	float maxDistVal=0;
	while (!sQ.empty()) {
		if (MaxHeapSize<sQ.size()) MaxHeapSize=sQ.size();
		StepEvent event=sQ.top();
		sQ.pop();
		AggPopSize++;
		
		int NodeID=event.node;
		int id=event.ClusID;
				
		if (NodeDists[id][NodeID]<MAX_DIST) continue;	// visited
		NodeDists[id][NodeID]=event.dist;
		DistMaps[id][NodeID]=event.dist;
		maxDistVal=max(event.dist,maxDistVal);
		
		AdjGrpAddr=getAdjListGrpAddr(NodeID);
		getFixedF(SIZE_A,Ref(AdjListSize),AdjGrpAddr);	// read # entries
		for (int z=0;z<AdjListSize;z++) {
			// read entry details
			getVarE(ADJNODE_A	,Ref(NewNodeID)	,AdjGrpAddr,z);
			getVarE(DIST_A		,Ref(eKdist)	,AdjGrpAddr,z);
			StepEvent newevent=event;	// copy ...
			newevent.node=NewNodeID;
			newevent.dist=event.dist+eKdist;
			if (NodeDists[id][NewNodeID]==MAX_DIST) sQ.push(newevent);	// propagation
		}
	}
	
	int GroupAddr=0;
	int PtGrpSize,PtGrpKey,PtGrpAddr,dummy,vJ;
	bool hasChanged=false;
	
	// scan adj. list
	for (int NodeID=0;NodeID<NodeNum;NodeID++) {
		AdjGrpAddr=getAdjListGrpAddr(NodeID);
		getFixedF(SIZE_A,Ref(AdjListSize),AdjGrpAddr);	// read # entries
		for (int z=0;z<AdjListSize;z++) {
			// read entry details
			getVarE(ADJNODE_A	,Ref(NewNodeID)	,AdjGrpAddr,z);
			getVarE(DIST_A		,Ref(eKdist)	,AdjGrpAddr,z);
			getVarE(PTKEY_A	,Ref(PtGrpKey)	,AdjGrpAddr,z);
			if (NodeID>NewNodeID) continue;	// as they are in ascending order ?
			
			// assume Ni<Nj
			if (PtGrpKey>=0||!isNNquery) // valid PtGrpKey  (-1 for invalid key)
				UpdatePoints(NodeID,NewNodeID,PtGrpKey,eKdist);
		}
	}
	printf("Heap: max_size %d, pop_size %d\n",MaxHeapSize,AggPopSize);
	printf("bestdist: %f\n",bestdist);
}

// memory requirement: r bitmap and r tempdiststore
void MultipleExpansion(StepQueue& sQ) {
	int MaxHeapSize=0,AggPopSize=0;
	while (!sQ.empty()) {
		if (MaxHeapSize<sQ.size()) MaxHeapSize=sQ.size();
		StepEvent event=sQ.top();
		sQ.pop();
		AggPopSize++;
		
		int NodeID=event.node;
		int id=event.ClusID;
		
		if (DistMaps[id].count(NodeID)>0) continue;
		DistMaps[id][NodeID]=event.dist;
		
		if (isSumDist) {
			if (event.dist>=bestdist) break;	// for S-dist
		} else {
			if (event.dist>=bestdist) break;	// for NM-dist
		}
		
		int AdjListSize,NewNodeID,AdjGrpAddr;
		float eKdist;
		AdjGrpAddr=getAdjListGrpAddr(NodeID);
		getFixedF(SIZE_A,Ref(AdjListSize),AdjGrpAddr);	// read # entries
		for (int z=0;z<AdjListSize;z++) {
			int PtGrpKey;
			
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
			if (!needUpdate) continue;
			
			getVarE(PTKEY_A	,Ref(PtGrpKey)	,AdjGrpAddr,z);
			if (PtGrpKey>=0||!isNNquery) {	// valid PtGrpKey  (-1 for invalid key)
				UpdatePoints(NodeID,NewNodeID,PtGrpKey,eKdist);
			}
		}
	}
	printf("Heap: max_size %d, pop_size %d\n",MaxHeapSize,AggPopSize);
	printf("bestdist: %f\n",bestdist);
	if (!sQ.empty()) 
		printf("topdist: %f\n",sQ.top().dist);
	int totalmapsize=0;
	for (int sp=0;sp<gQrySize;sp++) 
		totalmapsize+=DistMaps[sp].size();
	printf("totalmapsize: %d\n",totalmapsize);
}

// memory requirement: r bitmap, r tempdiststore minimized
void MAXAG_Expansion(StepQueue& sQ) {
	BitStore isVisited[gQrySize];
	for (int sp=0;sp<gQrySize;sp++) isVisited[sp].assign(NodeNum,false);
	int MaxHeapSize=0,AggPopSize=0;
	
	if (isSumDist) return;
	
	while (!sQ.empty()) {
		if (MaxHeapSize<sQ.size()) MaxHeapSize=sQ.size();
		StepEvent event=sQ.top();
		sQ.pop();
		AggPopSize++;
		
		int NodeID=event.node;
		int id=event.ClusID;
		
		if (isVisited[id][NodeID]) continue;
		isVisited[id][NodeID]=true;
		DistMaps[id][NodeID]=event.dist;
		
		if (event.dist>=bestdist) break;	// for NM-dist
		
		int AdjListSize,NewNodeID,AdjGrpAddr;
		float eKdist;
		AdjGrpAddr=getAdjListGrpAddr(NodeID);
		getFixedF(SIZE_A,Ref(AdjListSize),AdjGrpAddr);	// read # entries
		
		float maxedgedist=0;
		for (int z=0;z<AdjListSize;z++) {
			float vJ;
			int PtGrpKey,PtGrpAddr,PtGrpSize,dummy;
			
			getVarE(ADJNODE_A	,Ref(NewNodeID)	,AdjGrpAddr,z);
			getVarE(DIST_A		,Ref(eKdist)	,AdjGrpAddr,z);
			maxedgedist=max(eKdist,maxedgedist);
				
			StepEvent newevent=event;	// copy ...
			newevent.node=NewNodeID;		
			newevent.dist=event.dist+eKdist;
			if (isVisited[id][NewNodeID]==false) sQ.push(newevent);	// propagation
			
			// if error, add cond. to (... && not on same edge) [sp]
			bool needUpdate=true;
			for (int sp=0;sp<gQrySize;sp++) 
				if (isVisited[sp][NodeID]==false&&isVisited[sp][NewNodeID]==false)
					needUpdate=false;
			if (!needUpdate) continue;
			
			getVarE(PTKEY_A	,Ref(PtGrpKey)	,AdjGrpAddr,z);
			if (PtGrpKey>=0||!isNNquery) {	// valid PtGrpKey  (-1 for invalid key)
				bool hasChanged=false;
				int sid=0,eid=1;
				if (NewNodeID<NodeID) { sid=1; eid=0;}
				for (int sp=0;sp<gQrySize;sp++) {
					partdists[sp][0]=partdists[sp][1]=MAX_DIST;
					
					if (DistMaps[sp].count(NodeID)>0) 
						partdists[sp][sid]=DistMaps[sp][NodeID];
					else if (isVisited[sp][NodeID])
						partdists[sp][sid]=0;
					if (DistMaps[sp].count(NewNodeID)>0) 
						partdists[sp][eid]=DistMaps[sp][NewNodeID];
					else if (isVisited[sp][NewNodeID])
						partdists[sp][eid]=0;
					
					int qNi=gQryPts[sp].Ni,qNj=gQryPts[sp].Nj;
					onSameEdge[sp]=(NodeID==qNi&&NewNodeID==qNj)||(NodeID==qNj&&NewNodeID==qNi);
				}
				
				if (isNNquery) {
					PtGrpAddr=pointQuery(PtTree,PtGrpKey,dummy);
					getFixedF(SIZE_P,Ref(PtGrpSize),PtGrpAddr);
				 	for (int j=0;j<PtGrpSize;j++) {	// scan
						getVarE(PT_P,Ref(vJ),PtGrpAddr,j);
						UpdateOnePt(vJ,hasChanged,onSameEdge,eKdist);
					}
				} else 
					OptUpdate(hasChanged,onSameEdge,eKdist);
				if (hasChanged)	{
					printf("%d %d %f\n",NodeID,NewNodeID,eKdist);
					for (int sp=0;sp<gQrySize;sp++) {
						if (DistMaps[sp].count(NodeID)>0) 
							printf("%d 1: %f\n",sp,DistMaps[sp][NodeID]);
						if (DistMaps[sp].count(NewNodeID)>0) 
							printf("%d 2: %f\n",sp,DistMaps[sp][NewNodeID]);
					}
				}
			}
		}
		
		for (int sp=0;sp<gQrySize;sp++) {
			if (sp==id) continue;	// don't remove itself
			if (DistMaps[sp].count(NodeID)>0) {
				if (DistMaps[sp][NodeID]<event.dist-maxedgedist)
					DistMaps[sp].erase(NodeID);
			}
		}
	}
	printf("Heap: max_size %d, pop_size %d\n",MaxHeapSize,AggPopSize);
	printf("bestdist: %f\n",bestdist);
	if (!sQ.empty()) 
		printf("topdist: %f\n",sQ.top().dist);
	int totalmapsize=0;
	for (int sp=0;sp<gQrySize;sp++) 
		totalmapsize+=DistMaps[sp].size();
	printf("totalmapsize: %d\n",totalmapsize);
}


// Only find the necessary dist, not necessarily all pairs dist
void FindSoln(int Method) {	// "node" reused as the next nodeID instead !!!
	StepQueue sQ;
	StepEvent stepL,stepR;
	
	while (!sQ.empty()) sQ.pop();
	for (int s=0;s<gQrySize;s++) {
		float vP,eDist,eKdist;
		int Ni,Nj,AdjListSize,NewNodeID;
		
		vP=gQryPts[s].vP;
		Ni=gQryPts[s].Ni;
		Nj=gQryPts[s].Nj;
		
		// lookup eDist from adj. list. of Ni;
		int TmpAdjGrpAddr=getAdjListGrpAddr(Ni);
		getFixedF(SIZE_A,Ref(AdjListSize),TmpAdjGrpAddr);	// read # entries
		for (int z=0;z<AdjListSize;z++) {
			getVarE(ADJNODE_A	,Ref(NewNodeID)	,TmpAdjGrpAddr,z);
			getVarE(DIST_A		,Ref(eKdist)	,TmpAdjGrpAddr,z);
			if (NewNodeID==Nj) eDist=eKdist;
		}
		
		stepL.ClusID=stepR.ClusID=s;
		stepL.dist=vP;			stepL.node=Ni;
		stepR.dist=eDist-vP;	stepR.node=Nj;
		sQ.push(stepL);			sQ.push(stepR);	
	}
	
	// clear temp store, safe before use
	for (int sp=0;sp<gQrySize;sp++) DistMaps[sp].clear();
	
	// initialize variables
	bestdist=MAX_DIST;
	onSameEdge.assign(gQrySize,false);
	partdists=new (float*)[gQrySize];
	for (int i=0;i<gQrySize;i++) partdists[i]=new float[2];
	compdist.assign(gQrySize,MAX_DIST);
	bestcompdist.assign(gQrySize,MAX_DIST);
	
	if (Method==0) {
		printf("\n[Brute]\n");
		BruteExpansion(sQ);
	}
	if (Method==1) {
		printf("\n[Multiple]\n");
		MultipleExpansion(sQ);
	}
	if (Method==2) {
		isSumDist=false;
		printf("\n[MAXAG]\n");
		MAXAG_Expansion(sQ);
	}
	for (int s=0;s<gQrySize;s++)
		printf("%d: %f\n",s,bestcompdist[s]);
}

void GNN_Start() {
	DistMaps=new DISTMAP[gQrySize];
	FindSoln(0);
	//FindSoln(1);
	//FindSoln(2);
}

//typedef	priority_queue<BranchEvent,vector<BranchEvent>,BranchComp> BranchQueue;
// error for 10000 times: 2.804688
int main(int argc,char** argv) {
	ConfigType cr;
	AddConfigFromFile(cr,"config.prop");
	AddConfigFromCmdLine(cr,argc,argv);
	ListConfig(cr);
	
	string filename=getConfigStr("DBASE",cr);
	filename=filename+getConfigStr("data",cr);
	const char* fileprefix=filename.c_str();
	
	char queryf[255];
	sprintf(queryf,"query/%s.qry",getConfigStr("query",cr));
	ReadQueryFile(queryf);
  	 	
	InitClock();	// side effect: initialize the seeds for 2 random generators
	OpenDiskComm(fileprefix,DEFAULT_CACHESIZE);	
	InitClock();	// don't count file reading time
	
	//isSumDist=false;	// d_{max}
	//isNNquery=false;
	
	GNN_Start();
	
//	for (int j=0;j<10;j++) {
//		float delta=drand48();//4.956043;
//		float sum=0;
//		for (int i=0;i<8301;i++)
//			sum+=delta;
//		printf("%f %f %f\n",sum-8301*delta,sum,8301*delta);
//	}

	CloseDiskComm();	
	PrintElapsed();
	return 0;
}
