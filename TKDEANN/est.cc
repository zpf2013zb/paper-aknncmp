#include "utility.h"
#include "netshare.h"
#include "diskbased.h"
#include "matrix.h"

#define MAXREAL 1.0e20
typedef map<int,float> DISTMAP;

const int nterm=7;

struct smp_rec {
	float x,y;
	int NodeID;
	int cnt;
	float a,b;	// a: log intercept, b: log slope
	float poly_coeff[nterm];
};

smp_rec* sample;
int sampleCnt=0;


DISTMAP dm;	// only use cheap storage for range search
int* pcntm;
bool* isSample;
vector<float> pt_xcrd,pt_ycrd;

float getDist(float x1,float y1,float x2,float y2) {
	return float(sqrt((x2-x1)*(x2-x1)+(y2-y1)*(y2-y1)));
}

void poly_fit(int npt,double* xvec,double* yvec,float* coeff) {
	matrix y(npt,1);
	matrix X(npt,nterm);
	
	for (int i=0;i<npt;i++) {
		y(i+1,1)=yvec[i];
		
		double tval=1.0;
		for (int j=0;j<nterm;j++) {
			X(i+1,j+1)=tval;
			tval*=xvec[i];
		}
	}
	
	//y.print(); printf("\n");
	//X.print(); printf("\n");
	
	matrix XT=X.transpose();
	matrix a= (XT*X).inverse()*XT*y;
	//a.print();
	
	for (int j=0;j<nterm;j++) 
		coeff[j]=a(j+1,1);
}

float getPolyEst(float range,float* coeff) {
	double answer=0;
	double tval=1.0;
	for (int j=0;j<nterm;j++) {
		answer+=coeff[j]*tval;
		tval*= range/1000.0; // temp norm factor
	}
	if (answer<0) return 0;
	return (float)answer;
}

int RangeCount(int source_node,float radius,smp_rec* anchor=NULL) {
	StepEvent rootstep;
	static vector<float> tempd;
	int count=0;
	int MaxHeapSize=0,PopSize=0;
	int AdjListSize,NewNodeID;
	float eKdist,xVal,yVal;
	static StepQueue sQ;
	float startx,starty;
	bool isFirst=true;
	
	rootstep.ClusID=0;
	rootstep.dist=0;	
	rootstep.node=source_node;	// Ni=0
	rootstep.isSortAcc=true;
	while (sQ.size()>0) sQ.pop();	// clear sQ
	sQ.push(rootstep);	
	dm.clear();		// clear dm
	tempd.clear();
	
	while (!sQ.empty()) {
		if (MaxHeapSize<sQ.size()) MaxHeapSize=sQ.size();
		StepEvent event=sQ.top();
		sQ.pop();	PopSize++;
		if (event.dist>radius) break;		// terminate
		
		int NodeID=event.node;
		int id=event.ClusID;
		if (dm.count(NodeID)>0) continue;	// already found
		dm[NodeID]=event.dist;
		
		AdjListSize=AdjList[NodeID].size();	// read # entries
		xVal=xcrd[NodeID];
		yVal=ycrd[NodeID];
		
		if (isFirst) {
			isFirst=false;
			startx=xVal;
			starty=yVal;
		}
		
		for (int z=0;z<AdjListSize;z++) {
			NewNodeID=AdjList[NodeID][z];
			edge* e=EdgeMap[getKey(NodeID,NewNodeID)];
			eKdist=e->dist;
			
			StepEvent newevent=event;	// copy ...
			newevent.node=NewNodeID;
			newevent.dist=event.dist+eKdist;
			if (dm.count(NewNodeID)>0) continue;
			sQ.push(newevent);	// propagation
			
			for (int j=0;j<e->pts.size();j++) {	// scan
				float vJ;
				vJ=e->pts[j];
				if (event.dist+vJ<=radius) { // a bit incorrect !
					count++;
					if (anchor!=NULL)
						tempd.push_back(event.dist+vJ);
					float xpt=xcrd[NodeID]+(xcrd[NewNodeID]-xcrd[NodeID])*vJ/eKdist;
					float ypt=ycrd[NodeID]+(ycrd[NewNodeID]-ycrd[NodeID])*vJ/eKdist;
					float tempdist=getDist(startx,starty,xpt,ypt);
					
					//printf("%f\n",tempdist	);
				}
			}
		}
	}
	//printf("count %d | (ms,ps)=(%d,%d)\n",count,MaxHeapSize,PopSize);
	if (anchor!=NULL) {
		sort(tempd.begin(),tempd.end());
		
		int max_sz=200;
		int npt=1;
		double tempx[max_sz],tempy[max_sz];
		tempx[0]=0;	// first point
		tempy[0]=0;
		
		float sumx=0,sumx2=0,sumy=0,sumxy=0;
		int tnum=0;
		float lastx=0;
		const float incr=0.1;
		for (int i=0;i<tempd.size();i++) {
			float y=log( (float)(i+1) );
			float x=log( tempd[i] );
			if (x<=0||y<=0) {	// 
				//printf("*** %f %f\n",x,y);
				continue;
			}
			if (lastx!=0&&x<incr+lastx&&i+1!=tempd.size())
				continue;
			
			lastx=x;	// required for fair line fitting3
			sumx+=x;
			sumx2+=x*x;
			sumy+=y;
			sumxy+=x*y;
			//printf("%f %f\n",x,y);
			tnum++;
			
			if (npt<max_sz) { // avoid same point or just take uniform ...
				tempx[npt]= tempd[i]/1000.0;	// temp norm. factor
				tempy[npt]=i+1;
				npt++;
			}
		}
		float a=(sumy*sumx2-sumx*sumxy)/(tnum*sumx2-sumx*sumx);
		float b=(tnum*sumxy-sumx*sumy)/(tnum*sumx2-sumx*sumx);
		//printf("%f %f %d\n",a,b,count);
		anchor->a=a;
		anchor->b=b;
		
		poly_fit(npt,tempx,tempy,anchor->poly_coeff);
		//exit(0);
	}
	return count;
}

int getRandNID() {
	int Ni,Nj,dummy;
	int PtGrpKey=rand()%PtNum;
	int PtGrpAddr=pointQuery(PtTree,PtGrpKey,dummy);
	getFixedF(NI_P,Ref(Ni),PtGrpAddr);
	getFixedF(NJ_P,Ref(Nj),PtGrpAddr);
	
	// point based (biased) samples, distribution follows the points
	return ((rand()%2==0)?Ni:Nj);
	
	// node based samples
//	return (rand()%NodeNum);
}

// use two techniques to choose samples, one random, one choose far awayy

void ChooseSamples() {
	for (int i=0;i<sampleCnt;i++) {
		sample[i].NodeID=getRandNID();
		sample[i].x=xcrd[sample[i].NodeID];
		sample[i].y=ycrd[sample[i].NodeID];
		
		if (i>0) {
			int gNID=sample[i].NodeID;
			float gMinDist=0,gSumDist=0;
			
			for (int round=0;round<10;round++) {
				float curx,cury;
				int nid=getRandNID();
				curx=xcrd[nid];
				cury=ycrd[nid];
				
				float mdist=MAXREAL;
				for (int z=0;z<i;z++) {
					float tdist=getDist(curx,cury,sample[z].x,sample[z].y);
					if (tdist<mdist) mdist=tdist;
				}
								
				if (mdist>gMinDist) {	// i.e. mindist to others quite large !
					gMinDist=mdist;
					gNID=nid;
				}
			}
			
			// mindist to other furthest !
			sample[i].NodeID=gNID;
			sample[i].x=xcrd[sample[i].NodeID];
			sample[i].y=ycrd[sample[i].NodeID];
		}
		
		//printf("anchor: %d\n",pcntm[sample[i].NodeID]);
		isSample[sample[i].NodeID]=true;
	}
}

// new technique: build ref. point on-the-fly, useful when |S| >> 100

void test(int _smpCnt,float radius,float anchor_rad) {
	isSample=new bool[NodeNum];
	for (int i=0;i<NodeNum;i++) {
		isSample[i]=false;
	}
	
	/*for (int i=0;i<NodeNum;i++)
		if (pcntm[i]>0)
			printf("%d %d\n",i,pcntm[i]);*/
	
	vector<int> nQ;
	int qCnt=100;
	for (int z=0;z<qCnt;z++) {			// init queries
		int nid=getRandNID();
		nQ.push_back(nid);
	}
	
	// gen. samples
	sampleCnt=_smpCnt;
	sample=new smp_rec[sampleCnt];
	
	printf("\nsamples\n");
	ChooseSamples();
	for (int i=0;i<sampleCnt;i++) 
		sample[i].cnt=-1;	// on-the-fly testing 28/6/2005 18:49
	
	printf("\nqueries\n");
	float gSumErr_LPL=0,gSumErr_POLY=0,gSumErr_OPT=0,gSumAct=0;
	for (int z=0;z<qCnt;z++) {
		float curx,cury;
		int AdjGrpAddr=getAdjListGrpAddr(nQ[z]);
		curx=xcrd[nQ[z]];
		cury=ycrd[nQ[z]];
		
		int curcnt=RangeCount(nQ[z],radius);
		int est_lpl=0,est_poly=0,est_opt;
		
		float dist=MAXREAL;
		int id=-1;
		{
			for (int i=0;i<sampleCnt;i++) {
				float tdist=getDist(curx,cury,sample[i].x,sample[i].y);
				if (tdist<dist) {
					dist=tdist;
					id=i;
				}
			}
			if (sample[id].cnt<0) { // on-the-fly
				sample[id].cnt=0;
				sample[id].cnt=RangeCount(sample[id].NodeID,anchor_rad,&(sample[id]));
				if (anchor_rad!=radius)
					sample[id].cnt=RangeCount(sample[id].NodeID,radius);
			}
			
			est_opt=sample[id].cnt;  // using nearest witness
			est_poly=(int)getPolyEst(radius,sample[id].poly_coeff);
			
			float a=sample[id].a;
			float b=sample[id].b;
			est_lpl=(int)exp(a+b*log(radius));	
		}
		printf("%d: %d, %d %d %d (%f)\n",z,curcnt,est_lpl,est_poly,est_opt,dist);
		
		gSumErr_LPL+=ValAbs(est_lpl-curcnt);
		gSumErr_POLY+=ValAbs(est_poly-curcnt);
		gSumErr_OPT+=ValAbs(est_opt-curcnt);
		gSumAct+=curcnt;
	}
	
	printf("error (LPL POLY OPT): %0.4f %0.4f %0.4f\n",
			gSumErr_LPL/gSumAct,gSumErr_POLY/gSumAct,gSumErr_OPT/gSumAct);	// error defined in [APR99]
}


void ReadAllFiles() {
	// read point groups
	float pt_pos,edge_dist,xVal,yVal;
	int Ni,Nj,grp_sz;
	int pt_cnt=0;
	int AdjListSize,NewNodeID,AdjGrpAddr;
	
	xcrd.clear();
	ycrd.clear();
	assert(AdjList==NULL);

	AdjList=new FastArray<int>[NodeNum];	
	for (int NodeID=0;NodeID<NodeNum;NodeID++) {
		AdjGrpAddr=getAdjListGrpAddr(NodeID);
		getFixedF(XCRD_A,Ref(xVal),AdjGrpAddr);
		getFixedF(YCRD_A,Ref(yVal),AdjGrpAddr);
		xcrd.push_back(xVal);
		ycrd.push_back(yVal);
						
		getFixedF(SIZE_A,Ref(AdjListSize),AdjGrpAddr);	// read # entries
		for (int z=0;z<AdjListSize;z++) {
			getVarE(ADJNODE_A	,Ref(NewNodeID)	,AdjGrpAddr,z);
			getVarE(DIST_A		,Ref(edge_dist)	,AdjGrpAddr,z);
			AdjList[NodeID].push_back(NewNodeID);
			
			int em_key=getKey(NodeID,NewNodeID);
			if (EdgeMap.count(em_key)==0) {
				edge* e=new edge;
				EdgeMap[em_key]=e;
				e->Ni=min(NodeID,NewNodeID);	
				e->Nj=max(NodeID,NewNodeID);	// note that: Ni<Nj
				e->dist=edge_dist;
			}
		}
	}
	
	fseek(PtFile,0,SEEK_SET);
	while (!feof(PtFile)&&pt_cnt<PtNum) {
		fread(&(Ni),1,sizeof(int),PtFile);
		fread(&(Nj),1,sizeof(int),PtFile);
		fread(&(edge_dist),1,sizeof(float),PtFile);
		fread(&(grp_sz),1,sizeof(int),PtFile);
		
		// no repeated edge
		edge* e=EdgeMap[getKey(Ni,Nj)]; // already constructed
		for (int i=0;i<grp_sz;i++) {
			fread(&(pt_pos),1,sizeof(float),PtFile); 
			float xpt=xcrd[Ni]+(xcrd[Nj]-xcrd[Ni])*pt_pos/edge_dist;
			float ypt=ycrd[Ni]+(ycrd[Nj]-ycrd[Ni])*pt_pos/edge_dist;
			pt_xcrd.push_back(xpt);
			pt_ycrd.push_back(ypt);
			
			e->pts.push_back(pt_pos);
			if (pt_pos<edge_dist/2.0)
				pcntm[e->Ni]++;
			else 
				pcntm[e->Nj]++;
		}
		pt_cnt+=grp_sz;
	}
	printf("check cnt: %d\n",pt_cnt);
}

void CleanTempStruct() {
	for (int NodeID=0;NodeID<NodeNum;NodeID++) 
		AdjList[NodeID].clear();
	delete[] AdjList;
	
	EdgeMapType::iterator iter=EdgeMap.begin();
 	while (iter!=EdgeMap.end()) {
 		edge* e=iter->second;
 		e->pts.clear();
 		delete e;
 		iter++;
	}
	EdgeMap.clear();
}



int main(int argc,char** argv) {
	ConfigType cr;
	AddConfigFromFile(cr,"config.prop");
	AddConfigFromCmdLine(cr,argc,argv);
	//ListConfig(cr);
	
	string filename=getConfigStr("DBASE",cr);
	filename=filename+getConfigStr("data",cr);
	NNnum=getConfigInt("NNnum",cr);
	const char* fileprefix=filename.c_str();
  	
  	
  	OpenDiskComm(fileprefix,DEFAULT_CACHESIZE);
	InitClock();	// side effect: initialize the seeds for 2 random generators
	srand(0);	srand48(0);	// set ..
	
	RefreshCache();
	pcntm=new int[NodeNum];
	for (int i=0;i<NodeNum;i++)
		pcntm[i]=0;
	ReadAllFiles();
	
	int _smpCnt=getConfigInt("sc",cr);	//1000
	float _radius=getConfigFloat("e",cr);	//100
	float _nhood=getConfigFloat("nhood",cr);
	test(_smpCnt,_radius,_nhood);
	
	CleanTempStruct();
	CloseDiskComm();
	PrintElapsed();
	
	printf("\n");
	return 0;
}
