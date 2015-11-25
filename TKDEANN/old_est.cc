#include "utility.h"
#include "netshare.h"
#include "diskbased.h"

#define MAXREAL 1.0e20
typedef map<int,float> DISTMAP;

struct smp_rec {
	float x,y;
	int NodeID;
	int cnt;
	float eps,a,b;	// a: log intercept, b: log slope
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

int EuclidCnt(float curx,float cury,float radius) {
	int count=0;
	for (int i=0;i<PtNum;i++) 
		if (getDist(curx,cury,pt_xcrd[i],pt_ycrd[i])<=radius)
			count++;
	return count;
}

int RangeCount(StepQueue& sQ,float radius,smp_rec* anchor=NULL) {
	int count=0;
	int MaxHeapSize=0,PopSize=0;
	int AdjListSize,NewNodeID,AdjGrpAddr,PtGrpKey;
	float eKdist,xVal,yVal;
	
	dm.clear();		// clear dm
	
	float startx,starty;
	bool isFirst=true;
	float gEps=0;
	
	static vector<float> tempd;
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
		
		AdjGrpAddr=getAdjListGrpAddr(NodeID);
		getFixedF(SIZE_A,Ref(AdjListSize),AdjGrpAddr);	// read # entries
		xVal=xcrd[NodeID];
		yVal=ycrd[NodeID];
		
		if (isFirst) {
			isFirst=false;
			startx=xVal;
			starty=yVal;
		}
		
		for (int z=0;z<AdjListSize;z++) {
			getVarE(ADJNODE_A	,Ref(NewNodeID)	,AdjGrpAddr,z);
			getVarE(DIST_A		,Ref(eKdist)	,AdjGrpAddr,z);
			
			StepEvent newevent=event;	// copy ...
			newevent.node=NewNodeID;
			newevent.dist=event.dist+eKdist;
			if (dm.count(NewNodeID)>0) continue;
			sQ.push(newevent);	// propagation
			
			getVarE(PTKEY_A	,Ref(PtGrpKey)	,AdjGrpAddr,z);
			if (PtGrpKey>=0) {	// valid PtGrpKey  (-1 for invalid key)
				float vJ;
				int PtGrpAddr,PtGrpSize,dummy;
				
				PtGrpAddr=pointQuery(PtTree,PtGrpKey,dummy);
				getFixedF(SIZE_P,Ref(PtGrpSize),PtGrpAddr);
	 			for (int j=0;j<PtGrpSize;j++) {	// scan
					getVarE(PT_P,Ref(vJ),PtGrpAddr,j);	// PtGrpKey+j
					if (event.dist+vJ<=radius) {
						count++;
						if (anchor!=NULL)
							tempd.push_back(event.dist+vJ);
						float xpt=xcrd[NodeID]+(xcrd[NewNodeID]-xcrd[NodeID])*vJ/eKdist;
						float ypt=ycrd[NodeID]+(ycrd[NewNodeID]-ycrd[NodeID])*vJ/eKdist;
						float tempdist=getDist(startx,starty,xpt,ypt);
						gEps=max(gEps,tempdist);	// approx. 
						
						//printf("%f\n",tempdist	);
					}
				}
			}
		}
	}
	//printf("count %d | (ms,ps)=(%d,%d)\n",count,MaxHeapSize,PopSize);
	if (anchor!=NULL) {
		sort(tempd.begin(),tempd.end());
		
		float sumx=0,sumx2=0,sumy=0,sumxy=0;
		int tnum=0;
		float lastx=0;
		const float incr=0.05;
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
			printf("%f %f %f\n",tempd[i],x,y);
			tnum++;
		}
		//printf("\n");
		float a=(sumy*sumx2-sumx*sumxy)/(tnum*sumx2-sumx*sumx);
		float b=(tnum*sumxy-sumx*sumy)/(tnum*sumx2-sumx*sumx);
		printf("%f %f %d | %f %f %f %f\n",a,b,count,exp(a+b*log(100.0)),
				exp(a+b*log(500.0)),exp(a+b*log(1000.0)),exp(a+b*log(10000.0)));
		
		// find real network and real points !
		//gEps=gEps/count;
		
		anchor->a=a;
		anchor->b=b;
		anchor->eps=gEps;
	}
	return count;
}

void CountExpansion(StepQueue& sQ) {
	int MaxHeapSize=0,AggPopSize=0;
	int AdjListSize,NewNodeID,AdjGrpAddr,PtGrpKey;
	float eKdist;
	
	dm.clear();		// clear dm
	while (!sQ.empty()) {
		if (MaxHeapSize<sQ.size()) MaxHeapSize=sQ.size();
		StepEvent event=sQ.top();
		sQ.pop();	AggPopSize++;
		
		int NodeID=event.node;
		int id=event.ClusID;
		if (dm.count(NodeID)>0) continue;	// already found
		dm[NodeID]=event.dist;
		
		AdjGrpAddr=getAdjListGrpAddr(NodeID);
		getFixedF(SIZE_A,Ref(AdjListSize),AdjGrpAddr);	// read # entries
		
		for (int z=0;z<AdjListSize;z++) {
			getVarE(ADJNODE_A	,Ref(NewNodeID)	,AdjGrpAddr,z);
			getVarE(DIST_A		,Ref(eKdist)	,AdjGrpAddr,z);
			
			StepEvent newevent=event;	// copy ...
			newevent.node=NewNodeID;
			newevent.dist=event.dist+eKdist;
			if (dm.count(NewNodeID)>0) continue;
			sQ.push(newevent);	// propagation
			
			getVarE(PTKEY_A	,Ref(PtGrpKey)	,AdjGrpAddr,z);
			if (PtGrpKey>=0) {	// valid PtGrpKey  (-1 for invalid key)
				float vJ;
				int PtGrpAddr,PtGrpSize,dummy;
				
				PtGrpAddr=pointQuery(PtTree,PtGrpKey,dummy);
				getFixedF(SIZE_P,Ref(PtGrpSize),PtGrpAddr);
	 			for (int j=0;j<PtGrpSize;j++) {	// scan
					getVarE(PT_P,Ref(vJ),PtGrpAddr,j);	// PtGrpKey+j
					if (vJ<eKdist/2.0)
						pcntm[NodeID]++;
					else 
						pcntm[NewNodeID]++;
					//printf("%d %d %f %f\n",NodeID,NewNodeID,vJ,event.dist+vJ);
				}
			}
		}
	}
	printf("Heap: max_size %d, pop_size %d\n",MaxHeapSize,AggPopSize);
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
	StepQueue sQ;
	StepEvent rootstep;
	
	rootstep.ClusID=0;
	rootstep.dist=0;	
	rootstep.node=0;	// Ni=0
	rootstep.isSortAcc=true;
	
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
	StepQueue sQ;
	StepEvent rootstep;
	
	rootstep.ClusID=0;
	rootstep.dist=0;	
	rootstep.node=0;	// Ni=0
	rootstep.isSortAcc=true;
	
	isSample=new bool[NodeNum];
	pcntm=new int[NodeNum];
	for (int i=0;i<NodeNum;i++) {
		pcntm[i]=0;
		isSample[i]=false;
	}
	
	while (sQ.size()>0) sQ.pop();	// clear sQ
	sQ.push(rootstep);
	CountExpansion(sQ);	// find nodes adjacent to populated edges
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
	for (int i=0;i<sampleCnt;i++) {
		//if (i!=66) continue;	// testing
		
		sample[i].cnt=-1;	// on-the-fly testing 28/6/2005 18:49
		
		/*rootstep.node=sample[i].NodeID;
		while (sQ.size()>0) sQ.pop();	// clear sQ
		sQ.push(rootstep);
		sample[i].cnt=RangeCount(sQ,anchor_rad,&(sample[i]));*/
		//printf("%d (%f %f): %d\n",i,sample[i].x,sample[i].y,sample[i].cnt);
		
		/*if (i>=10)
			exit(0);*/
	}
	
	printf("\nqueries\n");
	float sum_mnerr=0,sum_err=0,gSumErr=0,gSumAct=0,gMaxErr=0;
	for (int z=0;z<qCnt;z++) {
		float curx,cury;
		rootstep.node=nQ[z];
		int AdjGrpAddr=getAdjListGrpAddr(nQ[z]);
		curx=xcrd[nQ[z]];
		cury=ycrd[nQ[z]];
		
		while (sQ.size()>0) sQ.pop();	// clear sQ
		sQ.push(rootstep);	
		int curcnt=RangeCount(sQ,radius);
		int estcnt=0;
		
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
				rootstep.node=sample[id].NodeID;
				while (sQ.size()>0) sQ.pop();	// clear sQ
				sQ.push(rootstep);
				sample[id].cnt=RangeCount(sQ,anchor_rad,&(sample[id]));
			}
			
			//estcnt=sample[id].cnt;  // using nearest witness
		
			float a=sample[id].a;
			float b=sample[id].b;
			estcnt=(int)exp(a+b*log(radius));	
		}
		
		//estcnt=(int)(EuclidCnt(curx,cury,sample[id].eps));	// strange: get very good accuracy ???
		printf("%d: %d %d (%f)\n",z,curcnt,estcnt,dist);
		
		gSumErr+=ValAbs(curcnt-estcnt);
		gSumAct+=curcnt;
		
		sum_mnerr+=(float)min(curcnt,estcnt)/max(curcnt,estcnt);
		if (curcnt>0) {
			float cur_err=ValAbs(curcnt-estcnt)/(float)curcnt;
			sum_err+=cur_err;
			gMaxErr=max(cur_err,gMaxErr);
		}
	}
	
	printf("error: %0.4f (%0.4f/%0.4f)\n",gSumErr/gSumAct,gSumErr,gSumAct);	// error defined in [APR99]
	printf("quality: %f %f %f\n",sum_mnerr/qCnt,sum_err/qCnt,gMaxErr);
}


void ReadPtFiles() {
	// read point groups
	float pt_pos,edge_dist;
	int Ni,Nj,grp_sz;
	int pt_cnt=0;
	
	fseek(PtFile,0,SEEK_SET);
	while (!feof(PtFile)&&pt_cnt<PtNum) {
		fread(&(Ni),1,sizeof(int),PtFile);
		fread(&(Nj),1,sizeof(int),PtFile);
		fread(&(edge_dist),1,sizeof(float),PtFile);
		fread(&(grp_sz),1,sizeof(int),PtFile);
		for (int i=0;i<grp_sz;i++) {
			fread(&(pt_pos),1,sizeof(float),PtFile); 
			float xpt=xcrd[Ni]+(xcrd[Nj]-xcrd[Ni])*pt_pos/edge_dist;
			float ypt=ycrd[Ni]+(ycrd[Nj]-ycrd[Ni])*pt_pos/edge_dist;
			pt_xcrd.push_back(xpt);
			pt_ycrd.push_back(ypt);
		}
		pt_cnt+=grp_sz;
	}
	printf("check cnt: %d\n",pt_cnt);
}


void ReadCrd() {
	float xVal,yVal;
	
	xcrd.clear();
	ycrd.clear();
	for (int i=0;i<NodeNum;i++) {
		int AdjGrpAddr=getAdjListGrpAddr(i);
		getFixedF(XCRD_A,Ref(xVal),AdjGrpAddr);
		getFixedF(YCRD_A,Ref(yVal),AdjGrpAddr);
		xcrd.push_back(xVal);
		ycrd.push_back(yVal);
	}
}
// second est. method: using actual Euclidean result +  Euc/net. approx. ratio


//#define FLOATZERO       (1e-20)
#define FLOATZERO       (0)

/*****************************************************************
this function performs the minimum squared error processing
on the input 2-dim pts. It is effecitve for small amount of noise
*****IMPORTANT*****
ensure number in array y are sorted in ascending order
*******************
para:
x: x-coords of the pts
y: 
n: # of pts
cnst: the intercept of the line
corcoef: the correlation coefficient
[rst, red] specifies the range of the array in which linearity is significant
returns the slope
Coded by Yufei Tao dec 22 2002
*****************************************************************/

float MSE1(float *_x, float *_y, int _n, float &_cnst, float &_corcoef, float &_rst, float &_red) {
	float xsum=0, ysum=0, xxsum=0, yysum=0, xysum=0;
	float min_corcoef=0.00001, min_per=0.2, min_pts=5, itatime=0;
	float minx=1e20, maxx=-1, miny=1e20, maxy=-2;
	
	
	//get the range of x[] and y[]
	for (int i=0; i<_n; i++) {
		if (_x[i]<minx) minx=_x[i];
		if (_x[i]>maxx) maxx=_x[i];
		if (_y[i]<miny) miny=_y[i];
		if (_y[i]>maxy) maxy=_y[i];
	}
	float ylen=maxy-miny;  //range length of y[]

	//get rid of those points within the lower 20% 
	int st=0, ed; //starting and ending subscripts

	for (int i=0; i<_n; i++) {
		if (fabs(_y[i]-maxy)/ylen>min_per) {
			st=i; break;
		}
	}
				
	ed=st;
	
	for (int i=_n-1; i>ed; i--) {
		if (fabs(_y[i]-miny)/ylen>min_per) {
			ed=i; break;
		}
	}

	if (ed-st+1<=min_pts) { //too few points to do mse
		_cnst=0; _corcoef=1;
		_rst=0; _red=0; return 0;
	}
	
	//init the components of mse 
	int cnt=0;
	for (int i=st; i<=ed; i++) {
		xsum+=_x[i]; ysum+=_y[i];
		xxsum+=_x[i]*_x[i]; yysum+=_y[i]*_y[i];
		xysum+=_x[i]*_y[i];
		cnt++;
	}

	//perform the minimum squared errors 
	float ssxy=xysum-1.0f/cnt*xsum*ysum;
	float ssxx=xxsum-1.0f/cnt*xsum*xsum;
	assert(fabs(ssxx)>FLOATZERO);
	float ssyy=yysum-1.0f/cnt*ysum*ysum;
	
	float corcoef;
	if (fabs(ssyy)<FLOATZERO)
		corcoef=1;
	else
		corcoef=sqrt(ssxy*ssxy/ssxx/ssyy);

	float slope, cnst;
	slope=ssxy/ssxx; cnst=ysum/cnt-slope*xsum/cnt;
	_cnst=cnst; _corcoef=corcoef;
	_rst=st; _red=ed;  //if you validate the above lines, then comment this line								
	printf("%.5f, %.5f\n", corcoef, sqrt(ssxy*ssxy/ssxx/ssyy));
									
	return slope;
}


/*****************************************************************
calculates the fractal dimensionality of data within a region

 para:
N: the cardinality
U: axis length
d: dimensionality
mbr: the minimum bounding rectangle that constraints the set of points to calculate
fractal from
corFDim: the correlated fractal dimensionality
base: the log base (1.5-2 ?)
q: =0 box-count fractal, =2 correlated fractal
returns the fractal dimensionality
*****************************************************************/


/*float calPartFractal(int _N, float _U, int _d, float *_mbr, float _base) {
	float *pt=NULL;	// serial of data points from file
	//first check how many points lie in mbr 
	printf("filtering...");
	int cnt=0;
	//		cnt set to certain values ?
	printf("%d points left\n", cnt);

	//init an array with cnt size 
	printf("making a new array...");
	int cpycnt=0;
	float *ppt=new float[cnt*_d];
	for (int i=0; i<_N; i++)
	{
		//copy this point to ppt
		memcpy(&ppt[cpycnt*_d], &pt[i*_d], sizeof(float)*_d);

		//scale the new point to fit [0, _U] 
		for (int j=0; j<_d; j++)
			ppt[cpycnt*_d+j]=(ppt[cpycnt*_d+j]-_mbr[2*j])/(_mbr[2*j+1]-_mbr[2*j])*_U;
		
		cpycnt++;
	}
	printf("done. %d objs left\n", cpycnt);

	//now get the frequencies 
	float base=_base;
	int st=0, ed=-log(0.001)/log(base)+10;  //we automatically control the steps to compute
//	int st=0, ed=200;

	float *logs2=new float[ed-st+1], *logx=new float[ed-st+1];
	cnt=0;
	int lastm=0;
	for (int i=st; i<=ed; i++)
	{
		float r=pow(base, -i), s2;
		int m=1.0f/r;
		if (m!=lastm)
		{
			s2=qPowerFreq(ppt, cpycnt, _U, _d, m, 2);
			logs2[cnt]=log(s2)/log(base); logx[cnt]=log(1.0f/m)/log(base);
			printf("x=%.5f, y=%.5f\n", logx[cnt], logs2[cnt]);
			cnt++;
			lastm=m;
		}
	}

	float cnst, corcoef, rst, red;
	float locfd=MSE1(logx, logs2, cnt, cnst, corcoef, rst, red);
	printf("local fd=%.3f, cnst=%.3f\n", locfd, sqrt(pow(_base, cnst)));
	printf("corcoef=%.5f\n", corcoef);

	delete []logs2; delete []logx;
	
	delete []ppt;
	delete []pt;
	return locfd;
}*/

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
	ReadCrd();
	ReadPtFiles();
	
	//test2();
	
	int _smpCnt=getConfigInt("sc",cr);	//1000
	float _radius=getConfigFloat("e",cr);	//100
	float _nhood=getConfigFloat("nhood",cr);
	test(_smpCnt,_radius,_nhood);
	
	CloseDiskComm();
	PrintElapsed();
		
	printf("\n");
	return 0;
}
