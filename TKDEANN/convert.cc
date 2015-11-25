#include "utility.h"
#include "netshare.h"
#include <list>

int nextNode=0;
FastArray<float> xcrd,ycrd;

const float DOMAIN_SIZE=10000;

struct PointType{
	float x,y;
};

struct PointComp {
	bool operator () (PointType a,PointType b) const
    { if (a.x>b.x) return true;
      else if (a.x==b.x&&a.y>b.y) return true;
      else return false; }
};

map<PointType,int,PointComp> CountMap;

int getNode(float x,float y) {
	int result=nextNode;
	
	PointType key;
	key.x=x;	key.y=y;
	if (CountMap.count(key)==0) {
		CountMap[key]=nextNode;
		xcrd.push_back(x);	ycrd.push_back(y);
		nextNode++;
	} else 
		result=CountMap[key];
	
	return result;
}

float getDist(float x1,float y1,float x2,float y2) {
	return float(sqrt((x2-x1)*(x2-x1)+(y2-y1)*(y2-y1)));
}

float getDist(int Ni,int Nj) {
	float x1=xcrd[Ni],y1=ycrd[Ni],x2=xcrd[Nj],y2=ycrd[Nj];
	float dist= float(sqrt((x2-x1)*(x2-x1)+(y2-y1)*(y2-y1)));
	return dist;
}

void ConvertMaproomData(char* map_prefix) {
	char edgef[255],nodef[255],rawf[255];
	sprintf(rawf,"data/%s.raw",map_prefix);
	sprintf(nodef,"data/%s.cnode",map_prefix);
	sprintf(edgef,"data/%s.cedge",map_prefix);
	
	FILE* fp=fopen(rawf,"r");	// already normalized
	FILE* cedge=fopen(edgef,"w");
	FILE* cnode=fopen(nodef,"w");
	
	int count=0;
	float eDist,totalLen=0;
	float x1,y1,x2,y2;
	int id,Ni,Nj;
	while (!feof(fp)) {
		fscanf(fp, "%d %f %f %f %f\n", &id, &x1, &y1, &x2, &y2);
		eDist= float(sqrt((x2-x1)*(x2-x1)+(y2-y1)*(y2-y1)));
		totalLen+=eDist;
		Ni=getNode(x1,y1);	Nj=getNode(x2,y2);
		fprintf(cedge,"%d %d %d %f\n",id,min(Ni,Nj),max(Ni,Nj),eDist);
		count++;
	}
	
	for (int i=0;i<nextNode;i++)
		fprintf(cnode,"%d %f %f\n",i,xcrd[i],ycrd[i]);
	
	printf("%d edge, %d nodes, total length %0.3f\n",count,nextNode,totalLen);
	
	fclose(cedge);
	fclose(cnode);
	fclose(fp);
}

BitStore isVisited;		// for 	GetConnectedNodes use
FastArray<int> Partition,OldToNew,NewToOld;


void SpatialPartitionNodes() {
	Partition.assign(NodeNum,0);
	for (int i=0;i<NodeNum;i++) {
		// assume domain is [10000,10000]
		// gridsz: 4*4  or for 8*8
		int maxNormCrd=8;		//  4 or 8
		float gridsz=10000.0/maxNormCrd;
		
		int xn=(int)(xcrd[i]/gridsz);
		int yn=(int)(ycrd[i]/gridsz);
		xn=min(xn,maxNormCrd-1);
		yn=min(yn,maxNormCrd-1);
		Partition[i]=xn*maxNormCrd+yn;
	}
}

int GetConnectedNodes() {
	isVisited.assign(NodeNum,false);
	
	list<int> seeds;
	//seeds.push_back( NodeNum/2 );	// push node 0
	seeds.push_back(0);
	int cnt=0;
	
	while (seeds.size()>0) { // seeds <> Empty
		int node=seeds.front();
		seeds.pop_front();
		
		if (isVisited[node]) continue;
		isVisited[node]=true;
		cnt++;
		
		for (int k=0;k<AdjList[node].size();k++) {
			int Nk=AdjList[node][k];
			if (!isVisited[Nk]) seeds.push_back(Nk);
		}
	}
	
	OldToNew.assign(NodeNum,-1);
	for (int i=0;i<NodeNum;i++) {
		if (isVisited[i]) {
			OldToNew[i]=NewToOld.size();
			NewToOld.push_back(i);
		}
	}
	return cnt;
}

void NormalizeCrd() {
	// perform normalization
	float minx,miny,maxx,maxy;
	minx=maxx=xcrd[0];	miny=maxy=ycrd[0];
	
	for (int i=0;i<NodeNum;i++) {
		if (xcrd[i]<minx) minx=xcrd[i];
		if (maxx<xcrd[i]) maxx=xcrd[i];
		if (ycrd[i]<miny) miny=ycrd[i];
		if (maxy<ycrd[i]) maxy=ycrd[i];
	}
	
	for (int i=0;i<NodeNum;i++) {
		xcrd[i]=(xcrd[i]-minx)/(maxx-minx)*DOMAIN_SIZE;
		ycrd[i]=(ycrd[i]-miny)/(maxy-miny)*DOMAIN_SIZE;
	}
}

void ReadNodeFile(const char* nodef,int& NodeCnt) {
	int id;
	float x,y;
	
	FILE* cnode=fopen(nodef,"r");
	CheckFile(cnode,nodef);
	NodeCnt=0;	// getKey depends on NodeNum so we have to init. it first
	
	while (!feof(cnode)) {
		fscanf(cnode,"%d %f %f\n", &id, &x, &y);
		xcrd.push_back(x);
		ycrd.push_back(y);
		NodeCnt++;
	}
	printf("Nodes read, ");		PrintElapsed();
	fclose(cnode);
}

void ReadEdgeFile(const char* edgef) {
	int id,Ni,Nj;
	float dist;
	
	FILE* cedge=fopen(edgef,"r");
	CheckFile(cedge,edgef);	
	
	AdjList=new FastArray<int>[NodeNum];
	while (!feof(cedge)) {
		fscanf(cedge, "%d %d %d %f\n", &id, &Ni, &Nj, &dist);
		AdjList[Ni].push_back(Nj);	AdjList[Nj].push_back(Ni);
	}
	printf("Edges read, ");		PrintElapsed();
	fclose(cedge);
}

void MakeConnectedGraph(char* map_prefix) {
	float dist,x1,y1,x2,y2;	
	char edgef[255],nodef[255],newedgef[255],newnodef[255];
		
	sprintf(edgef,"data/%s.cedge",map_prefix);
	sprintf(nodef,"data/%s.cnode",map_prefix);
	sprintf(newedgef,"data/%s.newedge",map_prefix);
	sprintf(newnodef,"data/%s.newnode",map_prefix);
	
	FILE* newedge=fopen(newedgef,"w");
	FILE* newnode=fopen(newnodef,"w");
	
	ReadNodeFile(nodef,NodeNum);
	ReadEdgeFile(edgef);
	int NewNodeNum=GetConnectedNodes();
	printf("%d connected nodes\n",NewNodeNum);
	
	for (int i=0;i<NewNodeNum;i++) {
		int OldIndex=NewToOld[i];
		xcrd[i]=xcrd[OldIndex];
		ycrd[i]=ycrd[OldIndex];
		AdjList[i]=AdjList[OldIndex];
		for (int k=0;k<AdjList[i].size();k++) {
			int OldValue=AdjList[i][k];
			AdjList[i][k]=OldToNew[OldValue];
		}
	}
	
	CheckFile(newedge,newedgef);	CheckFile(newnode,newnodef);
	int EdgeID=0; float totalLen=0;
	for (int i=0;i<NewNodeNum;i++) {
		fprintf(newnode,"%d %f %f\n",i,xcrd[i],ycrd[i]);
		x1=xcrd[i]; y1=ycrd[i];
		for (int k=0;k<AdjList[i].size();k++) {
			int AdjNode=AdjList[i][k];
			if (AdjNode>i) {
				x2=xcrd[AdjNode]; y2=ycrd[AdjNode];
				dist= float(sqrt((x2-x1)*(x2-x1)+(y2-y1)*(y2-y1)));
				totalLen+=dist;
				fprintf(newedge,"%d %d %d %f\n",EdgeID,i,AdjNode,dist);
				EdgeID++;
			}
		}
	}
	fclose(newedge);	fclose(newnode);
	printf("%d edge, %d nodes, total length %0.3f\n",EdgeID,NewNodeNum,totalLen);
}

void EdgeLengthCheck(char* map_prefix) {
	char edgef[255],nodef[255],newedgef[255];
		
	sprintf(edgef,"data/%s.cedge",map_prefix);
	sprintf(nodef,"data/%s.cnode",map_prefix);	
	sprintf(newedgef,"data/%s.newedge",map_prefix);
	
	FILE* cedge=fopen(edgef,"r");
	FILE* newedge=fopen(newedgef,"w");
	
	CheckFile(cedge,edgef);
	CheckFile(newedge,newedgef);
	
	int id,Ni,Nj;
	float dist,x1,y1,x2,y2;	
	
	ReadNodeFile(nodef,NodeNum);
	
	float totalLen=0;
	while (!feof(cedge)) {
		fscanf(cedge, "%d %d %d %f\n", &id, &Ni, &Nj, &dist);
		x1=xcrd[Ni]; y1=ycrd[Ni];
		x2=xcrd[Nj]; y2=ycrd[Nj];
		dist= float(sqrt((x2-x1)*(x2-x1)+(y2-y1)*(y2-y1)));
		totalLen+=dist;
		fprintf(newedge,"%d %d %d %f\n",id,Ni,Nj,dist);		
	}
	printf("Edges read, ");		PrintElapsed();
	
	printf("total length %0.3f\n",totalLen);
	
	fclose(cedge);
	fclose(newedge);
}

void FindNNnode(char* map_prefix,float px,float py) {
	char edgef[255],nodef[255];
	int id,NNnode=-1;
	float dist,x,y,NNdist=1.0e20,bx,by;	
	
	sprintf(nodef,"data/%s.cnode",map_prefix);
	FILE* cnode=fopen(nodef,"r");
	CheckFile(cnode,nodef);		
	while (!feof(cnode)) {
		fscanf(cnode,"%d %f %f\n", &id, &x, &y);
		dist=float(sqrt((px-x)*(px-x)+(py-y)*(py-y)));
		if (dist<NNdist) {
			NNdist=dist;	NNnode=id;
			bx=x;	by=y;
		}
	}
	printf("Nodes read, ");		PrintElapsed();	
	fclose(cnode);	
	printf("node %d (%f,%f) with dist %f\n",NNnode,bx,by,NNdist);
}

// assume unvisited firstId
int GroupConnectedNodes(int group_size,int firstId,FastArray<int>& NewToOldOrder,int cur_part) {
	StepQueue sQ;
	StepEvent stepZ;
	stepZ.node=firstId;	stepZ.dist=0;
	sQ.push(stepZ);	
	
	int cnt=0,for_cnt=0;	
	while (!sQ.empty()) {
		StepEvent event=sQ.top();
		sQ.pop();
		
		int node=event.node;
		if (isVisited[node]) continue;
		
		isVisited[node]=true;	NewToOldOrder.push_back(node);
		cnt++;
		if (cnt>=group_size) break;
		
		for (int k=0;k<AdjList[node].size();k++) {
			int Nk=AdjList[node][k];
			if (!isVisited[Nk]&&Partition[Nk]==cur_part) {
				StepEvent newevent=event;	// copy ...
				newevent.node=Nk;	// should be better !
				newevent.dist=event.dist+getDist(node,Nk);	// dij_FS
				
				//for_cnt++;	newevent.dist=for_cnt;	// BFS
				sQ.push(newevent);
			}
		}
	}	
	return cnt;
}

void GroupNodesByPage(char* map_prefix,int group_size) {
	int id,Ni,Nj;
	float dist,x1,y1,x2,y2;	
	char edgef[255],nodef[255];
	
	sprintf(edgef,"data/%s.cedge",map_prefix);
	sprintf(nodef,"data/%s.cnode",map_prefix);
	
	ReadNodeFile(nodef,NodeNum);
	//NormalizeCrd();		// only for converting a new map ???
	
	FILE* cedge=fopen(edgef,"r");
	CheckFile(cedge,edgef);
	AdjList=new FastArray<int>[NodeNum];
	while (!feof(cedge)) {
		fscanf(cedge, "%d %d %d %f\n", &id, &Ni, &Nj, &dist);
		AdjList[Ni].push_back(Nj);	AdjList[Nj].push_back(Ni);
	}
	printf("Edges read, ");		PrintElapsed();
	
	bool isSpatPart=false;
	
	if (isSpatPart)
		SpatialPartitionNodes();	//	partition the nodes using spatial coordinates first
	else
		Partition.assign(NodeNum,0);
	
	int GroupNum;
	isVisited.assign(NodeNum,false);
	FastArray<int> NewToOldOrder,OldToNewOrder;
	
	// grouping using DFS or dijstra
	for (int i=0;i<NodeNum;i++)
		if (!isVisited[i]) {
			GroupNum=GroupConnectedNodes(group_size,i,NewToOldOrder,Partition[i]);
			//printf("%d grouped nodes\n",GroupNum);
		}
	
	// control exp. test: random grouping (expected to have the worst result)
	/*for (int i=0;i<NodeNum;i++) NewToOldOrder.push_back(i);
	for (int i=0;i<NodeNum;i++) {
		int j=rand()%NodeNum;
		int tmp=NewToOldOrder[i];
		NewToOldOrder[i]=NewToOldOrder[j];
		NewToOldOrder[j]=tmp;
	}*/
	printf("Grouping finished, ");		PrintElapsed();
	
	FILE* newedge=fopen("data/grouped.cedge","w");
	FILE* newnode=fopen("data/grouped.cnode","w");
	
	OldToNewOrder.assign(NodeNum,-1);
	for (int i=0;i<NodeNum;i++) {
		int Nold=NewToOldOrder[i];
		OldToNewOrder[Nold]=i;
		fprintf(newnode,"%d %f %f\n",i,xcrd[Nold],ycrd[Nold]);
	}
	printf("Nodes written, ");		PrintElapsed();
	
	fseek(cedge,0,SEEK_SET);	// set cursor to the beginning again !
	while (!feof(cedge)) {
		fscanf(cedge, "%d %d %d %f\n", &id, &Ni, &Nj, &dist);
		//dist=getDist(xcrd[Ni],ycrd[Ni],xcrd[Nj],ycrd[Nj]);	// ???
		int Ri=OldToNewOrder[Ni],Rj=OldToNewOrder[Nj];
		Ni=min(Ri,Rj);	Nj=max(Ri,Rj);
		fprintf(newedge, "%d %d %d %f\n", id, Ni, Nj, dist);
	}	
	printf("Edges written, ");		PrintElapsed();
	fclose(cedge);
	fclose(newedge);	fclose(newnode);
}

struct ValuePair {
	float value;
	int id;
};

struct ValueComp {
	bool operator () (ValuePair left,ValuePair right) const
    { return left.value > right.value; }
};

struct IdComp {
	bool operator () (ValuePair left,ValuePair right) const
    { return left.id < right.id; }
};


FastArray<int> DomCenter;
FastArray<float> CenterDist;
FastArray<ValuePair> TmpPairs;

void FindBNodeDist(int cur_part,int cur_pos) {
	StepQueue sQ;
	StepEvent stepZ;
	stepZ.node=DomCenter[cur_pos];		stepZ.dist=0;
	sQ.push(stepZ);
	
	TmpPairs.clear();
	while (!sQ.empty()) {
		StepEvent event=sQ.top();
		sQ.pop();
		
		int node=event.node;
		if (isVisited[node]) continue;
		isVisited[node]=true;	// bitmap may be too expensive
		if (Partition[node]!=cur_part) continue;
		
		// safest
		for (int i=0;i<DomCenter.size();i++) {
			if (DomCenter[i]==node) {
				ValuePair tmp;
				tmp.value=event.dist;
				tmp.id=i;
				TmpPairs.push_back(tmp);
			}
		}
		
		for (int k=0;k<AdjList[node].size();k++) {
			int Nk=AdjList[node][k];
			if (!isVisited[Nk]) {
				StepEvent newevent=event;	// copy ...
				newevent.node=Nk;
				newevent.dist=event.dist+getDist(node,Nk);
				sQ.push(newevent);
			}
		}
	}
}

void WriteHiTiGraph(FILE* cmat,int* submapcnt,int maxCrdSq) {
	// count number of boundary nodes
	FastArray<int> CurB_Nodes[maxCrdSq];
	int bnodecnt[maxCrdSq];
	fill(bnodecnt,bnodecnt+maxCrdSq,0);	
	
	fwrite(&NodeNum,1,sizeof(int),cmat);
	fwrite(&maxCrdSq,1,sizeof(int),cmat);
	fwrite(&(Partition[0]),NodeNum,sizeof(int),cmat);
		
	for (int i=0;i<NodeNum;i++) {
		int my_part=Partition[i];
		FastArray<int>& CurAdjList=AdjList[i];
		int crash=0;
		for (int z=0;z<CurAdjList.size();z++) {
			int j=CurAdjList[z];
			if (Partition[j]!=my_part) crash++;
		}
		if (crash>0) {
			bnodecnt[my_part]++;
			CurB_Nodes[my_part].push_back(i);
		}
	}
	
	int totalspace=0;
	for (int part=0;part<maxCrdSq;part++) {
		sort(CurB_Nodes[part].begin(),CurB_Nodes[part].end());
	
		int bn=CurB_Nodes[part].size();
		fwrite(&bn,1,sizeof(int),cmat);		// # of border nodes
		fwrite(&(CurB_Nodes[part][0]),bn,sizeof(int),cmat);	// border node IDs		
				
		// reuse DomCenter now !
		DomCenter.assign(CurB_Nodes[part].begin(),CurB_Nodes[part].end());
		for (int cur_pos=1;cur_pos<DomCenter.size();cur_pos++) {
			isVisited.assign(NodeNum,false);
			FindBNodeDist(part,cur_pos);
			
			// write within links
			sort(TmpPairs.begin(),TmpPairs.end(),IdComp());
			//printf("[%d (%d)] ",cur_pos,TmpPairs.size());
			for (int t=0;t<cur_pos;t++) {
				fwrite(&(TmpPairs[t].value),1,sizeof(float),cmat);
				//printf("%d: %0.3f ",TmpPairs[t].id,TmpPairs[t].value);
			}
			//printf("\n");
		}
		printf("count(%d)=%d, bn=%d\n",part,submapcnt[part],bnodecnt[part]);
		totalspace+=bnodecnt[part]*(bnodecnt[part]-1)/2;
	}
	printf("space: %d\n",totalspace);
}

void ConcurrentExpansion(int maxCrdSq,int* submapcnt) {
	StepQueue sQ,tmpQ;
	for (int part=0;part<maxCrdSq;part++) {
		StepEvent stepZ;
		stepZ.node=DomCenter[part];
		stepZ.dist=0;
		stepZ.ClusID=part;
		sQ.push(stepZ);	
	}
	
	int cnt_limit=(int)(0.5*NodeNum/maxCrdSq);
	cnt_limit=INT_MAX;	// no effect
	
	int NumNodesVisited=0;
	isVisited.assign(NodeNum,false);
	
	for (int cases=0;cases<2;cases++) {
		while (!sQ.empty()) {
			StepEvent event=sQ.top();
			sQ.pop();
			
			int node=event.node;
			if (isVisited[node]) continue;
			
			if (cases==0&&submapcnt[event.ClusID]>=cnt_limit) {
				tmpQ.push(event);
				continue;
			}
			
			isVisited[node]=true;
			NumNodesVisited++;
			
			Partition[node]=event.ClusID;
			submapcnt[event.ClusID]++;
			
			for (int k=0;k<AdjList[node].size();k++) {
				int Nk=AdjList[node][k];
				if (!isVisited[Nk]) {
					StepEvent newevent=event;	// copy ...
					newevent.node=Nk;	// should be better !
					newevent.dist=event.dist+getDist(node,Nk);	// dij_FS
					sQ.push(newevent);
				}
			}
		}
		if (cases==0) {
			while (!tmpQ.empty()) {
				sQ.push(tmpQ.top());
				tmpQ.pop();
			}
		}
	}
	printf("%d \n",NumNodesVisited);
//	exit(0);
}

void BuildMatGraph(char* map_prefix) {
	int id,Ni,Nj;
	float dist,x1,y1,x2,y2;	
	char edgef[255],nodef[255],matf[255];
	
	sprintf(edgef,"data/%s.cedge",map_prefix);
	sprintf(nodef,"data/%s.cnode",map_prefix);
	
	ReadNodeFile(nodef,NodeNum);
	ReadEdgeFile(edgef);
	
	// assume domain is [10000,10000]
	// STR partition, sort in groups ordered by x-crd, then break a
	// group into several groups ordered by y-crd
	srand(0);	// avoid random effect
	int maxNormCrd=5;	//  N^(1/4)   8
	int maxCrdSq=maxNormCrd*maxNormCrd;
	FastArray<ValuePair> vAry;
	FastArray<ValuePair> wAry[maxNormCrd];
	int submapcnt[maxCrdSq];
	fill(submapcnt,submapcnt+maxCrdSq,0);
	
	Partition.assign(NodeNum,-1);
	DomCenter.assign(maxCrdSq,-1);
	// try both BFS and DFS
	
	for (int i=0;i<NodeNum;i++) {
		ValuePair tmp;
		tmp.id=i;
		tmp.value=xcrd[tmp.id];
		vAry.push_back(tmp);
	}
	sort(vAry.begin(),vAry.end(),ValueComp());
	for (int i=0;i<vAry.size();i++) {
		int xn=(i*maxNormCrd)/vAry.size();
		ValuePair tmp;
		tmp.id=vAry[i].id;
		tmp.value=ycrd[tmp.id];		
		wAry[xn].push_back(tmp);
	}
	vAry.clear();
	
	for (int xn=0;xn<maxNormCrd;xn++) {
		FastArray<ValuePair>& tAry=wAry[xn];
		sort(tAry.begin(),tAry.end(),ValueComp());
		
		// nearest to centroid
		for (int yn=0;yn<maxNormCrd;yn++) {
			float diff=0.5;		// default
			//diff=drand48();
			int pos=(int)((yn+diff)/maxNormCrd*tAry.size());
			// if error, then check here !
			DomCenter[xn*maxNormCrd+yn]=tAry[pos].id;
		}

		tAry.clear();
	}
	
	// must ensure that all subgraph are connected by concurrent expansion
	ConcurrentExpansion(maxCrdSq,submapcnt);
	int tsubcnt=0;
	printf("Checking ...\n");
	for (int i=0;i<maxCrdSq;i++) {
		printf("%d: %d\n",i,submapcnt[i]);
		tsubcnt+=submapcnt[i];
	}
	printf("tsubcnt: %d\n",tsubcnt);
	
	// write data
	sprintf(matf,"data/%s.cmat",map_prefix);
	FILE* cmat=fopen(matf,"w");
	CheckFile(cmat,matf);
	
	WriteHiTiGraph(cmat,submapcnt,maxCrdSq);	// ICDE'96 (TKDE'02)
	
	fclose(cmat);
}

int main(int argc,char** argv) {
	InitClock();
	
	char* prefix="cn_rr";
	ConvertMaproomData(prefix);
	
	// "SF","TG","NA","OL"
	//char* prefix="TG";
	//char* prefix="OL_dij";
	//char* prefix="SF_dnw";
	//char* prefix="cnrail";
	
	MakeConnectedGraph(prefix);
	//EdgeLengthCheck(prefix);
	
	// [NN-search application]
	//if (argc>=3) FindNNnode(prefix,atof(argv[1]),atof(argv[2]));
	
	// [Page ordering of the nodes for locality]
	// page size: 4K ; Avg fan out: 2.549 
	// have to be changed later as structure changes !!!
//	int ADJGRP_HEADSIZE=sizeof(int)+2*sizeof(float);
//	int ADJGRP_ITEMSIZE=2*sizeof(int)+sizeof(float);
//	int group_size=(int)(4096/(ADJGRP_HEADSIZE+2.549*ADJGRP_ITEMSIZE));
//	GroupNodesByPage(prefix,group_size);
	
	// [Build materialized graph]
	//BuildMatGraph(prefix);
	
	PrintElapsed();
	return 0;
}
