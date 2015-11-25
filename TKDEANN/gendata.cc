#include "btree.h"
#include "utility.h"
#include "netshare.h"
#include "rtree.h"
#include "rtacc.h"


int num_D;

float getDist(float x1,float y1,float x2,float y2) {
	return float(sqrt((x2-x1)*(x2-x1)+(y2-y1)*(y2-y1)));
}


#define MAXLEVEL 10

// to be initialized
char **cur_block;
int *cur_block_offset, *cur_node_maxkey, *cur_node_entries;
char *block;
int i_capacity,root_block_id,num_written_blocks,top_level;
int BlkLen;

int PtMaxKey=0;
bool PtMaxUsing=false;	// for solving the bug that was to find
bool IS_NODE_SIZE_GIVEN=false;
FastArray<float> xcrd,ycrd;

BTree* BT_initialize(char *treename) {
	BTree *bt;
	
	cur_block = new char*[MAXLEVEL]; //ptr to current node at each level
	cur_block_offset = new int[MAXLEVEL];
	cur_node_maxkey = new int[MAXLEVEL];
	cur_node_entries = new int[MAXLEVEL];
	top_level = 0;
		
	for (int i=0;i<MAXLEVEL; i++) cur_block[i] = NULL;
	block = new char[BlkLen];
	
	// number of cached file entries: 128
	bt = new BTree(treename, BlkLen, 128);
	i_capacity = (BlkLen - sizeof(char) - sizeof(int))/(sizeof(int)+sizeof(int));
	printf("i_capacity=%d\n", i_capacity);
	root_block_id = bt->root_ptr->block;
	num_written_blocks = 0;
	return  bt;
}

void BT_addentry(BTree *bt,int *top_level,int capacity,int level,int key,int *block_id,int RawAddr=0) {
	if (cur_block[level] == NULL) { //new node to be created
		if ((*top_level) < level) //new root
			*top_level = level;
		cur_block[level] = new char[BlkLen];
		
		char l = (char)level;
		memcpy(cur_block[level], &l, sizeof(char));
		
		cur_block_offset[level]=sizeof(char)+sizeof(int);
		cur_node_entries[level] = 0;
	}
	cur_node_maxkey[level]= key;
	if ((level==1)&&PtMaxUsing) cur_node_maxkey[level]=max(PtMaxKey,key);	// new change !!!	
	
	//copy key as new current entry and also the pointer to lower node
	memcpy(cur_block[level]+cur_block_offset[level], &key, sizeof(int));
	cur_block_offset[level]+=sizeof(int);
	
	//********* (Xblock_id for raw sequential page !)
	int Xblock_id=(level==1)?(RawAddr):(*block_id);
	memcpy(cur_block[level]+cur_block_offset[level], &Xblock_id, sizeof(int));
	
	cur_block_offset[level]+=sizeof(int);
	cur_node_entries[level]++;
	
	if (cur_node_entries[level] == capacity) { //node is full
		//copy capacity information
		memcpy(cur_block[level]+sizeof(char), &capacity, sizeof(int));
		//add maxkey of this node to upper level node
		bt->file->append_block(cur_block[level]);
		(*block_id)++;
		bt->num_of_inodes++;
		BT_addentry(bt, top_level, capacity, level+1, 
				cur_node_maxkey[level], block_id);
		delete [] cur_block[level];
		cur_block[level] = NULL;
	}
}

void BT_finalize(BTree* bt) {
	//flush non-empty blocks
	for (int level=1; level<= top_level; level++) {
		if (cur_block[level] != NULL) {
			//copy capacity information
			memcpy(cur_block[level]+sizeof(char), &cur_node_entries[level], sizeof(int));
			//add mbr of this node to upper level node
			if (level==top_level) {
				//root
				bt->file->write_block(cur_block[level], root_block_id);
				bt->num_of_inodes++;
				bt->root_ptr = NULL;
				bt->load_root();
				//printf("root written, id=%d\n", root_block_id);
			} else {
				bt->file->append_block(cur_block[level]);
				num_written_blocks++;
				bt->num_of_inodes++;
				BT_addentry(bt, &top_level, i_capacity, level+1, 
				cur_node_maxkey[level], &num_written_blocks);
			}
			delete [] cur_block[level];
			cur_block[level] = NULL;
		}
	}		
	delete [] block;
}

// Pt FlatFile Field:
//		Header:	Ni(int), Nj(int), dist(float), size(int)
//		Entry:	vP(float)
void makePtFiles(FILE *ptFile,char* treefile) {
	PtMaxUsing=true;
	BTree* bt=BT_initialize(treefile);		
	printf("making PtFiles\n");
	
	int RawAddr=0,key=0,size;	// changed
 	EdgeMapType::iterator iter=EdgeMap.begin();
 	while (iter!=EdgeMap.end()) {
 		edge* e=iter->second;
 		if (e->pts.size()>0) {	// do not index empty groups
			sort(e->pts.begin(),e->pts.end());
			
			RawAddr=ftell(ptFile);	// set addr to correct amt.
			size=e->pts.size();
			fwrite(&(e->Ni),1,sizeof(int),ptFile);
			fwrite(&(e->Nj),1,sizeof(int),ptFile);
			fwrite(&(e->dist),1,sizeof(float),ptFile);
			fwrite(&(size),1,sizeof(int),ptFile);
			fwrite(&(e->pts[0]),e->pts.size(),sizeof(float),ptFile);
			e->FirstRow=key;
			PtMaxKey=key+e->pts.size()-1;	// useful for our special ordering !
			
			//printf("(key,value)=(%d,%d)\n",key,RawAddr);
			BT_addentry(bt,&top_level,i_capacity,1,key,&num_written_blocks,RawAddr);
			key+=e->pts.size(); // changed, moved to after BT_addentry
		} else 
			e->FirstRow=-1;		// also later used by AdjFile
		iter++;
	}
	BT_finalize(bt);
	bt->UserField=num_D;
	delete bt;
	PtMaxUsing=false;
}



// change here when ...
#define MAXRECT  500000


int sort_center_y(const void *d1, const void *d2) {
    DATA *s1, *s2;
    float l1, u1, l2, u2;
    s1 = *((DATA **) d1);
    s2 = *((DATA **) d2);
    l1=s1->data[2];	u1=s1->data[3];
    l2=s2->data[2];	u2=s2->data[3];
    if (l1+u1<l2+u2)
        return -1;
    if (l1+u1>l2+u2)
        return 1;
    return 0;
}

int sort_center_x(const void *d1, const void *d2) {
    DATA *s1, *s2;
    float l1, u1, l2, u2;
    s1 = *((DATA **) d1);
    s2 = *((DATA **) d2);
    l1=s1->data[0];	u1=s1->data[1];
    l2=s2->data[0];	u2=s2->data[1];
    if (l1+u1<l2+u2)
        return -1;
    if (l1+u1>l2+u2)
        return 1;
    return 0;
}

void RT_addrect(RTree *rt, int *top_level, int blocklength, int capacity, int level, float *mbr, int *block_id, char **cur_block, int *cur_block_offset, float **cur_node_mbr, int *cur_node_mbrs) {
    //printf("RT_addrect: mbr=%f %f %f %f, blockid = %d\n", mbr[0], mbr[1], mbr[2], mbr[3], *block_id);
    if (cur_block[level] == NULL) { //new node to be created
        if ((*top_level) < level) //new root
            *top_level = level;
        cur_block[level] = new char[blocklength];
        
        //initialize node mbr
        for (int i=0; i<rt->dimension; i++) {
            cur_node_mbr[level][2*i] = MAXREAL;
            cur_node_mbr[level][2*i+1] = -MAXREAL;
        }
        bool son_is_data = (level == 1);
        memcpy(cur_block[level], &son_is_data, sizeof(bool));
        char l = (char)level;
        memcpy(cur_block[level]+sizeof(bool), &l, sizeof(char));
        cur_block_offset[level]=sizeof(bool)+sizeof(char)+sizeof(int);
        cur_node_mbrs[level] = 0;
    }

    for (int i=0; i<rt->dimension; i++) {
        if (mbr[2*i] < cur_node_mbr[level][2*i])
            cur_node_mbr[level][2*i] = mbr[2*i];
        if (mbr[2*i+1] > cur_node_mbr[level][2*i+1])
            cur_node_mbr[level][2*i+1] = mbr[2*i+1];
    }
    
    memcpy(cur_block[level]+cur_block_offset[level], mbr, 2*rt->dimension*sizeof(float));
    cur_block_offset[level]+=2*rt->dimension*sizeof(float);
    memcpy(cur_block[level]+cur_block_offset[level], block_id, sizeof(int));
    cur_block_offset[level]+=sizeof(int);
    cur_node_mbrs[level]++;
    //printf("level=%d, cur_node_mbrs[level]=%d\n", level, cur_node_mbrs[level]);

    if (cur_node_mbrs[level] == capacity) { //node is full
        //copy capacity information
        memcpy(cur_block[level]+sizeof(bool)+sizeof(char), &capacity, sizeof(int));
        //add mbr of this node to upper level node
        rt->file->append_block(cur_block[level]);
        //printf("intermediate block no %d written\n", *block_id);
        (*block_id)++;
        rt->num_of_inodes++;
        RT_addrect(rt, top_level, blocklength, capacity, level+1, cur_node_mbr[level], block_id, cur_block, cur_block_offset, cur_node_mbr, cur_node_mbrs);
        delete [] cur_block[level];
        cur_block[level] = NULL;
    }
}

RTree *RT_bulkload(char *treename, DATA **rects, long count, int blocklength) {
    int i,j,size,tmp;
    int offset;
    RTree *rt;
    float *node_mbr = new float[2*DIMENSION];
    int top_level = 0;
    int root_block_id;
	char level;
	
    char **cur_block = new char*[5]; //ptr to current node at each level
    int *cur_block_offset = new int[5];
    float **cur_node_mbr = new float*[5];
    int *cur_node_mbrs = new int[5];
    for (i=0; i<5; i++) {
        cur_block[i] = NULL;
        cur_node_mbr[i] = new float[2*DIMENSION];
    }
    char *block = new char[blocklength];
    
    rt = new RTree(treename, blocklength, 128, DIMENSION);
    int d_capacity = rt->root_ptr->capacity;
    int i_capacity = (blocklength - sizeof(bool) - sizeof(char) - sizeof(int))/(sizeof(int)+2*DIMENSION*sizeof(float));
    printf("d_capacity=%d, i_capacity=%d\n", d_capacity, i_capacity);
    
    // [LEL97, STR packing] sort the rects
    qsort(&rects[0],count,sizeof(DATA*),sort_center_x);
    int num_of_leaves = count/d_capacity + (count%d_capacity == 0 ? 0 : 1);
    int num_of_stripes = (int)sqrt((double)num_of_leaves);
    int count_in_stripe = num_of_leaves/num_of_stripes;
    printf("# leaves: %d, # stripes: %d, count/stripe: %d\n",num_of_leaves,num_of_stripes,count_in_stripe);
    
    // sort segments of rects w.r.t. center point on y-axis
    for (i=0; i<num_of_stripes; i++) {
    	int local_numitem=min(count_in_stripe*d_capacity,count-i*count_in_stripe*d_capacity);
        //printf("%d: %d\n", i*count_in_stripe*d_capacity, local_numitem);
        qsort(&rects[i*count_in_stripe*d_capacity],local_numitem,sizeof(DATA*),sort_center_y);
    }
    if (num_of_leaves%num_of_stripes) { //more rectangles
        int local_numitem=min(count_in_stripe*d_capacity,count-i*count_in_stripe*d_capacity);
        //printf("%d: %d\n", i*count_in_stripe*d_capacity, local_numitem);
        qsort(&rects[i*count_in_stripe*d_capacity],local_numitem,sizeof(DATA*), sort_center_y);
    }
	printf("sorted %d rect.\n", count);    
	//for (int i=0; i<count; i++)	rects[i]->print();	
	
    root_block_id = rt->root_ptr->block;
    rt->num_of_data = count;
    rt->num_of_dnodes = 0;
    int num_written_blocks = 0;

    int cur_rect = 0; //current position in rects array
    while(cur_rect<count) {
        for (i=0; i<DIMENSION; i++) {
            node_mbr[2*i] = MAXREAL;
            node_mbr[2*i+1] = -MAXREAL;
        }
        size = min(d_capacity, count-cur_rect);
        //printf("size=%d, capacity=%d\n",size,d_capacity);

        level = 0;
        memcpy(block, &level, sizeof(char));
        memcpy(block+sizeof(char), &size, sizeof(int));
        offset = sizeof(char) + sizeof(int);
        for (i=0; i<size; i++) {
            rects[cur_rect]->write_to_buffer(block + offset);
            // need to change ! (for rects)
            for (j=0; j<DIMENSION; j++) {
				if (rects[cur_rect]->data[2*j] < node_mbr[2*j])
                    node_mbr[2*j] = rects[cur_rect]->data[2*j];
                if (rects[cur_rect]->data[2*j+1] > node_mbr[2*j+1])
                    node_mbr[2*j+1] = rects[cur_rect]->data[2*j+1];
            }
            // remember to change here when size of data changes !!! (for rect.)
            offset+= 2*DIMENSION*sizeof(float)+2*sizeof(int); // check write_to_buffer for size !
            cur_rect++;
        }
        j = rt->file->append_block(block);
        //printf("written block %d, num_written_blocks=%d\n", j , num_written_blocks);
        rt->num_of_dnodes++;
        num_written_blocks++;
        //printf("block no %d written\n", num_written_blocks);
        RT_addrect(rt, &top_level, blocklength, i_capacity, 1, node_mbr, &num_written_blocks, cur_block, cur_block_offset, cur_node_mbr, cur_node_mbrs);
        //printf("node_mbr=%f %f %f %f\n", node_mbr[0], node_mbr[1], node_mbr[2], node_mbr[3]);
    }

    //flush non-empty blocks
    for (level=1; level<= top_level; level++) {
        if (cur_block[level] != NULL) {
            //copy capacity information
            memcpy(cur_block[level]+sizeof(bool)+sizeof(char), &cur_node_mbrs[level], sizeof(int));
            //add mbr of this node to upper level node
            if (level == top_level) {
                //root
                rt->file->write_block(cur_block[level], root_block_id);
                rt->num_of_inodes++;
                rt->root_is_data = FALSE;
                rt->root_ptr = NULL;
                rt->load_root();
                //rt->root_ptr->read_from_buffer(cur_block[level]);
                //printf("root written, id=%d\n", root_block_id);
            } else {
                rt->file->append_block(cur_block[level]);
                num_written_blocks++;
                rt->num_of_inodes++;
                RT_addrect(rt, &top_level, blocklength, i_capacity, level+1, cur_node_mbr[level], &num_written_blocks, cur_block, cur_block_offset, cur_node_mbr, cur_node_mbrs);
            }
            delete [] cur_block[level];
            cur_block[level] = NULL;
        }
    }
    delete [] block;
    return rt;
}

void makePtGrpRTree(char* rtfile) {
	int Ni,Nj;
	float x1,y1,x2,y2;
	long count = 0;
	DATA **rects = new DATA*[MAXRECT];
	
	// rtree for indexing each PtGroup
	printf("making PtGrpRTree\n");
	printf("|EdgeMap|: %d\n",EdgeMap.size());
	EdgeMapType::iterator iter=EdgeMap.begin();
 	while (iter!=EdgeMap.end()) {
 		edge* e=iter->second;
 		if (e->pts.size()>0) {
			Ni=e->Ni;	Nj=e->Nj;
			x1=xcrd[Ni];	y1=ycrd[Ni];
			x2=xcrd[Nj];	y2=ycrd[Nj];
			
			DATA *d = new DATA;
			d->data[0]=min(x1,x2);
			d->data[1]=max(x1,x2);
			d->data[2]=min(y1,y2);
			d->data[3]=max(y1,y2);
			//d->id=e->FirstRow;
			d->id=Ni;	d->id2=Nj;
			
			rects[count++] = d;
			if (count>MAXRECT) {
			    printf("memory limit reached, count=%d\n",count);
			    exit(1);
			}
		}
		iter++;
	}
	if (count==0) {
		printf("No data, count=%d\n",count);
		exit(1);
	}
	printf("loaded %d rect.\n", count);
	// rect sorted inside bulkload function
	RTree *rt = RT_bulkload(rtfile,rects,count,DEFAULT_BLOCKLENGTH);
	delete rt;
}

int ADJGRP_HEADSIZE=sizeof(int)+2*sizeof(float);
int ADJGRP_ITEMSIZE=2*sizeof(int)+sizeof(float);	

// Adj FlatFile Field:
//		Header:	xcrd(float), ycrd(float), size(int)
//		Entry:	Nk(int), eDist(float), PtGrpKey(int), PtSize(int)		changed
void makeAdjListFiles(FILE *alFile) {	
	printf("making alFiles, dependency on makePtFiles\n");
	
	int key=0,size,PtSize;
	fwrite(&NodeNum,1,sizeof(int),alFile);
	
	// slotted header info.
	int addr=sizeof(int)+sizeof(int)*NodeNum;
	for (int Ni=0;Ni<NodeNum;Ni++) {
		fwrite(&addr,1,sizeof(int),alFile);
		addr+=ADJGRP_HEADSIZE+AdjList[Ni].size()*ADJGRP_ITEMSIZE;
	}
	
	float distsum=0;
	for (int Ni=0;Ni<NodeNum;Ni++) {
		// write Ni's crd, added at 28/1/2004
		fwrite(&(xcrd[Ni]),1,sizeof(float),alFile);
		fwrite(&(ycrd[Ni]),1,sizeof(float),alFile);
		
		size=AdjList[Ni].size();
		fwrite(&(size),1,sizeof(int),alFile);
		
		for (int k=0;k<AdjList[Ni].size();k++) {
			int Nk=AdjList[Ni][k];	// Nk can be smaller or greater than Ni !!!
			edge* e=EdgeMap[getKey(Ni,Nk)];
			PtSize=e->pts.size();
			
			fwrite(&Nk,1,sizeof(int),alFile);
			fwrite(&(e->dist),1,sizeof(float),alFile);
			fwrite(&(e->FirstRow),1,sizeof(int),alFile); // use FirstRow for other purpose ... 
			//printf("(Ni,Nj,dataAddr)=(%d,%d,%d)\n",Ni,Nk,e->FirstRow);
			distsum+=e->dist;
		}
		key=Ni;
	}
	distsum=distsum/2;
	printf("total edge dist: %f\n",distsum);
}

void BuildBinaryStorage(const char* fileprefix) {
	FILE *ptFile,*edgeFile;
	char tmpFileName[255];
	BlkLen=getBlockLength();
	
	// make point group file and its bt file
	sprintf(tmpFileName,"%s.p_d",fileprefix);	remove(tmpFileName); // remove existing file
	ptFile=fopen(tmpFileName,"w+");
	sprintf(tmpFileName,"%s.p_bt",fileprefix);	remove(tmpFileName); // remove existing file
	makePtFiles(ptFile,tmpFileName);
	
	// make point group rt file // only for ANN storage model
	/*sprintf(tmpFileName,"%s.p_rt",fileprefix);	remove(tmpFileName); // remove existing file
	makePtGrpRTree(tmpFileName);*/
	
	// make adj list file
	sprintf(tmpFileName,"%s.al_d",fileprefix);	remove(tmpFileName); // remove existing file
	edgeFile=fopen(tmpFileName,"w+");
	makeAdjListFiles(edgeFile);
		
	fclose(ptFile);	fclose(edgeFile);
}

// future work: normalization of edge weight, reassign edge weight, divide into bands ...
// how about scanning whole file and putting edges with both nodeIDs in range ...
void ReadRealNetwork(char* prefix_name,int _NodeNum=0) {
	int id,Ni,Nj;
	float dist,x,y;
	
	// side effect: change nodes
	char edgef[255],nodef[255];	
	sprintf(edgef,"data/%s.cedge",prefix_name);
	sprintf(nodef,"data/%s.cnode",prefix_name);	
	
	FILE* cedge=fopen(edgef,"r");
	FILE* cnode=fopen(nodef,"r");
	CheckFile(cedge,edgef);
	CheckFile(cnode,nodef);	
	
	IS_NODE_SIZE_GIVEN=(_NodeNum>0);
	printf("Is |V| given ? %d\n",IS_NODE_SIZE_GIVEN);
	
	NodeNum=0;	// getKey depends on NodeNum so we have to init. it first
	while (!feof(cnode)) {
		fscanf(cnode,"%d %f %f\n", &id, &x, &y);
		xcrd.push_back(x);
		ycrd.push_back(y);
		NodeNum++;
	}
	if (IS_NODE_SIZE_GIVEN&&_NodeNum<=NodeNum) NodeNum=_NodeNum;
	printf("%d nodes read, ",NodeNum);		PrintElapsed();
	
	AdjList=new FastArray<int>[NodeNum];
	int EdgeNum=0;
	while (!feof(cedge)) {
		fscanf(cedge, "%d %d %d %f\n", &id, &Ni, &Nj, &dist);
		if (Ni<NodeNum&&Nj<NodeNum) {	// ignore edges outside the range
			//printf( "%d %d %d %f\n", id, Ni, Nj, dist);
			edge* e=new edge;
			
			// edge length distortion	1+ ( 0 0.125 0.25 0.5 0.75 1)
			dist=(1.0+1*drand48())*dist;
			
			e->dist=dist;
			e->Ni=min(Ni,Nj);	e->Nj=max(Ni,Nj);	// enforce the constriant Ni<Nj
			AdjList[Ni].push_back(Nj);	AdjList[Nj].push_back(Ni);
		   	EdgeMap[getKey(Ni,Nj)]=e;	// should be ok
		   	EdgeNum++;
	   	}
	}
	printf("%d edges read, ",EdgeNum);		PrintElapsed();
	fclose(cedge);	fclose(cnode);
}

void WriteTextFile(const char* filename) {
	FILE* f=fopen(filename,"w");
	fprintf(f,"%d\n",NodeNum);
	
	//fprintf(f,"// Edge Distance And Point List\n");
	EdgeMapType::iterator p=EdgeMap.begin();
	while (p!=EdgeMap.end()) {
		edge* e=p->second;
		sort(e->pts.begin(),e->pts.end());
		fprintf(f,"%d %d %f %d\n",e->Ni,e->Nj,e->dist,e->pts.size());
		
		for (int pCnt=0;pCnt<e->pts.size();pCnt++) {
			fprintf(f,"%0.3f ",e->pts[pCnt]);
			if ((pCnt+1)%100==0) fprintf(f,"\n");
		}
		if (e->pts.size()>0) fprintf(f,"\n");
		delete e;
		p++;
	}
	fclose(f);
}

void WritePointCrdFile(const char* filename) {
	// not clear => sample size too small  
	int SAMPLE_SIZE=2000;
	int INC=num_D/SAMPLE_SIZE;
	if (INC==0) INC=1;
	
	int ScanNum=0;	
	float x1,y1,x2,y2,x,y,eDist;
	FILE* f=fopen(filename,"w");
	EdgeMapType::iterator p=EdgeMap.begin();
	while (p!=EdgeMap.end()) {		
		edge* e=p->second;
		sort(e->pts.begin(),e->pts.end());
		x1=xcrd[e->Ni];	y1=ycrd[e->Ni];
		x2=xcrd[e->Nj];	y2=ycrd[e->Nj];		
		eDist=e->dist;
		for (int pCnt=0;pCnt<e->pts.size();pCnt++) {
			float vP=e->pts[pCnt];
			ScanNum++;
			if (ScanNum%INC!=0) continue;
			if (eDist==0) {x=x1;y=y1;} else {
				x=x1+(vP/eDist)*(x2-x1);
				y=y1+(vP/eDist)*(y2-y1);
			}
			fprintf(f,"%f %f\n",x,y);
		}
		p++;
	}
	fclose(f);
}

void ConnectedGraphCheck() {
	BitStore isVisited;
	isVisited.assign(NodeNum,false);
	
	StepQueue sQ;
	StepEvent stepX;
	stepX.node=0;
	sQ.push(stepX);
	
	int cnt=0;
	while (!sQ.empty()) {	// beware of infty loops
		StepEvent event=sQ.top();
		sQ.pop();
		if (isVisited[event.node]) continue;
		isVisited[event.node]=true;
		cnt++;
		int Ni=event.node;
		for (int k=0;k<AdjList[Ni].size();k++) {
			int Nk=AdjList[Ni][k];	// Nk can be smaller or greater than Ni !!!
			if (!isVisited[Nk]) {
				StepEvent newevent=event;
				newevent.node=Nk;
				sQ.push(newevent);
			}
		}
	}
	printf("connected component: %d\n",cnt);
}


float FACTOR=1;

void SpreadPoints(float STEP_SIZE,int startNode,int _NumPoint) {
	// 1*eps -> FACTOR*eps as NumPoint -> 0
	// i.e. factor>=1 must hold
	float SLOPE=(FACTOR-1)*STEP_SIZE/_NumPoint;
	float STAT_STEPSIZE=0,STAT_MAXDIST=0,STAT_MAXAMT=0;
	int NumPoint=_NumPoint;
	BitStore isVisited;
	isVisited.assign(NodeNum,false);
	
	
	int Ni,Nj;
	StepQueue sQ;
	
	Ni=startNode;	Nj=AdjList[Ni][0];
	edge* e=EdgeMap[getKey(Ni,Nj)];
	e->pts.push_back(0);	// anchor pt. of the cluster
	num_D+=1;	// set num_D to correct value
	NumPoint--;
	
	StepEvent step;
	step.node=Ni;
	step.dist=0;	step.accDist=STEP_SIZE;
	sQ.push(step);
	
	while (!sQ.empty()) {	// beware of infty loops
		StepEvent event=sQ.top();
		sQ.pop();
		
		int NodeID=event.node;
		if (isVisited[NodeID]) continue;
		isVisited[NodeID]=true;
		
		FastArray<int>& NAdjList=AdjList[NodeID];
		for (int z=0;z<NAdjList.size();z++) {
			StepEvent newevent=event;
			int NewNodeID=NAdjList[z];
			if (isVisited[NewNodeID]) continue;
			
			edge* eK=EdgeMap[getKey(NodeID,NewNodeID)];			
			
			if (newevent.accDist>eK->dist) 
				newevent.accDist-=eK->dist;
			else {
				float START_P,END_P;
				bool isLeft=(NodeID>NewNodeID);
				int INC_S;
				
				if (isLeft) {	// R to L
					START_P=eK->dist-newevent.accDist;	END_P=0;	INC_S=-1;
				} else {	// L to R
					START_P=0+newevent.accDist;	END_P=eK->dist;		INC_S=1;
				}
				newevent.accDist=0;
				while (INC_S*(END_P-START_P)>=0) {
					eK->pts.push_back(START_P);
					num_D+=1;	// set num_D to correct value
					
					float NEW_STEP=FACTOR*STEP_SIZE-SLOPE*NumPoint;		
					float STEP_AMT=NEW_STEP*(0.50+drand48());	// [0.5R,1.5R]
					if (STAT_MAXAMT<STEP_AMT) STAT_MAXAMT=STEP_AMT;					
					START_P+=STEP_AMT*INC_S;
					NumPoint--;
					if (NumPoint==0) goto SPREAD_EXIT;
					STAT_STEPSIZE=max(STAT_STEPSIZE,NEW_STEP);
				}
				if (isLeft)
					newevent.accDist=-START_P;
				else 
					newevent.accDist=START_P-eK->dist;
			}
			
			newevent.node=NewNodeID;
			newevent.dist=event.dist+eK->dist;
			sQ.push(newevent);	// propagation
			STAT_MAXDIST=max(STAT_MAXDIST,newevent.dist);
		}
	}

	SPREAD_EXIT:
	printf("R_MAX: %f/%f, dist_MAX: %f\n",STAT_STEPSIZE,STAT_MAXAMT,STAT_MAXDIST);
	if (NumPoint>0) {
		printf("Error: %d points of this cluster not generated\n",NumPoint);
		exit(0);
	}
}

void BlowupNetwork(float STEP_SIZE) {
	int Ni,Nj;
	StepQueue sQ;
	
	Ni=0;	Nj=AdjList[Ni][0];
	edge* e=EdgeMap[getKey(Ni,Nj)];
	e->pts.push_back(0);	// anchor pt. of the cluster
	num_D=1;	// set num_D to correct value
	
	StepEvent step;
	step.node=Ni;
	step.dist=0;	step.accDist=STEP_SIZE;
	sQ.push(step);
	
	BitStore isVisited;
	isVisited.assign(NodeNum,false);
	
	while (!sQ.empty()) {	// beware of infty loops
		//if (num_D>=50000) break;	// limit the number of points
		StepEvent event=sQ.top();
		sQ.pop();
		
		int NodeID=event.node;
		if (isVisited[NodeID]) continue;
		isVisited[NodeID]=true;
		
		FastArray<int>& NAdjList=AdjList[NodeID];
		for (int z=0;z<NAdjList.size();z++) {
			StepEvent newevent=event;
			int NewNodeID=NAdjList[z];
			if (isVisited[NewNodeID]) continue;
			
			edge* eK=EdgeMap[getKey(NodeID,NewNodeID)];
			
			if (newevent.accDist>eK->dist) 
				newevent.accDist-=eK->dist;
			else {
				float START_P,END_P;
				bool isLeft=(NodeID>NewNodeID);
				int INC_S;
				
				if (isLeft) {	// R to L
					START_P=eK->dist-newevent.accDist;	END_P=0;	INC_S=-1;
				} else {	// L to R
					START_P=0+newevent.accDist;		END_P=eK->dist;	INC_S=1;
				}
				newevent.accDist=0;
				while (INC_S*(END_P-START_P)>=0) {
					eK->pts.push_back(START_P);
					num_D++;
					START_P+=STEP_SIZE*INC_S;
				}
				if (isLeft)
					newevent.accDist=-START_P;
				else 
					newevent.accDist=START_P-eK->dist;
			}
			newevent.node=NewNodeID;
			newevent.dist=event.dist+eK->dist;
			sQ.push(newevent);	// propagation
		}
	}
	printf("Blowup procedure completed\n");
}

void ReadRealPoints(char* rawptfn,char* nodertrfn) {	// full path 
	printf("Need to be changed from point rt to rect(ptgrp) rt !");
	exit(0);

	RTree *rt = new RTree(nodertrfn,128);
    rt->file->isReadOnly=true;	// don't really write blocks to disk
    rt->load_root();
    
	float xc,yc,x1,y1,x2,y2,x4,y4,dist;
	int id,Ni,Nj,bestNode;
	
	FILE* cpoint=fopen(rawptfn,"r");
	CheckFile(cpoint,rawptfn);
	num_D=0;
	FastArray<float> ptxc,ptyc;
	
	while (!feof(cpoint)) {
		fscanf(cpoint,"%d %f %f %f %f\n",&id,&x1,&y1,&x2,&y2);
		ptxc.push_back((x1+x2)/2);
		ptyc.push_back((y1+y2)/2);
		num_D++;
	}
	fclose(cpoint);
	printf("%d points read, ",num_D);		PrintElapsed();
	
	DATA NNpt,querypt;
    BranchQueue INN_queue;
	BranchEvent event;
	int rt_totna=0;
	
	for (int p=0;p<num_D;p++) {
		event.pnode=NULL;
		event.son=-1;
		event.mindist=0;	// always search from root node
		INN_queue.push(event);
		getNextNN(rt,INN_queue,rt_totna,&querypt,&NNpt);
		
		xc=querypt.data[0]=ycrd[p];
		yc=querypt.data[1]=xcrd[p];
		while (!INN_queue.empty()) INN_queue.pop();		// clear queue
		
		// can find nearest node, how about its position on node
		// nearest edge and projected dist. on it
		Ni=NNpt.id;
		edge* bestE=NULL;
		float best_ndist=MAX_DIST,best_vP=0;
		for (int j=0;j<AdjList[Ni].size();j++) {
			int Nj=AdjList[Ni][j];	// Nk can be smaller or greater than Ni !!!
			edge* e=EdgeMap[getKey(Ni,Nj)];
			dist=e->dist;
			float x1=xcrd[Ni],y1=ycrd[Ni],x2=xcrd[Nj],y2=ycrd[Nj];
			float vP=((xc-x1)*(x2-x1)+(yc-y1)*(y2-y1))/dist;
			vP=max(0,vP);		vP=min(vP,dist);
			x4=x1+(vP/dist)*(x2-x1);
			y4=y1+(vP/dist)*(y2-y1);
			float nd=getDist(x4,y4,xc,yc);
			if (nd<best_ndist) {
				best_ndist=nd;	best_vP=vP;		bestE=e;
			}
		}
		bestE->pts.push_back(best_vP);
		
		if (p%1000==0) {
			printf("point %d ",p); PrintElapsed();
		}
	}
	delete rt;
}


void GenFarClusters(int CLUS_NUM,int CLUS_SZ,float STEP_SIZE) {
	int seedID[CLUS_NUM];
	
	for (int i=0;i<CLUS_NUM;i++) {
		do {
			seedID[i]=rand()%NodeNum;
		} while (seedID[i]==0);	// exclude the case for i>0 !
		
		
		/*if (i>0) {
			float gmMinDist=0;
			for (int round=0;round<10;round++) {
				int nid=rand()%NodeNum;				
				float mdist=MAXREAL;
				for (int z=0;z<i;z++) {
					float tdist=getDist(xcrd[nid],ycrd[nid],xcrd[seedID[z]],ycrd[seedID[z]]);
					if (tdist<mdist) mdist=tdist;
				}
								
				if (mdist>gmMinDist) {	// i.e. mindist to others quite large !
					gmMinDist=mdist;
					seedID[i]=nid;		// mindist to other furthest !
				}
			}
		}*/
		
		SpreadPoints(STEP_SIZE,seedID[i],CLUS_SZ);	// default 50000 points
	}	
}


void ReadRealDataPoints(char* rawptfn) {	// full path 
	float xc,yc,x1,y1,x2,y2,x4,y4,dist;
	int id,Ni,Nj,bestNode;
	
	FILE* cpoint=fopen(rawptfn,"r");
	CheckFile(cpoint,rawptfn);
	num_D=0;
	
	while (!feof(cpoint)) {
		num_D++;
		fscanf(cpoint,"%d %f %f %f %f\n",&id,&x1,&y1,&x2,&y2);
		xc=(x1+x2)/2;
		yc=(y1+y2)/2;
		
		// can find nearest node, how about its position on node
		// nearest edge and projected dist. on it
		Ni=0;	//NNpt.id;
		float bestEucDist=MAXREAL;
		for (int idx=0;idx<NodeNum/10;idx++) { // cheaper: only test a set of random nodes
			int i=rand()%NodeNum;
			float simp_dist=fabs(xcrd[i]-xc)+fabs(ycrd[i]-yc);
			if (simp_dist<bestEucDist) {
				bestEucDist=simp_dist;
				Ni=i;
			}
		}
		
		edge* bestE=NULL;
		float best_ndist=MAX_DIST,best_vP=0;
		for (int j=0;j<AdjList[Ni].size();j++) {
			int Nj=AdjList[Ni][j];	// Nk can be smaller or greater than Ni !!!
			edge* e=EdgeMap[getKey(Ni,Nj)];
			dist=e->dist;
			float x1=xcrd[Ni],y1=ycrd[Ni],x2=xcrd[Nj],y2=ycrd[Nj];
			float vP=((xc-x1)*(x2-x1)+(yc-y1)*(y2-y1))/dist;
			vP=max(0,vP);		vP=min(vP,dist);
			x4=x1+(vP/dist)*(x2-x1);
			y4=y1+(vP/dist)*(y2-y1);
			float nd=getDist(x4,y4,xc,yc);
			if (nd<best_ndist) {
				best_ndist=nd;	best_vP=vP;		bestE=e;
			}
		}
		bestE->pts.push_back(best_vP);

		
		if (num_D%1000==0) {
			printf("point %d ",num_D); PrintElapsed();
		}
	}
	fclose(cpoint);
	printf("%d points read, ",num_D);		PrintElapsed();
}

/* Generation parameters:
	[Complusory]
	|V|			int, number of nodes
	|D|			int, number of points
	
	[Optional]
	Step 		float, step size for spreading point uniformly */
void GenDataMain(int argc, char *argv[]) {
	if (argc<4||argc>5) {
		printf("Usage: %s <file> <|V|> <|D|> [G/W]\n",argv[0]);
		exit(0);
	}
	
  	char filename[255],fileprefix[255];
	int _NodeNum = atoi(argv[2]);	// not used now !
	num_D = atoi(argv[3]);
	
	float step_p=1.0;
	if (argc>=5)	step_p=atof(argv[4]);
	
	sprintf(fileprefix,"map/%s",argv[1]);
	sprintf(filename,"%s.map",fileprefix);
	
	// [optimized maps, should use one of these]
	char* mapname="SF_dnw";
	
	ReadRealNetwork(mapname,_NodeNum);
	ConnectedGraphCheck();	// checking whether connected graph	
	
	//ReadRealDataPoints("data/cn_hspt.txt");
	float avg_length=0;	
	EdgeMapType::iterator iter=EdgeMap.begin();
 	while (iter!=EdgeMap.end()) {
 		edge* e=iter->second;
 		avg_length+=e->dist;
		iter++;
	}
	avg_length=avg_length/EdgeMap.size();
	{	num_D=0;	// assume real real points not invoked !
		
		BlowupNetwork(step_p*avg_length);
		
//		float STEP_SIZE=step_p*avg_length;
//		int CLUS_NUM=100;
//		int CLUS_SZ=1000;
//		FACTOR=1;
//		GenFarClusters(CLUS_NUM,CLUS_SZ,STEP_SIZE);
		printf("%d points generated\n",num_D);
	}
	BuildBinaryStorage(fileprefix);	// writing the network and points into disk
	
	char CrdFile[255];
	sprintf(CrdFile,"visual/%s.pdat",argv[1]);
	WritePointCrdFile(CrdFile);
	//WriteTextFile(filename);
}

/* Generation parameters:
	[Complusory]
	|Q|			int, number of query points
	P%			Subnetwork size (# nodes) as percentage of nodes in whole network*/
void GenQueryMain(int argc, char *argv[]) {
	FastArray<int> NodeChoices;
	StepQueue sQ;
	BitStore isVisited;
	StepEvent step;
	
	if (argc!=5) {
		printf("Usage: %s <query-file> <map-files> <|Q|> <A%%>\n",argv[0]);
		exit(0);
	}
	
	int seedNum=atoi(argv[3]);
	float pcent=atof(argv[4]);
	ReadRealNetwork(argv[2]);
	step.node=rand()%NodeNum;
	step.dist=0;
	sQ.push(step);
	isVisited.assign(NodeNum,false);
	
	while (!sQ.empty()) {	// beware of infty loops
		StepEvent event=sQ.top();
		sQ.pop();
		
		int NodeID=event.node;
		if (isVisited[NodeID]) continue;
		isVisited[NodeID]=true;
		
		if (NodeChoices.size()>=NodeNum*pcent) break;
		NodeChoices.push_back(NodeID);
		
		FastArray<int>& NAdjList=AdjList[NodeID];
		for (int z=0;z<NAdjList.size();z++) {
			StepEvent newevent=event;
			int NewNodeID=NAdjList[z];
			if (isVisited[NewNodeID]) continue;
			edge* eK=EdgeMap[getKey(NodeID,NewNodeID)];
			newevent.node=NewNodeID;
			newevent.dist=event.dist+eK->dist;
			sQ.push(newevent);	// propagation
		}
	}
	
	// write query points to the file
	// additional function to sort query points for better network locality ?
	char queryf[255];
	sprintf(queryf,"query/%s.qry",argv[1]);
	FILE* cquery=fopen(queryf,"w");
	CheckFile(cquery,queryf);
	fprintf(cquery,"%d\n",seedNum);
	for (int s=0;s<seedNum;s++) {
		int Ni,Nj,Ri,Rj;
		int pick=rand()%NodeChoices.size();
		pick=(int)((s+0.5)*NodeChoices.size()/seedNum);
		
		Ri=NodeChoices[pick];
		FastArray<int>& NAdjList=AdjList[Ri];
		Rj=NAdjList[rand()%NAdjList.size()];
		Ni=min(Ri,Rj);	Nj=max(Ri,Rj);
		edge* eK=EdgeMap[getKey(Ni,Nj)];
		float vP=drand48()*eK->dist;
		fprintf(cquery,"%d %d %f\n",Ni,Nj,vP);
	}
	fclose(cquery);
}

int main(int argc, char *argv[]) {
	InitClock();	// side effect: init. seeds for random generators
	
	ConfigType cr;
	AddConfigFromFile(cr,"config.prop");
	AddConfigFromCmdLine(cr,argc,argv);
	
	const char* generator=getConfigStr("generator",cr);
	printf("generator: %s\n",generator);
	srand(0);	srand48(0);		// fix the query points to be generated !!!
	
	if (strcmp(generator,"data")==0)
		GenDataMain(argc,argv);
	
	if (strcmp(generator,"query")==0)
		GenQueryMain(argc,argv);
	
	PrintElapsed();
	return 0;
}
