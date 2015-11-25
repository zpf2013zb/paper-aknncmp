#include "rtree.h"
#include "utility.h"
#include "diskbased.h"

// change here when ...
#define MAXRECT  3125000

// later change to ID of the point group ??
// sequentially access a point file, and build this
void LoadRectByPointGroup(char *nodef,DATA **rects, long *count) {
	int id,Ni,Nj,GroupAddr,PtGrpSize,r;
	float x,y,x1,y1,x2,y2,dist;
	DATA *d;
	
	FILE* cnode=fopen(nodef,"r");
	CheckFile(cnode,nodef);
	
	while (!feof(cnode)) {
		fscanf(cnode,"%d %f %f\n", &id, &x, &y);
		xcrd.push_back(x);	ycrd.push_back(y);
	}	
	
	r=0;
	GroupAddr=0;
	while (r<PtNum) {	// for rect. R-tree
		// read header	method: scanning  (use cache for linear scan ?)
		getFixedF(NI_P,Ref(Ni),GroupAddr);
		getFixedF(NJ_P,Ref(Nj),GroupAddr);
		getFixedF(SIZE_P,Ref(PtGrpSize),GroupAddr);
		
		x1=xcrd[Ni];	y1=ycrd[Ni];
		x2=xcrd[Nj];	y2=ycrd[Nj];
		
		DATA *d = new DATA;
		d->data[0]=min(x1,x2);
		d->data[1]=max(x1,x2);
		d->data[2]=min(y1,y2);
		d->data[3]=max(y1,y2);
		d->id=Ni;	d->id2=Nj;
		
		rects[(*count)++] = d;
		if (*count>MAXRECT) {
		    printf("memory limit reached, count=%d\n", *count);
		    exit(1);
		}
		
		//printf("%d %d %d\n",Ni,Nj,PtGrpSize);
		r+=PtGrpSize;
		GroupAddr+=PTGRP_HEADSIZE+PtGrpSize*PTGRP_ITEMSIZE;
	}
	fclose(cnode);
}

// need to be changed: read pt grp and make rt file !!!
void LoadRectsByEdges(char *nodef,char *edgef,DATA **rects, long *count) {
	int id,Ni,Nj;
	float x,y,x1,y1,x2,y2,dist;
	DATA *d;
	
	FILE* cnode=fopen(nodef,"r");
	CheckFile(cnode,nodef);
	
	while (!feof(cnode)) {
		fscanf(cnode,"%d %f %f\n", &id, &x, &y);
		xcrd.push_back(x);	ycrd.push_back(y);
	}	
	
	FILE* cedge=fopen(edgef,"r");
	CheckFile(cedge,edgef);	
	while (!feof(cedge)) {	// for rect. R-tree
		fscanf(cedge, "%d %d %d %f\n", &id, &Ni, &Nj, &dist);
		x1=xcrd[Ni];	y1=ycrd[Ni];
		x2=xcrd[Nj];	y2=ycrd[Nj];
		
		DATA *d = new DATA;
		d->data[0]=min(x1,x2);
		d->data[1]=max(x1,x2);
		d->data[2]=min(y1,y2);
		d->data[3]=max(y1,y2);
		d->id=Ni;	d->id2=Nj;
		
		rects[(*count)++] = d;
		if (*count>MAXRECT) {
		    printf("memory limit reached, count=%d\n", *count);
		    exit(1);
		}
	}
	fclose(cnode);	fclose(cedge);
}

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
//	for (int i=0; i<count; i++)	rects[i]->print();	
//	return 0;
	
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
            offset+= 2*DIMENSION*sizeof(float)+2*sizeof(int);
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

int main(int argc, char* argv[]) {
	char nodef[255],edgef[255],ptprefix[255],rtfile[255];
    DATA **rects = new DATA*[MAXRECT];
    long count = 0;
    
    if (argc<2||argc>3)
    	printf("Usage: %s <map> [point] \n", argv[0]);
  	else {
		sprintf(nodef,"data/%s.cnode",argv[1]);
    	if (argc==2) {
    		sprintf(edgef,"data/%s.cedge",argv[1]);
    		sprintf(rtfile,"data/%s.ed_rt",argv[1]);
    		LoadRectsByEdges(nodef,edgef,rects,&count);
    	} else {
    		sprintf(ptprefix,"map/%s",argv[2]);
    		sprintf(rtfile,"map/%s.p_rt",argv[2]);
    		OpenDiskComm(ptprefix,DEFAULT_CACHESIZE);
    		LoadRectByPointGroup(nodef,rects,&count);
    		CloseDiskComm();
    	}
        printf("loaded %d rect.\n", count);
        
        // rect sorted inside bulkload function
        printf("output R-tree file: %s\n",rtfile);
        remove(rtfile);			// remove file if exists
		RTree *rt = RT_bulkload(rtfile,rects,count,DEFAULT_BLOCKLENGTH);
        delete rt;
    }
   return 0;
}
