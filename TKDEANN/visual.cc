#include "utility.h"
#include "netshare.h"

FastArray<float> xcrd,ycrd;

void NETtoPLOT(char* edgef,char* nodef,char* visualf) {
	// side effect: change nodes
	FILE* cedge=fopen(edgef,"r");
	FILE* cnode=fopen(nodef,"r");
	FILE* cview=fopen(visualf,"w");
	CheckFile(cedge,edgef);
	CheckFile(cnode,nodef);
	CheckFile(cview,visualf);
	
	int id,Ni,Nj;
	float dist,x,y;	
	
	// read node info with format: NodeId xcrd ycrd
	NodeNum=0;
	while (!feof(cnode)) {
		fscanf(cnode,"%d %f %f\n", &id, &x, &y);
		xcrd.push_back(x);	ycrd.push_back(y);
		NodeNum++;
	}
	printf("Nodes read, ");		PrintElapsed();
	
	// read edge info with format: EdgeId Ni Nj eDist
	while (!feof(cedge)) {
		fscanf(cedge, "%d %d %d %f\n", &id, &Ni, &Nj, &dist);
		fprintf(cview,"%f %f\n%f %f\n\n",xcrd[Ni],ycrd[Ni],xcrd[Nj],ycrd[Nj]);
	}
	fclose(cedge);
	fclose(cnode);
	fclose(cview);
}

void RAWtoPLOT(char* rawedgef,char* visualf) {
	FILE* crawedge=fopen(rawedgef,"r");
	FILE* cview=fopen(visualf,"w");
	CheckFile(crawedge,rawedgef);
	CheckFile(cview,visualf);
	float id,x1,y1,x2,y2,x,y;
	while (!feof(crawedge)) {
		fscanf(crawedge, "%d %f %f %f %f\n", &id, &x1, &y1, &x2,&y2);
		x=(x1+x2)/2;	y=(y1+y2)/2;
		//fprintf(cview,"%f %f\n",x,y);	// to point
		fprintf(cview,"%f %f\n%f %f\n\n",x1,y1,x2,y2);	// to line
	}
	fclose(crawedge);
	fclose(cview);
}

int main(int argc, char *argv[]) {
	InitClock();	// side effect: init. seeds for random generators
	if (argc!=2) {
		printf("Usage: %s <map-prefix>\n",argv[0]);
		return 0;
	}
	
//	RAWtoPLOT("data/cnrrline.txt","visual/cnrrline.net");
//	exit(0);
	
	char edgef[255],nodef[255],map_prefix[255],visualf[255];
	sprintf(map_prefix,"%s",argv[1]);
	sprintf(edgef,"data/%s.cedge",map_prefix);
	sprintf(nodef,"data/%s.cnode",map_prefix);
	sprintf(visualf,"visual/%s.net",map_prefix);
	NETtoPLOT(edgef,nodef,visualf);
	
	PrintElapsed();
	return 0;
}
