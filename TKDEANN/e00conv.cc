//              e00 Converter v1.0
//
// This program converts the .e00 (Arc/Info Export File Format) to generic data format
//
// Coded by Zhang Jun (Bobby)  Jan. 27, 2001
//
// The input file follows the .e00 data format
// 
// The format of the output file
// <id, x1, y1, x2, y2> (only at the end of writing ! )

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#define MAXOBJ 1000000
#define LEN 1000

#define UNIVERSE 10000

#define PROCESS_ARC true // whether of not process ARC (arc) objects
#define PROCESS_LAB true // whether of not process LAB (label) objects
#define min(a, b) (((a) < (b))? (a) : (b)  )
#define max(a, b) (((a) > (b))? (a) : (b)  )

int *id;
float *mbr;
	
void PolylineToLines(FILE* i_fp,int& count,int num_coordinates) {
	// read the starting point
	float oldx, oldy,x, y;
	fscanf(i_fp, "%e%e", &oldx, &oldy);
	for (int i = 0; i < num_coordinates - 1; i++) {
		fscanf(i_fp, "%e%e", &x, &y);
		id[count] = count;
		mbr[4*count] = oldx;
		mbr[4*count+1] = x;
		mbr[4*count+2] = oldy;
		mbr[4*count+3] = y;
		count++;
		oldx = x;	oldy = y;
	}
}

void PolylineToPoints(FILE* i_fp,int& count,int num_coordinates) {
	// read the starting point
	float oldx, oldy;
	
	// new version
	for (int i = 0; i < num_coordinates; i++) {
		fscanf(i_fp, "%e%e", &oldx, &oldy); 
		id[count] = count;
		mbr[4*count] = oldx;
		mbr[4*count+1] = oldx;
		mbr[4*count+2] = oldy;
		mbr[4*count+3] = oldy;
		count++;
	}
	
	/*fscanf(i_fp, "%e%e", &oldx, &oldy); 
	if (num_coordinates == 1) {
		id[count] = count;
		mbr[4*count] = oldx;
		mbr[4*count+1] = oldx;
		mbr[4*count+2] = oldy;
		mbr[4*count+3] = oldy;
		count++;
	} else {
		// compute the mid point of each edge
		float x, y;
		for (int i = 0; i < num_coordinates - 1; i++) {
			fscanf(i_fp, "%e%e", &x, &y);
			id[count] = count;
			mbr[4*count] = (oldx + x) / 2;
			mbr[4*count+1] = (oldx + x) / 2;
			mbr[4*count+2] = (oldy + y) / 2;
			mbr[4*count+3] = (oldy + y) / 2;
			count++;
			oldx = x;	oldy = y;
		}
	}*/
}

int main(int argc, char* argv[]) {
	char i_fname[300] = "";
	char o_fname[300] = "";
	
	printf("e00 Converter v1.0\n");
	if (argc!=3) {
		printf("%s <input> <output>",argv[0]);
		return 0;
	}

	char* INPUT_FILENAME=argv[1];
	char* OUTPUT_FILENAME=argv[2];

	strcpy(i_fname, INPUT_FILENAME);
	strcpy(o_fname, OUTPUT_FILENAME);
	
	FILE* i_fp = NULL;
	FILE* o_fp = NULL;
	
	i_fp = fopen(i_fname, "r");
	o_fp = fopen(o_fname, "w");

	char EXP[LEN], filename[LEN];
	int compression;
	fscanf(i_fp, "%s%d%s", EXP, &compression, filename);
	if (strcmp(EXP, "EXP") != 0) {
		printf("This is not a valid .e00 file.\n.e00 file should start with 'EXP'.\n");
		exit(-1);
	}
	printf("The original filename: %s\n", filename);
	if (compression!=0) {
		printf("Compression is not supported by this converter.\n");
		exit(-1);
	}
	
	int count = 0;
	int oldcount = count;
	id = new int[MAXOBJ];
	mbr = new float[4*MAXOBJ];
	int i, read;
	char format[LEN];
	
	while (true) {
		read = fscanf(i_fp, "%s", format);
		if (read!=1) break;
		
		bool ARC_processed = !PROCESS_ARC;
		bool LAB_processed = !PROCESS_LAB;
		
		oldcount = count;
		if (strcmp(format, "ARC") == 0 && !ARC_processed) {
			ARC_processed = true;
			printf("ARC is being processed...\n");
			
			int precision;
			read = fscanf(i_fp, "%d", &precision);
			
			if (precision == 2)
				printf("Single precision\n");
			else if (precision == 3)
				printf("Double precision\n");
			else {
				printf("This is not a valid .e00 file.\nPrecision is invalid.\n");
				exit(-1);
			}
			
			int coveragesn, coverageid, fromnode, tonode, leftpolygon, 
				rightpolygon, num_coordinates;
			
			while (true) {
				read = fscanf(i_fp, "%d%d%d%d%d%d%d", &coveragesn, &coverageid, 
					&fromnode, &tonode, &leftpolygon, &rightpolygon, &num_coordinates);
					
				if (read != 7) {
					printf("This is not a valid .e00 file.\nARC parameters are invalid.\n");
					exit(-1);
				}

				if (coveragesn == -1) break;
				/* use the center of the line as a point */	
				
				PolylineToPoints(i_fp,count,num_coordinates);	// polyline => points
				//PolylineToLines(i_fp,count,num_coordinates); // polyline => lines
			}
			printf("ARC was processed. %d objects were read.\n", count - oldcount);
		} else if (strcmp(format, "LAB") == 0 && !LAB_processed) {
			LAB_processed = true;
			printf("LAB is being processed...\n");
			int precision;
			read = fscanf(i_fp, "%d", &precision);
			
			if (precision == 2)
				printf("Single precision\n");
			else if (precision == 3)
				printf("Double precision\n");
			else {
				printf("This is not a valid .e00 file.\nPrecision is invalid.\n");
				exit(-1);
			}

			int coverageid, polygonid;
			float x, y, d0, d1, d2, d3;
			while (true) {
				read = fscanf(i_fp, "%d%d%f%f%f%f%f%f", &coverageid, &polygonid, 
					&x, &y, &d0, &d1, &d2, &d3);

				if (read != 8 && read != 2 && (read == 2 && coverageid != -1)) {
					printf("This is not a valid .e00 file.\nARC parameters are invalid.\n");
					exit(-1);
				}
				if (coverageid == -1) break;
				id[count] = count;
				mbr[4*count] = x;
				mbr[4*count+1] = x;
				mbr[4*count+2] = y;
				mbr[4*count+3] = y;	
				count++;
			}
			printf("LAB was processed. %d objects were read.\n", count - oldcount);
		}
	}

	// find the min and max for normalization
	float minx = 1e30f,maxx = -1e30f;
	float miny = 1e30f,maxy = -1e30f;
	for (i = 0; i < count; i++) { // note: we can have lines ( i.e. not in MBR format)
		float curminx=min(mbr[4*i],mbr[4*i+1]);
		float curmaxx=max(mbr[4*i],mbr[4*i+1]);
		float curminy=min(mbr[4*i+2],mbr[4*i+3]);
		float curmaxy=max(mbr[4*i+2],mbr[4*i+3]);
	
		if (curminx < minx)  minx = curminx;
		if (curmaxx > maxx) maxx = curmaxx;
		if (curminy < miny) miny = curminy;
		if (curmaxy > maxy) maxy = curmaxy;
	}
	
	// perform normalization
	for (i = 0; i < count; i++) {
		mbr[4*i] =   (mbr[4*i] - minx ) / (maxx - minx) * UNIVERSE;
		mbr[4*i+1] = (mbr[4*i+1] - minx ) / (maxx - minx) * UNIVERSE;
		mbr[4*i+2] = (mbr[4*i+2] - miny ) / (maxy - miny) * UNIVERSE;
		mbr[4*i+3] = (mbr[4*i+3] - miny ) / (maxy - miny) * UNIVERSE;
	}
	
	// output the objects in the order: x1,y1,x2,y2
	for (i = 0; i < count; i++) 
		fprintf(o_fp, "%d %f %f %f %f\n",id[i],mbr[4*i],mbr[4*i+2],mbr[4*i+1],mbr[4*i+3]);
	printf("norm. : %f %f %f %f\n",minx,miny,maxx,maxy);
	printf("Total %d objects were processed.\n", count);	
	delete []mbr;	delete []id;
	fclose(o_fp);	fclose(i_fp);
	return 0;
}
