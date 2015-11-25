#include "rtree.h"

#define CHECK_MBR

//////////////////////////////////////////////////////////////////////////////
// C-functions
//////////////////////////////////////////////////////////////////////////////
void check_mbr(int dimension, float *mbr)
{
}

float area(int dimension, float *mbr) {
    float sum = 1.0;
    for (int i = 0; i < dimension; i++)
	sum *= mbr[2*i+1] - mbr[2*i];
    return sum;
}

float margin(int dimension, float *mbr) {
    float *ml, *mu, *m_last, sum;
    sum = 0.0;
    m_last = mbr + 2*dimension;
    ml = mbr;
    mu = ml + 1;
    while (mu < m_last) {
		sum += *mu - *ml;
		ml += 2;
		mu += 2;
    }
    return sum;
}

bool inside(float &p, float &lb, float &ub) {
    return (p >= lb && p <= ub);
}

bool inside(float *v, float *mbr, int dimension) {
    for (int i = 0; i < dimension; i++)
		if (!inside(v[i], mbr[2*i], mbr[2*i+1]))
	    	return FALSE;
    return TRUE;
}

float overlap(int dimension, float *r1, float *r2)
// calcutales the overlapping area of r1 and r2
// calculate overlap in every dimension and multiplicate the values
{
    float sum;
    float *r1pos, *r2pos, *r1last, r1_lb, r1_ub, r2_lb, r2_ub;

    sum = 1.0;
    r1pos = r1; r2pos = r2;
    r1last = r1 + 2 * dimension;

    while (r1pos < r1last)
    {
	r1_lb = *(r1pos++);
	r1_ub = *(r1pos++);
	r2_lb = *(r2pos++);
	r2_ub = *(r2pos++);

        // calculate overlap in this dimension

        if (inside(r1_ub, r2_lb, r2_ub))
        // upper bound of r1 is inside r2
	{
            if (inside(r1_lb, r2_lb, r2_ub))
            // and lower bound of r1 is inside
                sum *= (r1_ub - r1_lb);
            else
                sum *= (r1_ub - r2_lb);
	}
	else
	{
            if (inside(r1_lb, r2_lb, r2_ub))
	    // and lower bound of r1 is inside
		sum *= (r2_ub - r1_lb);
	    else
	    {
		if (inside(r2_lb, r1_lb, r1_ub) &&
		    inside(r2_ub, r1_lb, r1_ub))
	        // r1 contains r2
		    sum *= (r2_ub - r2_lb);
		else
		// r1 and r2 do not overlap
		    sum = 0.0;
	    }
	}
    }

    return sum;
}

void enlarge(int dimension, float **mbr, float *r1, float *r2) {
	// enlarge r in a way that it contains s
    *mbr = new float[2*dimension];
    for (int i = 0; i < 2*dimension; i += 2) {
		(*mbr)[i]   = min(r1[i],   r2[i]);
		(*mbr)[i+1] = max(r1[i+1], r2[i+1]);
    }
#ifdef CHECK_MBR
    check_mbr(dimension,*mbr);
#endif
}

bool section(int dimension, float *mbr1, float *mbr2) {
    for (int i = 0; i < dimension; i++) {
		if (mbr1[2*i]>mbr2[2*i+1]||mbr1[2*i+1]<mbr2[2*i])
	    	return FALSE;
    }
    return TRUE;
}

int sort_lower_mbr(const void *d1, const void *d2) {
    SortMbr *s1, *s2;
    float erg;
    int dimension;

    s1 = (SortMbr *) d1;
    s2 = (SortMbr *) d2;
    dimension = s1->dimension;
    erg = s1->mbr[2*dimension] - s2->mbr[2*dimension];
    if (erg < 0.0)
		return -1;
    else if (erg == 0.0)
		return 0;
    else
		return 1;
}

int sort_upper_mbr(const void *d1, const void *d2) {
    SortMbr *s1, *s2;
    float erg;
    int dimension;
    
    s1 = (SortMbr *) d1;
    s2 = (SortMbr *) d2;
    dimension = s1->dimension;
    erg = s1->mbr[2*dimension+1] - s2->mbr[2*dimension+1];
    if (erg < 0.0)
		return -1;
    else if (erg == 0.0)
		return 0;
    else
		return 1;
}

int sort_center_mbr(const void *d1, const void *d2) {
    SortMbr *s1, *s2;
    int i, dimension;
    float d, e1, e2;

    s1 = (SortMbr *) d1;
    s2 = (SortMbr *) d2;
    dimension = s1->dimension;

    e1 = e2 = 0.0;
    for (i = 0; i < dimension; i++) {
        d = ((s1->mbr[2*i] + s1->mbr[2*i+1]) / 2.0) - s1->center[i];
        e1 += d*d;
        d = ((s2->mbr[2*i] + s2->mbr[2*i+1]) / 2.0) - s2->center[i];
        e2 += d*d;
    }
    
    if (e1 < e2)
		return -1;
    else if (e1 == e2)
		return 0;
    else
		return 1;
}

// changed for rectangle DATA
float MINDIST(DATA* p, float *bounces) {
    float diff,summe = 0.0;
    for (int i = 0; i < DIMENSION; i++) {
    	diff=0;
		if (p->data[2*i+1] < bounces[2*i])
		    diff = bounces[2*i]-p->data[2*i+1];
		else if (p->data[2*i] > bounces[2*i +1])
			diff = p->data[2*i]-bounces[2*i+1];
		summe += pow(diff,2);
    }
    return(summe);
}


//////////////////////////////////////////////////////////////////////////////
// DATA
//////////////////////////////////////////////////////////////////////////////

DATA::DATA() {
	data = new float[2*DIMENSION];
	id=-1;	id2=-1;
}

DATA::~DATA() {
    delete [] data;
}

DATA & DATA::operator = (DATA &_d) {
    memcpy(data, _d.data, 2* sizeof(float) * DIMENSION);
    id=_d.id;	id2=_d.id2;
    return *this;
}

void DATA::read_from_buffer(char *buffer) {
    int i;
    i = 2*DIMENSION*sizeof(float);
    memcpy(data, buffer, i);
    memcpy(&id, &buffer[i], sizeof(int));
    i+=sizeof(int);
    memcpy(&id2, &buffer[i], sizeof(int));
}

void DATA::write_to_buffer(char *buffer) {
    int i;
    i = 2*DIMENSION*sizeof(float);
    memcpy(buffer, data, i);
    memcpy(&buffer[i], &id, sizeof(int));
    i+=sizeof(int);
    memcpy(&buffer[i], &id2, sizeof(int));
}

int DATA::get_size() {
    return (2 * DIMENSION * sizeof(float) + 2*sizeof(int));
}

void DATA::print() {
	printf("%d %d ( ",id,id2);
	for (int i = 0; i < DIMENSION; i++)
		printf("[%7.5f,%7.5f]",data[2*i],data[2*i+1]);
	printf(")\n");
}

float *DATA::get_mbr() {
	float *f = new float[2*DIMENSION];
    for (int i = 0; i < DIMENSION; i++) {
		f[2*i] = data[2*i];
		f[2*i+1] = data[2*i+1];
	}
    return f;
}

//////////////////////////////////////////////////////////////////////////////
// DirEntry
//////////////////////////////////////////////////////////////////////////////

DirEntry ::DirEntry(int _dimension, RTree  *rt) {
    if (_dimension != 0) {
        dimension = _dimension;
        my_tree = rt;
        bounces = new float[2*dimension];
    }
    son_ptr = NULL;
}

DirEntry ::~DirEntry() {
    delete [] bounces;
    if (son_ptr != NULL)
	delete son_ptr;
}

DirEntry & DirEntry ::operator = (DirEntry &_d) {
    dimension = _d.dimension;
    son = _d.son;
    son_ptr = _d.son_ptr;
    son_level = _d.son_level;
    son_is_data = _d.son_is_data;
    memcpy(bounces, _d.bounces, sizeof(float) * 2 * dimension);
    return *this;
}

bool DirEntry ::is_inside(float *v) {
    for (int i = 0; i < dimension; i++) {
	if (v[i] < bounces[2*i] || v[i] > bounces[2*i + 1])   
	    return FALSE;
    }
    return TRUE;
}

SECTION DirEntry ::section(float *mbr) {
    bool inside,overlap;
    int i;

    overlap = TRUE;
    inside = TRUE;
    for (i = 0; i < dimension; i++) {
		if (mbr[2*i]     > bounces[2*i + 1] ||
		    mbr[2*i + 1] < bounces[2*i])
		    overlap = FALSE;
		if (mbr[2*i]     < bounces[2*i] ||
		    mbr[2*i + 1] > bounces[2*i + 1])
		    inside = FALSE;
    }
    if (inside)
		return INSIDE;
    else if (overlap)
		return OVERLAP;
    else
		return S_NONE;
}

void DirEntry ::read_from_buffer(char *buffer) {
    int i;
    i = 2*dimension*sizeof(float);
    memcpy(bounces, buffer, i);
    memcpy(&son, &buffer[i], sizeof(int));
}

void DirEntry ::write_to_buffer(char *buffer) {
    int i;
    i = 2*dimension*sizeof(float);
    memcpy(buffer, bounces, i);
    memcpy(&buffer[i], &son, sizeof(int));
}

int DirEntry ::get_size() {
    return 2*dimension * sizeof(float) + sizeof(int);
}

RTNode * DirEntry ::get_son() {
    if (son_ptr == NULL) {
		if (son_is_data) {
			//printf("data page\n");
		    son_ptr = new RTDataNode (my_tree, son);
		} else {
			//printf("dir page\n");
		    son_ptr = new RTDirNode (my_tree, son);
		}
    }
    return son_ptr;
}

bool DirEntry ::section_circle(DATA *center, float radius) {
	float r2;
	r2=radius*radius;
	if ((r2-MINDIST(center,bounces )) < 1.0e-8)
		return TRUE;
	else
		return FALSE;
}


//////////////////////////////////////////////////////////////////////////////
// RTNode
//////////////////////////////////////////////////////////////////////////////

RTNode ::~RTNode() {
}

RTNode ::RTNode(RTree *rt) {
    my_tree = rt;
    dimension = rt->dimension;
    num_entries = 0;
    block = -1;
}


int RTNode ::split(float **mbr, int **distribution) {
#ifdef SHOWMBR
	split_000++;
#endif

    bool lu;
    int i, j, k, l, s, n, m1, dist, split_axis;
    SortMbr *sml, *smu;
    float minmarg, marg, minover, mindead, dead, over,
	*rxmbr, *rymbr;

    // how much nodes are used?
    n = get_num();

    // nodes have to be filled at least 40%
    m1 = (int) ((float)n * 0.40);

    // sort by lower value of their rectangles
    // Indexarray aufbauen und initialisieren
    sml = new SortMbr[n];
    smu = new SortMbr[n];
    rxmbr = new float[2*dimension];
    rymbr = new float[2*dimension];

    // choose split axis
    minmarg = MAXREAL;
    for (i = 0; i < dimension; i++) {  // for each axis
        for (j = 0; j < n; j++) {
            sml[j].index = smu[j].index = j;
            sml[j].dimension = smu[j].dimension = i;
            sml[j].mbr = smu[j].mbr = mbr[j];
        }

        // Sort by lower and upper value perpendicular axis_i
      	qsort(sml, n, sizeof(SortMbr), sort_lower_mbr);
        qsort(smu, n, sizeof(SortMbr), sort_upper_mbr);

        marg = 0.0;
        // for all possible distributions of sml
        for (k = 0; k < n - 2*m1 + 1; k++) {
            // now calculate margin of R1
		    // initialize mbr of R1
		    for (s = 0; s < 2*dimension; s += 2) {
				rxmbr[s] =    MAXREAL;
				rxmbr[s+1] = -MAXREAL;
		    }
            for (l = 0; l < m1+k; l++) {
                // calculate mbr of R1
				for (s = 0; s < 2*dimension; s += 2) {
		            rxmbr[s] =   min(rxmbr[s],   sml[l].mbr[s]);
		            rxmbr[s+1] = max(rxmbr[s+1], sml[l].mbr[s+1]);
                }
            }
	    	marg += margin(dimension, rxmbr);

            // now calculate margin of R2
		    // initialize mbr of R2
		    for (s = 0; s < 2*dimension; s += 2) {
				rxmbr[s] =    MAXREAL;
				rxmbr[s+1] = -MAXREAL;
			}
		
            for ( ; l < n; l++) {
                // calculate mbr of R1
				for (s = 0; s < 2*dimension; s += 2) {
	            	rxmbr[s] =   min(rxmbr[s],   sml[l].mbr[s]);
	            	rxmbr[s+1] = max(rxmbr[s+1], sml[l].mbr[s+1]);
                }
            }
	    	marg += margin(dimension, rxmbr);
        }

        // for all possible distributions of smu
       	for (k = 0; k < n - 2*m1 + 1; k++) {
            // now calculate margin of R1
		    // initialize mbr of R1
		    for (s = 0; s < 2*dimension; s += 2) {
				rxmbr[s] =    MAXREAL;
				rxmbr[s+1] = -MAXREAL;
		    }
		    
            for (l = 0; l < m1+k; l++) {
                // calculate mbr of R1
				for (s = 0; s < 2*dimension; s += 2) {
	            	rxmbr[s] =   min(rxmbr[s],   smu[l].mbr[s]);
	            	rxmbr[s+1] = max(rxmbr[s+1], smu[l].mbr[s+1]);
                }
            }
	    	marg += margin(dimension, rxmbr);

            // now calculate margin of R2
		    // initialize mbr of R2
		    for (s = 0; s < 2*dimension; s += 2) {
				rxmbr[s] =    MAXREAL;
				rxmbr[s+1] = -MAXREAL;
		    }
            for ( ; l < n; l++) {
                // calculate mbr of R1
				for (s = 0; s < 2*dimension; s += 2) {
		            rxmbr[s] =   min(rxmbr[s],   smu[l].mbr[s]);
		            rxmbr[s+1] = max(rxmbr[s+1], smu[l].mbr[s+1]);
                }
            }
	    	marg += margin(dimension, rxmbr);
        }

        // actual margin better than optimum?
        if (marg < minmarg) {
            split_axis = i;
            minmarg = marg;
        }
    }

    printf("split axis = %d\n", split_axis);
    
    // choose best distribution for split axis
    for (j = 0; j < n; j++) {
		sml[j].index = smu[j].index = j;
		sml[j].dimension = smu[j].dimension = split_axis;
		sml[j].mbr = smu[j].mbr = mbr[j];
    }

    // Sort by lower and upper value perpendicular split axis
    qsort(sml, n, sizeof(SortMbr), sort_lower_mbr);
    qsort(smu, n, sizeof(SortMbr), sort_upper_mbr);

    minover = MAXREAL;
    mindead = MAXREAL;
    // for all possible distributions of sml and snu
    for (k = 0; k < n - 2*m1 + 1; k++) {
        // lower sort
		// now calculate margin of R1
		// initialize mbr of R1
        dead = 0.0;
		for (s = 0; s < 2*dimension; s += 2) {
		    rxmbr[s] =    MAXREAL;
		    rxmbr[s+1] = -MAXREAL;
		}
		for (l = 0; l < m1+k; l++) {
		    // calculate mbr of R1
		    for (s = 0; s < 2*dimension; s += 2) {
				rxmbr[s] =   min(rxmbr[s],   sml[l].mbr[s]);
				rxmbr[s+1] = max(rxmbr[s+1], sml[l].mbr[s+1]);
		    }
			dead -= area(dimension, sml[l].mbr);
		}
		dead += area(dimension, rxmbr);

	// now calculate margin of R2
	// initialize mbr of R2
	for (s = 0; s < 2*dimension; s += 2) {
	    rymbr[s] =    MAXREAL;
       	    rymbr[s+1] = -MAXREAL;
	}
	for ( ; l < n; l++) {
	    // calculate mbr of R1
	    for (s = 0; s < 2*dimension; s += 2) {
			rymbr[s] =   min(rymbr[s],   sml[l].mbr[s]);
			rymbr[s+1] = max(rymbr[s+1], sml[l].mbr[s+1]);
	    }
		dead -= area(dimension, sml[l].mbr);
	}
        dead += area(dimension, rymbr);
	over = overlap(dimension, rxmbr, rymbr);

        if ((over < minover) || (over == minover) && dead < mindead) {
            minover = over;
            mindead = dead;
            dist = m1+k;
            lu = TRUE;
        }

        // upper sort
	// now calculate margin of R1
	// initialize mbr of R1
        dead = 0.0;
	for (s = 0; s < 2*dimension; s += 2) {
	    rxmbr[s] =    MAXREAL;
	    rxmbr[s+1] = -MAXREAL;
	}
	for (l = 0; l < m1+k; l++) {
	    // calculate mbr of R1
	    for (s = 0; s < 2*dimension; s += 2) {
			rxmbr[s] =   min(rxmbr[s],   smu[l].mbr[s]);
			rxmbr[s+1] = max(rxmbr[s+1], smu[l].mbr[s+1]);
	    }
		dead -= area(dimension, smu[l].mbr);
	}
        dead += area(dimension, rxmbr);

	// now calculate margin of R2
	// initialize mbr of R2
	for (s = 0; s < 2*dimension; s += 2) {
	    rymbr[s] =    MAXREAL;
	    rymbr[s+1] = -MAXREAL;
	}
	for ( ; l < n; l++) {
	    // calculate mbr of R1
	    for (s = 0; s < 2*dimension; s += 2) {
			rymbr[s] =   min(rymbr[s],   smu[l].mbr[s]);
			rymbr[s+1] = max(rymbr[s+1], smu[l].mbr[s+1]);
	    }
            dead -= area(dimension, smu[l].mbr);
	}
        dead += area(dimension, rxmbr);
	over = overlap(dimension, rxmbr, rymbr);

        if ((over < minover) || (over == minover) && dead < mindead) {
            minover = over;
            mindead = dead;
            dist = m1+k;
            lu = FALSE;
        }
    }

    // calculate best distribution
    *distribution = new int[n];
    for (i = 0; i < n; i++) {
        if (lu)
            (*distribution)[i] = sml[i].index;
        else
            (*distribution)[i] = smu[i].index;
    }

    delete sml;
    delete smu;
    delete rxmbr;
    delete rymbr;


    return dist;
}


//////////////////////////////////////////////////////////////////////////////
// RTDirNode
//////////////////////////////////////////////////////////////////////////////


RTDirNode ::RTDirNode(RTree  *rt)
    : RTNode (rt)
{
    int i;
    char *b;
    int header_size;
    DirEntry  * d;

    d = new DirEntry (dimension, rt);
    header_size = sizeof(char) + sizeof(char) + sizeof(int);
    capacity = (rt->file->get_blocklength() - header_size) / d->get_size();
    delete d;

    entries = new DirEntry [capacity];
    for (i=0; i<capacity;i++)
    {
        entries[i].dimension = dimension;
        entries[i].my_tree = my_tree;
        entries[i].bounces = new float[2*dimension];
    }

    b = new char[rt->file->get_blocklength()];
    block = rt->file->append_block(b);
    delete [] b;
    rt->num_of_inodes ++;
    dirty = TRUE;
}

RTDirNode ::RTDirNode(RTree  *rt, int _block)
    : RTNode (rt)
{
    int i;
    char *b;
    int header_size;
    DirEntry  * d;

    d = new DirEntry (dimension, rt);
    header_size = sizeof(char) + sizeof(char) + sizeof(int);
    capacity = (rt->file->get_blocklength() - header_size) / d->get_size();
    delete d;

    entries = new DirEntry [capacity];
    for (i=0; i<capacity;i++)
    {
        entries[i].dimension = dimension;
        entries[i].my_tree = my_tree;
        entries[i].bounces = new float[2*dimension];
    }

    block = _block;
    b = new char[rt->file->get_blocklength()];
    rt->file->read_block(b, block);
    read_from_buffer(b);
    delete [] b;

    dirty = FALSE;
}

RTDirNode ::~RTDirNode() {
    char *b;
    if (dirty) {
		b = new char[my_tree->file->get_blocklength()];
		write_to_buffer(b);
		my_tree->file->write_block(b, block);
        delete [] b;
    }
    delete [] entries;
}

void RTDirNode ::read_from_buffer(char *buffer) {
    int i, j, s;

    memcpy(&son_is_data, buffer, sizeof(bool));
    j = sizeof(bool);

    memcpy(&level, &buffer[j], sizeof(char));
    j += sizeof(char);

    memcpy(&num_entries, &buffer[j], sizeof(int));
    j += sizeof(int);

    s = entries[0].get_size();
    for (i = 0; i < num_entries; i++)
    {
	entries[i].read_from_buffer(&buffer[j]);
	entries[i].son_is_data = son_is_data;
	j += s;
    }
}

void RTDirNode ::write_to_buffer(char *buffer) {
    int i, j, s;
        
    memcpy(buffer, &son_is_data, sizeof(bool));
    j = sizeof(bool);

    memcpy(&buffer[j], &level, sizeof(char));
    j += sizeof(char);

    memcpy(&buffer[j], &num_entries, sizeof(int));
    j += sizeof(int);

    s = entries[0].get_size();
    for (i = 0; i < num_entries; i++) {
		entries[i].write_to_buffer(&buffer[j]);
		j += s;
    }
}

void RTDirNode ::print()
{
    int i, n;

    n = get_num();
    for (i = 0; i < n ; i++)
    {
        printf("(%4.1lf, %4.1lf, %4.1lf, %4.1lf)\n",
	       entries[i].bounces[0],
	       entries[i].bounces[1],
	       entries[i].bounces[2],
	       entries[i].bounces[3]);
    }
    printf("level %d\n", level);
}

float * RTDirNode ::get_mbr() {
    int i, j, n;
    float *mbr;

    mbr = new float[2*dimension];
    for (i = 0; i < 2*dimension; i ++ )

	mbr[i] = entries[0].bounces[i];

    n = get_num();
    for (j = 1; j < n; j++)
    {
	for (i = 0; i < 2*dimension; i += 2)
	{
	    mbr[i]   = min(mbr[i],   entries[j].bounces[i]);
	    mbr[i+1] = max(mbr[i+1], entries[j].bounces[i+1]);
        }
    }
    return mbr;
}

void RTDirNode ::enter(DirEntry  *de) {
    // ist ein Einfuegen ueberhaupt moeglich?
    if (get_num() > (capacity-1))
        error("RTDirNode::enter: called, but node is full", TRUE);

    // Eintrag an erste freie Stelle kopieren
    entries[num_entries] = *de;

    // jetzt gibts einen mehr
    num_entries++;

    // Eintrag freigeben, aber keine Nachfolgeknoten loeschen !!!!
    de->son_ptr = NULL;
    delete de;
}

void RTDirNode ::split(RTDirNode  *sn) {
    int i, *distribution, dist, n;
    float **mbr_array;
    DirEntry  *new_entries1, *new_entries2;

#ifdef SHOWMBR
	split_000++;
#endif

    n = get_num();

    // mbr_array allokieren und belegen
    mbr_array = new floatptr[n];
    for (i = 0; i < n; i++)
       	mbr_array[i] = entries[i].bounces;

    // Verteilung berechnen
    dist = RTNode ::split(mbr_array, &distribution);

    // mbr_array freigeben
    delete [] mbr_array;

    // neues Datenarray erzeugen
    // --> siehe Konstruktor
    new_entries1 = new DirEntry [capacity];
    for (i=0; i<capacity;i++)
    {
        new_entries1[i].dimension = dimension;
        new_entries1[i].my_tree = my_tree;
        new_entries1[i].bounces = new float[2*dimension];
    }
    new_entries2 = new DirEntry [capacity];
    for (i=0; i<capacity;i++)
    {
        new_entries2[i].dimension = dimension;
        new_entries2[i].my_tree = my_tree;
        new_entries2[i].bounces = new float[2*dimension];
    }

    for (i = 0; i < dist; i++)
    {
       	new_entries1[i] = entries[distribution[i]];
    }

    for (i = dist; i < n; i++)
    {
       	new_entries2[i-dist] = entries[distribution[i]];
    }

    for (i = 0; i < n; i++)
    {
       	entries[i].son_ptr = NULL;
       	sn->entries[i].son_ptr = NULL;
    }
    delete [] entries;
    delete [] sn->entries;

    entries = new_entries1;
    sn->entries = new_entries2;

    num_entries = dist;
    sn->num_entries = n - dist;  // muss wegen Rundung so bleiben !!
}

int RTDirNode ::choose_subtree(float *mbr) {
    int i, j, n, follow, minindex, *inside, inside_count, *over;
    float *bmbr, old_o, o, omin, a, amin, f, fmin;

    n = get_num();

    // faellt d in einen bestehenden Eintrag ?
    inside_count = 0;
    inside = new int[n];
    over = new int[n];
    for (i = 0; i < n; i++)
    {
	switch (entries[i].section(mbr))
	{
	case INSIDE:
	    inside[inside_count++] = i;
	    break;
	}
    }

    if (inside_count == 1)
    // Fall 1: Rechteck faellt genau in einen Eintrag
	follow = inside[0];
    else if (inside_count > 1)
    // Fall 2: Rechteck faellt in mehrere Eintrag --> nimm den
    // mit der geringsten Flaeche (volumen!!!)
    {
	fmin = MAXREAL;
	for (i = 0; i < inside_count; i++)
	{
	    f = area(dimension, entries[inside[i]].bounces);
	    if (f < fmin)
	    {
		minindex = i;
      		fmin = f;
       	    }
       	}
	follow = inside[minindex];
    }
    else
    // Fall 3: Rechteck faellt in keinen Eintrag -->
    // fuer Knoten, die auf interne Knoten zeigen:
    // nimm den Eintrag, der am geringsten vergroessert wird;
    // bei gleicher Vergroesserung:
    // nimm den Eintrag, der die geringste Flaeche hat
    //
    // fuer Knoten, die auf Datenknoten zeigen:
    // nimm den, der die geringste Ueberlappung verursacht
    // bei gleicher Ueberlappung:
    // nimm den Eintrag, der am geringsten vergroessert wird;
    // bei gleicher Vergroesserung:
    // nimm den Eintrag, der die geringste Flaeche hat
    {
       	if (son_is_data)
	{
            omin = MAXREAL;
	    fmin = MAXREAL;
	    amin = MAXREAL;
	    for (i = 0; i < n; i++)
	    {
		// berechne die mbr, wenn mbr in entries[i] eingefuegt wird
		enlarge(dimension, &bmbr, mbr, entries[i].bounces);

		// calculate area and area enlargement
		a = area(dimension, entries[i].bounces);
		f = area(dimension, bmbr) - a;

		// calculate overlap before enlarging entry_i
		old_o = o = 0.0;
		for (j = 0; j < n; j++)
		{
		    if (j != i)
		    {
			old_o += overlap(dimension,
					 entries[i].bounces,
					 entries[j].bounces);
			o += overlap(dimension,
				     bmbr,
				     entries[j].bounces);
		    }
	        }
	        o -= old_o;

	        // is this entry better than the former optimum ?
	        if ((o < omin) ||
		    (o == omin && f < fmin) ||
		    (o == omin && f == fmin && a < amin))
	        {
	       	    minindex = i;
		    omin = o;
		    fmin = f;
		    amin = a;
	        }
	        delete [] bmbr;
	    }
        }
        else
        {
	    fmin = MAXREAL;
	    amin = MAXREAL;
	    for (i = 0; i < n; i++)
	    {
	        // berechne die mbr, wenn mbr in entries[i] eingefuegt wird
	        enlarge(dimension, &bmbr, mbr, entries[i].bounces);

	        // calculate area and area enlargement
	        a = area(dimension, entries[i].bounces);
	        f = area(dimension, bmbr) - a;

	        // is this entry better than the former optimum ?
	        if ((f < fmin) ||
		    (f == fmin && a < amin))
	        {
	       	    minindex = i;
		    fmin = f;
	            amin = a;
	        }
	        delete [] bmbr;
	    }
        }

	enlarge(dimension, &bmbr, mbr, entries[minindex].bounces);
	memcpy(entries[minindex].bounces, bmbr,
	       dimension * 2 * sizeof(float));
	follow = minindex;
	delete [] bmbr;

	dirty = TRUE;
    }

    delete [] inside;
    delete [] over;

    return follow;
}

R_OVERFLOW RTDirNode ::insert(DATA *d, RTNode  **sn) {
    int follow;
    RTNode  *succ, *new_succ;
    RTDirNode  *brother;
    DirEntry  *de;
    R_OVERFLOW ret;
    float *mbr,*nmbr;

    // Einfuegepfad waehlen
    mbr = d->get_mbr();
    follow = choose_subtree(mbr);
    delete [] mbr;

    // Sohn laden
    succ = entries[follow].get_son();

    // insert d into son
    ret = succ->insert(d, &new_succ);
    if (ret != NONE)
    // if anything happend --> update bounces of entry "follow"
    {
        mbr = succ->get_mbr();
        memcpy(entries[follow].bounces, mbr, sizeof(float) * 2*dimension);
        delete [] mbr;
    }

    // recalculate # of succeeders in the tree
    if (ret == SPLIT) {
        if (get_num() == capacity)
     	    error("RTDataNode::insert: maximum capacity violation", TRUE);

         // create and insert new entry
        de = new DirEntry (dimension, my_tree);
		nmbr = new_succ->get_mbr();
        memcpy(de->bounces, nmbr, 2*dimension*sizeof(float));
	delete [] nmbr;
        de->son = new_succ->block;
        de->son_ptr = new_succ;
        de->son_is_data = son_is_data;
        enter(de);

        if (get_num() == (capacity - 1)) {
	    brother = new RTDirNode (my_tree);
	    brother->son_is_data = son_is_data;
	    brother->level = level;
	    split(brother);
            *sn = brother;
	    printf("ddddddddiiiiiiirr: splitting node %d, creating %d\n",
		   block, (*sn)->block);

            ret = SPLIT;
	}
        else
      	    ret = NONE;
    }
    // must write page
    dirty = TRUE;

    return ret;
}

void RTDirNode ::point_query(float *p) {
    int i, n;
    RTNode  *succ;

    n = get_num();
    for (i = 0; i < n; i++) {
		if (entries[i].is_inside(p)) {
            succ = entries[i].get_son();
            succ->point_query(p);
        }
    }
}

void RTDirNode ::rangeQuery(float *mbr) {
    int i, n;
    SECTION s;
    RTNode  *succ;

#ifdef ZAEHLER
    page_access++;
#endif

    n = get_num();
    for (i = 0; i < n; i++) {
        s = entries[i].section(mbr);
        if (s == INSIDE || s == OVERLAP) {
            succ = entries[i].get_son();
            succ->rangeQuery(mbr);
        }
    }
}


//////////////////////////////////////////////////////////////////////////////
// RTDataNode
//////////////////////////////////////////////////////////////////////////////
RTDataNode ::RTDataNode(RTree  *rt)
    : RTNode (rt)
{
    char *b;
    int i;
    int header_size;
    DATA * d;

    level = 0;

    d = new DATA;
    header_size = sizeof(char) + sizeof(int);
    capacity = (rt->file->get_blocklength() - header_size) / d->get_size();
    delete d;

    data = new DATA[capacity];		// serious bug !!!
    //for (i=0; i<capacity;i++) data[i].data = new float[DIMENSION];

    b = new char[rt->file->get_blocklength()];
    block = rt->file->append_block(b);
    delete [] b;
    rt->num_of_dnodes ++;
    dirty = TRUE;
}

RTDataNode ::RTDataNode(RTree  *rt, int _block)
    : RTNode (rt)
{
    int i;
    char *b;
    int header_size;
    DATA *d;

    d = new DATA;
    header_size = sizeof(char) + sizeof(int);
    capacity = (rt->file->get_blocklength() - header_size) / d->get_size();
    delete d;

    data = new DATA[capacity];		// serious bug !!!
    //for (i=0; i<capacity;i++) data[i].data = new float[dimension];

    block = _block;
    b = new char[rt->file->get_blocklength()];
    rt->file->read_block(b, block);
    read_from_buffer(b);
    delete [] b;

    dirty = FALSE;
}

RTDataNode ::~RTDataNode()
{
    char *b;

    if (dirty)
    {
	b = new char[my_tree->file->get_blocklength()];
	write_to_buffer(b);
	my_tree->file->write_block(b, block);
        delete [] b;
    }

    delete [] data;
}

void RTDataNode ::read_from_buffer(char *buffer) {
    int i, j, s;

    memcpy(&level, buffer, sizeof(char));
    j = sizeof(char);

    level = 0;
    
    memcpy(&num_entries, &buffer[j], sizeof(int));
    j += sizeof(int);

    s = data[0].get_size();
    for (i = 0; i < num_entries; i++)
    {
	data[i].read_from_buffer(&buffer[j]);
	j += s;
    }
}

void RTDataNode ::write_to_buffer(char *buffer) {
    int i, j, s;

    level = 0;
    memcpy(&level, buffer, sizeof(char));
    j = sizeof(char);

    memcpy(&buffer[j], &num_entries, sizeof(int));
    j += sizeof(int);

    s = data[0].get_size();
    for (i = 0; i < num_entries; i++)
    {
	data[i].write_to_buffer(&buffer[j]);
	j += s;
    }
}


void RTDataNode ::print() {
    int i, n;

    n = get_num();
    for (i = 0; i < n ; i++)
    {
        printf("(%4.1lf, %4.1lf)\n",
	       data[i].data[0],
	       data[i].data[1]);
    }
    printf("level %d\n", level);
}

float * RTDataNode ::get_mbr() {
    int i, j, n;
    float *mbr, *tm;

    mbr = data[0].get_mbr();
    n = get_num();
    for (j = 1; j < n; j++)
    {
	tm = data[j].get_mbr();
	for (i = 0; i < 2*dimension; i += 2)
	{
		mbr[i]   = min(mbr[i],   tm[i]);
		mbr[i+1] = max(mbr[i+1], tm[i+1]);
       	}
       	delete [] tm;
    }
    return mbr;
}

void RTDataNode ::split(RTDataNode  *sn) {
    int i, *distribution, dist, n;
    float **mbr_array;
    DATA *new_data1, *new_data2;

#ifdef SHOWMBR
	split_000++;
#endif

    n = get_num();

    // mbr_array allokieren und belegen
    mbr_array = new floatptr[n];
    for (i = 0; i < n; i++) {
        mbr_array[i] = data[i].get_mbr();
    }

    dist = RTNode ::split(mbr_array, &distribution);

    new_data1 = new DATA[capacity];
    //for (i=0; i<capacity;i++)	// serious bug !
    //    new_data1[i].data = new float[2*dimension];
    
    new_data2 = new DATA[capacity];
    //for (i=0; i<capacity;i++)	// serious bug !
    //    new_data2[i].data = new float[2*dimension];

    for (i = 0; i < dist; i++)
        new_data1[i] = data[distribution[i]];

    for (i = dist; i < n; i++)
        new_data2[i-dist] = data[distribution[i]];

    delete [] data;
	delete [] ((RTDataNode *)sn)->data;

    data = new_data1;
    sn->data = new_data2;

    num_entries = dist;
    sn->num_entries = n - dist;

    for (i = 0; i < n; i++)
	delete [] mbr_array[i];

    delete [] mbr_array;
}

R_OVERFLOW RTDataNode ::insert(DATA *d, RTNode  **sn) {
    int i, last_cand;
    float *mbr, *center;
    SortMbr *sm;
    DATA *nd, *new_data;

    if (get_num() == capacity)
	error("RTDataNode::insert: maximum capacity violation", TRUE);

    data[get_num()] = *d;
    num_entries ++;

    // Plattenblock zum Schreiben markieren
    dirty = TRUE;

    if (get_num() == (capacity - 1)) {
        if (my_tree->re_level[0] == FALSE)
	// there was no reinsert on level 0 during this insertion
        {
            // calculate center of page
            mbr = get_mbr();
            center = new float[2*dimension];
            for (i = 0; i < dimension; i++)
                 center[i] = (mbr[2*i] + mbr[2*i+1]) / 2.0;

            new_data = new DATA[capacity];
            //for (i=0; i<capacity;i++)
            //    new_data[i].data = new float[dimension];

		    sm = new SortMbr[num_entries];
		    for (i = 0; i < num_entries; i++)
		    {
				sm[i].index = i;
				sm[i].dimension = dimension;
				sm[i].mbr = data[i].get_mbr();
				sm[i].center = center;
			}

            // sort by distance of each center to the overall center
            qsort(sm, num_entries, sizeof(SortMbr), sort_center_mbr);

            last_cand = (int) ((float)num_entries * 0.30);

            // copy the nearest candidates to new array
            for (i = 0; i < num_entries - last_cand; i++)
	        new_data[i] = data[sm[i].index];

            // insert candidates into reinsertion list
            for ( ; i < num_entries; i++)
            {
                nd = new DATA;
                *nd = data[sm[i].index];
                my_tree->re_data_cands.push_front(nd);
            }

            // free and copy data array
            delete [] data;
	    data = new_data;

       	    for (i = 0; i < num_entries; i++)
	       	delete [] sm[i].mbr;
	    delete sm;
	    delete [] mbr;
	    delete [] center;
	    my_tree->re_level[0] = TRUE;

	    // correct # of entries
	    num_entries -= last_cand;

	    // must write page
	    dirty = TRUE;

		return REINSERT;
	} else {
	    *sn = new RTDataNode (my_tree);
	    (*sn)->level = level;
	    split((RTDataNode *) *sn);
	    printf("splitting node %d, creating %d\n", block, (*sn)->block);
	}
	return SPLIT;
    }
    else
        return NONE;
}

void RTDataNode ::point_query(float *p)
{
    int i, n,j;
    float *dmbr;

#ifdef ZAEHLER
    page_access++;
#endif

    n = get_num();
    for (i = 0; i < n; i++) {
        dmbr = data[i].get_mbr();
        if (section(dimension, p, dmbr))
	{
/*		printf("( ");
		for( j = 0; j < dimension; j++)
			printf(" %f",dmbr[j]);
		printf(" ) \n");
*/	}
	delete [] dmbr;
    }
}

void RTDataNode ::rangeQuery(float *mbr) {
    int i, n,j,k;
    float *dmbr;
    DATA *point, *center;

#ifdef ZAEHLER
    page_access++;
#endif

    n = get_num();
    for (i = 0; i < n; i++) {
        dmbr = data[i].get_mbr();
        if (section(dimension, dmbr, mbr))	{
		    point = new DATA;
		    *point = data[i];
			//res->insert(point);
		}
        delete [] dmbr;
    }
}

// RTree

RTree ::RTree(char *fname, int _b_length, int cache_size, int _dimension) {
    file = new CachedBlockFile(fname, _b_length, cache_size);

    // header allokieren und lesen
    header = new char [file->get_blocklength()];

    read_header(header);

    dimension = _dimension;
    root = 0;
    root_ptr = NULL;
    root_is_data = TRUE;
    num_of_data = num_of_inodes = num_of_dnodes = 0;

    root_ptr = new RTDataNode (this);
    root = root_ptr->block;
}


RTree ::RTree(char *fname, int cache_size) {
    int j;

    file = new CachedBlockFile(fname, 0, cache_size);

    header = new char [file->get_blocklength()];
    file->read_header(header);
    read_header(header);
    root_ptr = NULL;

    if (get_num() == 0) {
		root_ptr = new RTDataNode (this);
		root = root_ptr->block;
		root_ptr->level = 0;
	} else
		load_root();
}

RTree ::~RTree() {
    int j;

    write_header(header);

    file->set_header(header);
    delete [] header;

    if (root_ptr != NULL) {
        delete root_ptr;
        root_ptr = NULL;
    }
    delete file;

	if (!(file->isReadOnly))
    	printf("saved R-Tree containing %d internal, %d data nodes and %d data\n",
	   		num_of_inodes, num_of_dnodes, num_of_data);
}

void RTree ::read_header(char *buffer) {
    int i;

    memcpy(&dimension, buffer, sizeof(dimension));
    i = sizeof(dimension);

    memcpy(&num_of_data, &buffer[i], sizeof(num_of_data));
    i += sizeof(num_of_data);

    memcpy(&num_of_dnodes, &buffer[i], sizeof(num_of_dnodes));
    i += sizeof(num_of_dnodes);

    memcpy(&num_of_inodes, &buffer[i], sizeof(num_of_inodes));
    i += sizeof(num_of_inodes);

    memcpy(&root_is_data, &buffer[i], sizeof(root_is_data));
    i += sizeof(root_is_data);

    memcpy(&root, &buffer[i], sizeof(root));
    i += sizeof(root);

    user_header = &buffer[i];
}

void RTree ::write_header(char *buffer)
{
    int i;

    memcpy(buffer, &dimension, sizeof(dimension));
    i = sizeof(dimension);

    memcpy(&buffer[i], &num_of_data, sizeof(num_of_data));
    i += sizeof(num_of_data);

    memcpy(&buffer[i], &num_of_dnodes, sizeof(num_of_dnodes));
    i += sizeof(num_of_dnodes);

    memcpy(&buffer[i], &num_of_inodes, sizeof(num_of_inodes));
    i += sizeof(num_of_inodes);

    memcpy(&buffer[i], &root_is_data, sizeof(root_is_data));
    i += sizeof(root_is_data);

    memcpy(&buffer[i], &root, sizeof(root));
    i += sizeof(root);

    user_header = &buffer[i];
}

void RTree ::load_root() {
    if (root_ptr == NULL) {
		if (root_is_data)
	    	root_ptr = new RTDataNode (this, root);
        else
	    	root_ptr = new RTDirNode (this, root);
    }
}

void RTree ::insert(DATA* d) {
    int i, j;
    RTNode  *sn;
    RTDirNode  *nroot_ptr;
    int nroot;
    DirEntry  *de;
    R_OVERFLOW split_root;
    DATA *d_cand, *dc;
    float *nmbr;

    // load root into memory
    load_root();

    // no overflow occured until now
    re_level = new bool[root_ptr->level+1];
    for (i = 0; i <= root_ptr->level; i++)
        re_level[i] = FALSE;
        
    // insert d into re_data_cands as the first entry to insert
    // make a copy of d because it shouldnt be erased later
    dc = new DATA;
    *dc = *d;
    re_data_cands.push_front(dc);

    j = -1;
    while (re_data_cands.size()>0) {
        // first try to insert data, then directory entries
		d_cand = (DATA*)re_data_cands.front();
        if (d_cand != NULL) {
            // since erase deletes the according data element of the
            // list, we should make a copy of the data before erasing it
            dc = new DATA;
            *dc = *d_cand;
            re_data_cands.pop_front();

	    	// start rekursive insert with root
		    split_root = root_ptr->insert(dc, &sn);
		    delete dc;
        } else
	    	error("RTree::insert: inconsistent list re_data_cands", TRUE);

		if (split_root == SPLIT) {
			// insert has lead to split --> new root-page having to sons
			// root and sn
		    nroot_ptr = new RTDirNode (this);
		    nroot_ptr->son_is_data = root_is_data;
		    nroot_ptr->level = root_ptr->level + 1;
		    nroot = nroot_ptr->block;
	
		    de = new DirEntry (dimension, this);
		    nmbr = root_ptr->get_mbr();
		    memcpy(de->bounces, nmbr, 2*dimension*sizeof(float));
		    delete [] nmbr;
		    de->son = root_ptr->block;
		    de->son_ptr = root_ptr;
		    de->son_is_data = root_is_data;
		    nroot_ptr->enter(de);
	
		    de = new DirEntry (dimension, this);
		    nmbr = sn->get_mbr();
		    memcpy(de->bounces, nmbr, 2*dimension*sizeof(float));
		    delete [] nmbr;
		    de->son = sn->block;
		    de->son_ptr = sn;
		    de->son_is_data = root_is_data;
		    nroot_ptr->enter(de);
	
		    root = nroot;
			root_ptr = nroot_ptr;
	        root_is_data = FALSE;
        }
        j++;
    }
    num_of_data++;
    
    // free lists
    delete [] re_level;
}

void RTree ::point_query(float *p) {
    // load root node into main memory
    load_root();
    root_ptr->point_query(p);
}

void RTree ::rangeQuery(float *mbr) {
    load_root();
    root_ptr->rangeQuery(mbr);
}

