#ifndef __RTREE
#define __RTREE

#include "blk_file.h"
#include <math.h>
using namespace std;
#include <deque>

// defines
#define MAXREAL 9.99e20
#define MAX_DIMENSION 256
const int DIMENSION=2;

// C-prototypes
class DATA;

void error(char *t, bool ex);

// enums / types
enum SECTION {OVERLAP, INSIDE, S_NONE};
enum R_OVERFLOW {SPLIT, REINSERT, NONE};
typedef float *floatptr;

struct SortMbr {
    int dimension;
    float *mbr;
    float *center;
    int index;
};

// DATA class for rectangle data
class DATA
{
public:
    float *data;                       // Vector
    float *get_mbr();                  // returns MBR of the object
    float get_area()                   // returns the area of the MBR
    { return 0.0; }                    // for vectors always 0.0
    void read_from_buffer(char *buffer);
                                       // reads data from buffer
    void write_to_buffer(char *buffer);
                                       // writes data to buffer
    int get_size();                    // returns amount of needed space
    int id,id2;
    void print();

    virtual DATA & operator = (DATA &_d);

    DATA();
    virtual ~DATA();
};

class RTNode;

class DirEntry {
    friend class RTDirNode ;
    friend class RTree ;

public:
    RTree  *my_tree;               // pointer to my R-tree
    int son;                            // block # of son
    RTNode  *son_ptr;              // pointer to son if in main mem.
    bool son_is_data;                   // TRUE, if son is a data page
    int dimension;                      // dimension of the box
    float *bounces;                     // pointer to box
    int son_level;                      // level of the node pointed to
                                        // son of this entry
    bool is_inside(float *v);           // tests, if v is inside the box
    SECTION section(float *mbr);        // tests, if mbr intersects the box
    RTNode * get_son();            // returns pointer to the son,
                                        // loads son if necessary

    bool section_circle(DATA *center, float radius);

    void read_from_buffer(char *buffer);// reads data from buffer
    void write_to_buffer(char *buffer); // writes data to buffer

    int get_size();                     // returns amount of needed buffer space

    virtual DirEntry  & operator = (DirEntry  &_d);

    DirEntry(int dimension = 0,RTree  *rt = NULL);
    
    virtual ~DirEntry();
};

class RTNode {
    friend class RTree ;
public:
    RTree  *my_tree;               // pointer to R-tree
    int capacity;                       // max. # of entries
    int dimension;
    int num_entries;                    // # of used entries
    bool dirty;                         // TRUE, if node has to be written

    int block;                          // disc block
    char level;                         // level of the node in the tree

    int get_num()                       // returns # of used entries
    { return num_entries;}

    virtual void read_from_buffer(char *buffer) = 0;
                                        // reads data from buffer
    virtual void write_to_buffer(char *buffer) = 0;
                                        // writes data to buffer

    virtual bool is_data_node() = 0;	// returns TRUE, if "this" is
                                        // *RTDataNode

    virtual float *get_mbr() = 0;       // returns mbr enclosing whole page
    virtual void print()                // prints rectangles
    { }

    int split(float **mbr, int **distribution);
                                        // splits an Array of mbr's into 2
                                        // and returns in distribution
                                        // which *m mbr's are to be
                                        // moved
    virtual void rangeQuery(float *mbr) = 0;

    virtual R_OVERFLOW insert(DATA *d, RTNode  **sn) = 0;
                                        // inserts d recursivly, if there
                                        // occurs a split, FALSE will be
                                        // returned and in sn a
                                        // pointer to the new node

    virtual void point_query(float *p) = 0;  // prints all entries equal to p
    
    RTNode(RTree  *rt);
    RTNode(RTree  *rt, int _block);
    virtual ~RTNode();
};


class RTDirNode: public RTNode {
    friend class RTree ;
public:
    DirEntry  *entries;            // array of entries
    bool son_is_data;                   // TRUE, if son is a data page
    void read_from_buffer(char *buffer);// reads data from buffer
    void write_to_buffer(char *buffer); // writes data to buffer
    
    bool is_data_node() {return FALSE;}

    void print();                       // prints rectangles

    float *get_mbr();                   // returns mbr enclosing whole page

    void enter(DirEntry  *de);     // inserts new entry

    void split(RTDirNode  *sn);    // splits directory page to *this and sn

    int choose_subtree(float *brm);     // chooses best subtree for insertion
    
    void rangeQuery(float *mbr);

    R_OVERFLOW insert(DATA *d, RTNode  **sn);
                                        // inserts d recursivly, if there
                                        // occurs a split, FALSE will be
                                        // returned and in sn a
                                        // pointer to the new node
    void point_query(float *p);         // prints all entries equal to p

    RTDirNode(RTree  *rt);
    RTDirNode(RTree  *rt, int _block);
    virtual ~RTDirNode();
};


class RTDataNode: public RTNode {
public:
    DATA *data;                         // array of data
    void read_from_buffer(char *buffer);// reads data from buffer
    void write_to_buffer(char *buffer); // writes data to buffer

    bool is_data_node() {return TRUE;}
    
    void print();                       // prints vector

    float *get_mbr();                   // returns mbr enclosing whole page

    void split(RTDataNode  *sn);   // splits data page into sn *this
    
    void rangeQuery(float *mbr);

    R_OVERFLOW insert(DATA *d, RTNode  **sn);
                                        // inserts d recursivly, if there
                                        // occurs a split, FALSE will be
                                        // returned and in sn a
                                        // pointer to the new node

    void point_query(float *mbr);       // prints all entries equal to p

    RTDataNode(RTree  *rt);
    RTDataNode(RTree  *rt, int _block);
    virtual ~RTDataNode();
};

class RTree
{
public:
    friend class RTNode ;
    friend class RTDirNode ;
    friend class RTDataNode ;

    int root;                            // block # of root node
    RTNode  *root_ptr;              // root-node
    bool root_is_data;                   // TRUE, if root is a data page
    int dimension;                       // dimension of the data's

    int num_of_data;	                 // # of stored data
    int num_of_dnodes;	                 // # of stored data pages
    int num_of_inodes;	                 // # of stored directory pages

    bool *re_level;                      // if re_level[i] is TRUE,
                                         // there was a reinsert on level i
    deque<void*> re_data_cands;         // data entries to reinsert

    CachedBlockFile *file;	         // storage manager for harddisc blocks

    void load_root();                    // loads root_node into memory
    char *header;
protected:
    char *user_header;
    void read_header(char *buffer);      // reads Rtree header
    void write_header(char *buffer);     // writes Rtree header
public:
    int get_num()                        // returns # of stored data
    { return num_of_data; }

    void insert(DATA *d);                // inserts new data into tree

    void point_query(float *p);          // print all entries equal to p
    
    void rangeQuery(float *mbr);
    
    RTree(char *fname,int _b_length, int cache_size,int _dimension);
    RTree(char *fname, int cache_size);
    virtual ~RTree();
};

#endif // __RTREE


