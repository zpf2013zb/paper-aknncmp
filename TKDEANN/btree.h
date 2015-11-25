#ifndef __BTree
#define __BTree

#include "blk_file.h"

class BTDirNode {
public:
	int capacity;    	// max. # of entries
	int num_entries;  	// # of used entries
	bool dirty;        	// TRUE, if node has to be written
	int block;         	// disc block
	char level;			// level of the node in the BTree
	
    void read_from_buffer(char *buffer);	// reads data from buffer
    void write_to_buffer(char *buffer); 	// writes data to buffer
    
    BTDirNode(CachedBlockFile *cf);
    BTDirNode(CachedBlockFile *cf, int _block);
    virtual ~BTDirNode();
    
	static const int EntrySize=sizeof(int)+sizeof(int);		// size of an entry
private:
	CachedBlockFile *my_file;	// storage manager of the object creator
    void read_entry_from_buffer(char *buffer,int pos);	// reads an entry from buffer
    void write_entry_to_buffer(char *buffer,int pos);	// writes an entry to buffer
    
	// array of entries
    int* keys;
    int* sons;
};


class BTree {
public:
	int UserField;
    int root;				// block # of root node
    BTDirNode*	root_ptr;	// root-node
    int num_of_inodes;		// # of stored directory pages
    CachedBlockFile *file;	// storage manager for harddisc blocks
    void load_root();   	// loads root_node into memory
    char *header;
protected:
    char *user_header;
    void read_header(char *buffer);  	// reads RBTree header
    void write_header(char *buffer); 	// writes RBTree header
public:
    BTree(char *fname,int _b_length, int cache_size);
    BTree(char *fname, int cache_size);
    virtual ~BTree();
};

#endif // __BTree
