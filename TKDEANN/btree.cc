#include "BTree.h"

// BTDirNode

void BTDirNode::read_entry_from_buffer(char *buffer,int pos) {
	memcpy(&keys[pos], buffer, sizeof(int));
	memcpy(&sons[pos], buffer+sizeof(int), sizeof(int));
}

void BTDirNode::write_entry_to_buffer(char *buffer,int pos) {
    memcpy(buffer, &keys[pos], sizeof(int));
    memcpy(buffer+sizeof(int), &sons[pos], sizeof(int));
}

// creates a brand new BT directory node
BTDirNode::BTDirNode(CachedBlockFile *cf) { 
	// from parent class
	my_file=cf;		num_entries = 0;	block = -1;
	
	// header of page keeps node info 
	// level + num_entries
	int header_size = sizeof(char) + sizeof(int);
	capacity = (my_file->get_blocklength() - header_size) / (EntrySize);
	
	// Initialize entries
	keys=new int[capacity];		sons=new int[capacity];
	
	// create new block for the node
	char* b = new char[my_file->get_blocklength()];
	block = my_file->append_block(b);
	delete [] b;
	
	dirty = true;	// must be written to disk before destruction
}

// reads an existing BT directory node
BTDirNode::BTDirNode(CachedBlockFile *cf, int _block) {	
	// from parent class
	my_file=cf;		num_entries = 0;	block = -1;
	
	// header of page keeps node info 
	// level + num_entries
	int header_size = sizeof(char) + sizeof(int);
	capacity = (my_file->get_blocklength() - header_size) / (EntrySize);
	
	// Initialize entries
	keys=new int[capacity];		sons=new int[capacity];
	
	// now load block and read BTDirNode data from it
	block = _block;
	char* b = new char[my_file->get_blocklength()];
	my_file->read_block(b, block);    
	read_from_buffer(b);
	delete [] b;
	
	dirty = FALSE;		// not dirty yet
}

BTDirNode::~BTDirNode() {
    if (dirty){
		// Update changes on disk
		char* b = new char[my_file->get_blocklength()];
		write_to_buffer(b);
		my_file->write_block(b, block);
        delete [] b;
    }
    
    // destroy entries
    if (keys!=NULL) delete[] keys;
    if (sons!=NULL) delete[] sons;
}

void BTDirNode::read_from_buffer(char *buffer) {
    int j=0;

    // first read header info
    memcpy(&level, &buffer[j], sizeof(char));
    j += sizeof(char);
    
    memcpy(&num_entries, &buffer[j], sizeof(int));
    j += sizeof(int);

    // then read entries
    for (int i = 0; i < num_entries; i++) {
		read_entry_from_buffer(&buffer[j],i);
		j += EntrySize;
    }
}

void BTDirNode::write_to_buffer(char *buffer) {
    int j=0;

    // first, write header info
    memcpy(&buffer[j], &level, sizeof(char));
    j += sizeof(char);
    memcpy(&buffer[j], &num_entries, sizeof(int));
    j += sizeof(int);

    // then, write entries
    for (int i = 0; i < num_entries; i++) {
		write_entry_to_buffer(&buffer[j],i);
		j += EntrySize;
    }
}

// BTree

// Construction of a new BTree
BTree::BTree(char *fname, int _b_length, int cache_size) {	
    file = new CachedBlockFile(fname, _b_length, cache_size);
    
    // first block is header
    header = new char [file->get_blocklength()];

    root = 0;
    root_ptr = NULL;
    num_of_inodes = 0;
    
    //create an (empty) root
    root_ptr = new BTDirNode (this->file);
    this->num_of_inodes++;		// 27/1/2004 changed
    root = root_ptr->block;
}

// load an existing BTree
BTree::BTree(char *fname, int cache_size) { 
    file = new CachedBlockFile(fname, 0, cache_size);	// load file

    // read header
    header = new char [file->get_blocklength()];
    file->read_header(header);
    read_header(header);
    root_ptr = NULL;
}

BTree::~BTree() {
    write_header(header);
    file->set_header(header);
    delete [] header;

    if (root_ptr!=NULL) {
        delete root_ptr;
        root_ptr=NULL;
    }
    delete file;
    //printf("saved Tree containing %d internal nodes\n",num_of_inodes);
}

void BTree::read_header(char *buffer) {
    int i = 0;

    memcpy(&num_of_inodes, &buffer[i], sizeof(num_of_inodes));
    i += sizeof(num_of_inodes);

    memcpy(&root, &buffer[i], sizeof(root));
    i += sizeof(root);

	memcpy(&UserField, &buffer[i], sizeof(UserField));
    i += sizeof(UserField);

    user_header = &buffer[i];
}

void BTree::write_header(char *buffer) {
    int i = 0;

    memcpy(&buffer[i], &num_of_inodes, sizeof(num_of_inodes));
    i += sizeof(num_of_inodes);

    memcpy(&buffer[i], &root, sizeof(root));
    i += sizeof(root);

	memcpy(&buffer[i], &UserField, sizeof(UserField));
    i += sizeof(UserField);

    user_header = &buffer[i];
}

void BTree::load_root() {
    if (root_ptr==NULL)
		root_ptr=new BTDirNode (this->file, root);
}

