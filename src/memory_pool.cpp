#include "memory_pool.h"


memory_pool::memory_pool(uint64_t pool_size)
{
	this->current_pool = 0;
	//this->max_pools = max_pools;
	this->mem_pool = new std::vector<char *>();
	//this->mem_pool = (char **) std::calloc(max_pools, sizeof(char *));
	//this->base = (uint64_t *) std::malloc(max_pools * sizeof(uint64_t));
	this->base = new std::vector<uint64_t>();
	//this->base[0] = 0;
	this->base->push_back(0);

	//if (this->mem_pool == NULL) terror("Could not allocate memory pools");
	//this->mem_pool[0] = (char *) std::calloc(pool_size, sizeof(char));
	char * a_pool = (char *) std::calloc(pool_size, sizeof(char));
	this->mem_pool->push_back(a_pool);
	//if (this->mem_pool[0] == NULL) terror("Could not allocate initial memory pool");
	this->max_pools = 1;
	this->pool_size = pool_size;
}

void * memory_pool::request_bytes(uint64_t n_bytes)
{
	void * ptr;
	if (this->base->at(this->current_pool) + n_bytes >= this->pool_size) {

		this->current_pool++;

		if(this->current_pool == this->max_pools){
			//if(this->current_pool == this->max_pools) terror("Reached maximum number of pools. Exiting.");
			//this->mem_pool[this->current_pool] = (char *) std::calloc(this->pool_size, sizeof(char));
			char * a_pool = (char *) std::calloc(this->pool_size, sizeof(char));
			if(a_pool == NULL) throw "Could not allocate new memory pool";
			this->mem_pool->push_back(a_pool);
			//if (this->mem_pool[this->current_pool] == NULL) terror("Could not allocate memory pool");
			this->base->push_back(0);
			this->max_pools++;
		}
		
	}
	
	ptr = &this->mem_pool->at(this->current_pool)[0] + this->base->at(this->current_pool);
	this->base->at(this->current_pool) = this->base->at(this->current_pool) + n_bytes;
	
	return ptr;
}

void memory_pool::reset_n_bytes(uint64_t bytes){
	//Makes no checks, assuming you allocated some a priori
	if(bytes >= this->base->at(current_pool)){
		this->base->at(current_pool) = this->base->at(current_pool) - bytes;
	}
	/*
	if(bytes >= this->base[current_pool]){
		this->base[current_pool] = this->base[current_pool] - bytes;
	}
	*/
}

void memory_pool::full_reset(){
	
	for (uint64_t i = 0; i <= this->current_pool; i++) {
		memset(this->mem_pool->at(i), 0, this->pool_size);
		/*
		for(uint64_t j=0;j<5;j++){
			if(&this->mem_pool->at(0)[0] == NULL) printf("Is good\n");
		}
		*/
		this->base->at(i) = 0;
	}
	this->current_pool = 0;
}

memory_pool::~memory_pool(){

	//For the whole list
	for(uint64_t i=0;i<=this->current_pool;i++){
		std::free(this->mem_pool->at(i));
	}
	delete this->mem_pool;
	delete this->base;

	/*
	for (uint64_t i = 0; i <= this->current_pool; i++) {
		free(this->mem_pool[i]);
	}
	free(this->mem_pool);
	free(this->base);
	*/
}