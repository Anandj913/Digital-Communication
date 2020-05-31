#include "header.h"

char state[total_state];
float **tranMat, **PretranMat;

struct code
{
	int length;
	char *seq;
};
struct heapnode
{
	int symbol;
	float prob;

	struct heapnode *l, *r;
};
struct Heap
{
	int size;
	int capacity;
	struct heapnode **H;
};
struct heapnode* newnode(int data, float prob)
{
	struct heapnode *A = (struct heapnode*)malloc(sizeof(struct heapnode));
	A->symbol = data;
	A->prob = prob;
	A->l = A->r = NULL;
	return A;
}
void swap(struct heapnode**A, struct heapnode**B)
{
	struct heapnode* t = *A;
	*A = *B;
	*B = t;
}
struct Heap * createHeap(int capacity)
{
	struct Heap *A = (struct Heap*)malloc(sizeof(struct Heap));
	A->capacity = capacity;
	A->size = 0;
	A->H = (struct heapnode**)malloc(capacity*sizeof(struct heapnode*));
	return A;
	
}
void heapify(struct Heap *heap, int idx)
{
	int s = idx;
	int l = 2*idx +1;
	int r = 2*idx + 2;
	if(l < heap->size && heap->H[l]->prob <heap->H[s]->prob)
		s = l;
	if(r < heap->size && heap->H[r]->prob <heap->H[s]->prob)
		s = r;
	if(s!=idx)
	{
		swap(&heap->H[s], &heap->H[idx]);
		heapify(heap, s);
	}
}
struct heapnode* findmin(struct Heap* H)
{
	struct heapnode* t = H->H[0];
	H->H[0] = H->H[H->size-1];
	H->size = H->size -1;
	heapify(H, 0);
	return t;
}
void insert(struct Heap *H, struct heapnode *node)
{
	int i = H->size;
	H->size = H->size+1;
	while(i && node->prob < H->H[(i-1)/2]->prob)
	{
		H->H[i] = H->H[(i-1)/2];
		i = (i-1)/2;
	}
	H->H[i] = node;
}
void makeheap(struct Heap *H)
{
	int n = H->size -1;
	int i;
	for(i = (n-1)/2; i>=0; --i)
		heapify(H, i);
}
struct Heap* createandbuildHeap(float *prob, int size)
{
	struct Heap *H = createHeap(size);
	for(int i=0; i<size; i++)
	{
		H->H[i] = newnode(i, prob[i]);
	}
	H->size = size;
	makeheap(H);
	return H;
}
struct heapnode* makeHuffTree(float *prob, int size)
{
	struct Heap *H = createandbuildHeap(prob, size);
	struct heapnode *l, *r, *t;
	while(H->size != 1)
	{
		l = findmin(H);
		r = findmin(H);
		t = newnode(-1, l->prob + r->prob);
		t->l = l;
		t->r = r;
		insert(H, t);
	}
	return findmin(H);
}
void makecode(struct code **A, struct heapnode *root, char arr[], int size)
{
	if(root->l != NULL)
	{
		arr[size] = '0';
		makecode(A, root->l, arr, size+1);
	}
	if(root->r != NULL)
	{
		arr[size] = '1';
		makecode(A, root->r, arr, size+1);
	}

	if(root->l == NULL && root->r == NULL)
	{
		A[root->symbol] = (struct code*)malloc(sizeof(struct code));
		A[root->symbol]->length = size;
		A[root->symbol]->seq = (char *)malloc(size*sizeof(char));
		for(int i=0; i<size; i++)
		{
			A[root->symbol]->seq[i] = arr[i];
		}
	}
}
struct code** assigncode(struct heapnode* root, int size)
{
	struct code **A = (struct code**)malloc(size*sizeof(struct code*));
	char arr[10];
	makecode(A, root, arr, 0);
	return A;
} 
// To make transition Probability matrix for Markov process with number of state given by total_state 
float ** maketranMat( int state)
{
	float **tMat = (float **)malloc(state*sizeof(float *));
	int i,j;
	float sum,k;
    for(i=0; i<state; i++)
	{
		tMat[i] = (float *)calloc(state, sizeof(float));
	}
	for(i=0; i<state; i++)
	{
		sum=1;
		for(j=0; j<state-1; j++)
		{
			if(uniform_transition_mat == 1)
				k = 1.0/total_state;
			else
				k = sum*((float)rand()/(float)(1.5*RAND_MAX)); //maximum sum*2/3

			tMat[i][j] = k;
			sum = sum - k;
		}
		tMat[i][state-1] = sum; 
	}
	return tMat;
}

// To find commulative sum of all state in TPM
float ** makepretranMat(float **tranMat, int state)
{
	int i,j,k;
	float sum;
	float **PMat = (float **)malloc(state*sizeof(float *));
	for(i=0; i<state; i++)
	{
		PMat[i] = (float *)calloc(state, sizeof(float));
	}
	for(i=0;i<state; i++)
	{
		PMat[i][0] = tranMat[i][0];
		for(j=1; j<state; j++)
		{
			PMat[i][j] = PMat[i][j-1] + tranMat[i][j];
		}
	}
	return PMat;
}
// Return index of next state depending upon current state and TMP
int next(float **PretranMat, int current_state)
{
	float k = (float)rand()/(float)RAND_MAX;
	int j=0;
	while(k > PretranMat[current_state][j])
	{
		j++;
	}
	return j;
}
// Return index of 1st state 
int first()
{
	float *A = (float *)malloc(sizeof(float)*total_state);
	int i,j=0;
	for(i=0; i<total_state; i++)
	{
		A[i] = (float)(i+1)/total_state;
	}
	float k = (float)rand()/(float)RAND_MAX;
	while(k > A[j])
	{
		j++;
	}
	return j;
}
int main()
{
	long int i,j,k;
	long long int total=0;
	long long int org_total = 0;
	int current_state, next_state;
	FILE *A, *B;
	for(i=0; i<total_state; i++)
	{
		state[i] = first_state + i;
	}
	// Markov Chain Coading 
	printf("\n|-----Starting TPM Generation -----|\n");
	if(uniform_transition_mat == 1)
	{
		printf("|--TPM is made uniformally and is of order: %d --|\n", total_state);
	}
	else
	{
		printf("|--TPM is made non-uniformally and is of order: %d --|\n", total_state);
	}
	tranMat = maketranMat(total_state);
	PretranMat = makepretranMat(tranMat, total_state);
	printf("\n|--TPM generation complete TPM generated:\n");
	for(i=0; i<total_state; i++)
	{
		printf("\t");
		for(j=0;j<total_state; j++)
		{
			printf("%.5f ", tranMat[i][j]);
		}
		printf("\n");
	}
	A = fopen(file_name_for_seq_of_state,"w");

	printf("\n|-----Starting Sequence Generation of length: %d -----|\n", num_symbol);
	current_state = first();
	fprintf(A,"%c",state[current_state]);
	for(i=1; i<num_symbol; i++)
	{
		next_state = next(PretranMat, current_state);
		fprintf(A,"%c",state[next_state]);
		current_state = next_state;
	}
	fclose(A);
	printf("\n|-----Sequence Generation complete check file name ");
	printf(file_name_for_seq_of_state);
	printf(" for generated Sequence-----|\n");

	// Markov chain Coading complete 
	//Huffman Coading based on TPM
	printf("\n|-----Starting Huffmann coading with number of states: %d -----|\n", total_state);
	struct code ***codes = (struct code***)malloc(total_state*sizeof(struct code**));
	printf("\n|-----Preparing Huffman codes for each state based on TPM-----|");
	for(i=0; i<total_state; i++)
	{
		struct heapnode *root = makeHuffTree(tranMat[i], total_state);
		codes[i] = assigncode(root, total_state);
	}
	printf("\n|--Huffman code assignment complete---|");
	printf("\n|--Writing codes in file ");
	printf(file_name_for_huffcode);
	printf(" --|\n");
	A = fopen(file_name_for_seq_of_state,"r");
	B = fopen(file_name_for_huffcode, "w");
	char c, *arr;
	int r = log2(total_state);
	fscanf(A, "%c", &c);
	convert_bin(B, c-first_state, r);
	total = total + r;
	current_state = c-first_state;
	for(i=1; i<num_symbol; i++)
	{
		fscanf(A, "%c", &c);
		next_state = c-first_state;
		j = codes[current_state][next_state]->length;
		arr = codes[current_state][next_state]->seq;
		for(k=0; k<j; k++)
		{
			fprintf(B, "%c", arr[k]);
		}
		current_state = next_state;
		total = total +j;
	}
	printf("|--Code Writing complete--|\n");
	org_total = num_symbol*r;
	printf("\n|-----Printing bit per symbol-----|\n");
	printf("|- Bits per symbol: %.2f\n", (float)total/num_symbol);
	printf("\n|----- Huffman coading complete check file named ");
	printf(file_name_for_huffcode);
	printf(" for code-----|\n");
	fclose(A);
	fclose(B);
	free(tranMat);
	free(PretranMat);
	free(codes);
	return 0;
}
