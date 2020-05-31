#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define num_symbol 1000000 // i.e L
#define total_state 8 //i.e M
#define first_state 'A'
#define file_name_for_seq_of_state "state.txt"
#define file_name_for_LebZ_coading "LZencode.txt"
#define file_name_for_huffcode "huffmanencode.txt"
#define uniform_transition_mat 0 // make it 1 if you want uniform TPM else it will be random
#define W 10000

void convert_bin(FILE *A, int n, int c)
{
	int i,j;
	for (i = c-1; i >= 0; i--)
	{
		j = n >> i;

		if (j & 1)
		  fprintf(A,"%c", '1');
		else
		  fprintf(A,"%c", '0');
	}
}
void trival_coading(FILE *B, int n, int c)
{
	fprintf(B,"%c", '1');
	convert_bin(B, n, c);
}