#include "header.h"

int match(char *seq, int index, FILE *A, long long int *total)
{
	int i=0;
	int n=0, l=0, maxn=1, maxl=0, temp;
	while(i!= W)
	{
		if(seq[i] == seq[index])
		{
			temp = i;
			l = i;
			n = 0;
			while(seq[l] == seq[index+n])
			{
				l++;
				n++;
				if(index + n == num_symbol)
					break;
			}
			if(maxn <= n)
			{
				maxn = n;
				maxl = temp;
			}
		}
		i++;
	}
	if(maxn!=1)
	{
		l = W - maxl;
		int rl = log2(l);
		int rn = log2(maxn);
		rl++;
		rn++;
		for(i=0; i<rl;i++)
		{
			fprintf(A, "%c", '0');
		}
		convert_bin(A, l, rl);
		for(i=0; i<rn; i++)
		{
			fprintf(A, "%c", '0');
		}
		convert_bin(A, maxn, rn);
		*total = *total + 2*rn + 2*rl;
	}
	return maxn;
}

int main()
{

	FILE *A, *B;
	int i;
	//Leb-Z encoding
	B = fopen(file_name_for_LebZ_coading,"w");
	A = fopen(file_name_for_seq_of_state,"r");
	long long int total=0;
	long long int org_total = 0;
	//assign normal coading for symbols in Window
	int r = log2(total_state); //bits for representing each symbol
	org_total = (r+1)*num_symbol;
	char c;
	printf("\n|-----Running LZ encoding on seq of data stored in ");
	printf(file_name_for_seq_of_state);
	printf(" -----|\n");
	printf("|- Window Size: %d\n", W);
	printf("|- Total number of symbol in seq: %d\n", num_symbol);
	printf("\n|-----Encoding Window symbols by trival coading-----|\n");
	for(i=0; i<W; i++)
	{
		fscanf(A, "%c", &c);
		trival_coading(B, c-first_state, r);
		total = total + r + 1;
	}
	fclose(A);
	//assign code for l and n in lebZ coading
	printf("\n|-----Encoding rest symbols by LZ coading-----|\n");
	A = fopen(file_name_for_seq_of_state,"r");
	char *seq = (char *)malloc((num_symbol + 1)*sizeof(char));
	fgets(seq,num_symbol + 1,A);
	long int start = W;
	while(start < num_symbol)
	{
		int k;
		k = match(seq, start, B, &total);
		if(k != 1)
		{
			start = start + k;
		}
		else
		{
			c = seq[start];
			trival_coading(B, c-first_state, r);
			total = total + r + 1;
			start++;
		}
	}
	free(seq);
	fclose(A);
	fclose(B);
	printf("\n|-----Printing bit per symbol-----|\n");
	printf("|- Bits per symbol: %.2f\n", (float)total/num_symbol);
	printf("\n|----- LZ coading complete check file named ");
	printf(file_name_for_LebZ_coading);
	printf(" for code-----|\n");
	return 0;
}