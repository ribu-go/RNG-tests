#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <time.h>
#include <math.h>

#define L 255
												//Quantiles for alpha equal to
float alph0 = 2.327,		//0.01
		  alph1 = 1.645,		//0.05
			alph2 = 1.282;		//0.1

void zeroize (unsigned char *arr, int len) {
	
	for (int i = 0; i < len; i++)
		arr[i] = 0;
}

void zeroize_int (int *arr, int len) {
	
	for (int i = 0; i < len; i++)
		arr[i] = 0;
}

void shift_arr (char *arr, int len) {

	for (int i = 0; i < len - 1; i++)
		arr[i] = arr[i + 1];

	
}
/*for all generators defined below 'len' stands for number
of bytes in generated output sequence*/

/*built-in C random number generator*/
unsigned char *built_in_gen (int len) {

	time_t t;
	unsigned char *seq = (unsigned char *) malloc(len * sizeof(char));
	
	zeroize(seq, len);
	srand((unsigned) time(&t));

	for (int i = 0; i < len; i++) {
		seq[i] = rand() % 256;
	}
	
	return seq;
}

/*Lehmer_High and Lehmer_Low generator; mode is defined through input*/
unsigned char *lehmer (unsigned int len, uint32_t initval, char mode) {

	unsigned char *seq = (unsigned char *) malloc(len * sizeof(char));
	uint64_t m = 1UL << 32,
			 a = (1 << 16) + 1,
			 c = 119,
			 x_curr, x_prev = initval,
			 mask_l = (1 << 9) - 1,
			 mask_h = mask_l << 24;
			 
	zeroize(seq, len);

	for (int i = 0; i < len; i++) {

			x_curr = (a * x_prev + c) % m;
				
			//lehmer_low is default mode
			if (mode == 'h')
				seq[i] = (x_curr & mask_h) >> 24;
			else
				seq[i] = (x_curr & mask_l);
				
			x_prev = x_curr;
	}
	
	return seq;
}

/*L20 LSR generator*/
unsigned char *L20 (unsigned int len, char initval) {
	
	unsigned char *seq = (unsigned char *) malloc(len * sizeof(char));
	int state = initval, aux;
	
	zeroize(seq, len);

	for (int i = 0; i < len; i++) {
		
		for (int j = 0; j < 8; j++) {
			seq[i] |= (state & 1) << j;
			aux = (state ^ (state >> 11) ^ (state >> 15) ^ (state >> 17)) & 1;
			state = (state >> 1) | (aux << (20 - 1));
		}
	}
	
	return seq;
}

/*L89 LSR generator*/
unsigned char *L89 (unsigned int len, int initval) {

	unsigned char *seq = (unsigned char *) malloc(len * sizeof(char));
	uint64_t state[2], aux;
	
	state[0] = initval;
	state[1] = 0;
	zeroize(seq, len);

	for (int i = 0; i < len; i++) {

		for (int j = 0; j < 8; j++) {
			seq[i] |= (state[0] & 1) << j;
			aux = (state[0] ^ (state[0] >> 51)) & 1;
			
			state[0] = (state[0] >> 1) | (state[1] << 63);
			state[1] = (state[1] >> 1) | (aux << (25 - 1));
		}
	}

	return seq;
}

/*Geffe generator*/
unsigned char *Geffe (int len, int initval) {
	
	unsigned char *seq = (unsigned char *) malloc(len * sizeof(char));
	int state11 = initval,
		state9 = initval >> 11,
		state10 = initval >> 22,
		aux;

	zeroize(seq, len);
		
	for (int i = 0; i < len; i++) {
		
		for (int j = 0; j < 32; j++) {
			seq[i] |= (((state11 & state10) ^ ((1 ^ state10) & state9)) & 1) << j;
			
			aux = (state11 ^ (state11 >> 2)) & 1;
			state11 = (state11 >> 1) | (aux << (11 - 1));
			
			aux = (state9 ^ (state9 >> 1) ^ (state9 >> 3) ^ (state9 >> 4)) & 1;
			state9 = (state9 >> 1) | (aux << (9 - 1));
			
			aux = (state10 ^ (state10 >> 3)) & 1;
			state10 = (state10 >> 1) | (aux << (10 - 1));
		}
	}
	
	return seq;
}

/*Wolfram generator*/
unsigned char *Wolfram (int len, int initval) {
	
	unsigned char *seq = (unsigned char *) malloc(len * sizeof(char));
	int r_curr, r_prev = initval, aux1, aux2;
	
	zeroize(seq, len);

	for (int i = 0; i < len; i++) {
		
		for (int j = 0; j < 8; j++) {
			
			aux1 = (r_prev << 1) | (r_prev >> 31);
			aux2 = (r_prev >> 1) | (r_prev << 31);
			r_curr = aux1 ^ (r_prev | aux2);
			seq[i] |= (r_prev & 1) << j;
			r_prev = r_curr;
		}
	}
	
	return seq;
}

/*Librarian random number generator*/
unsigned char *librarian (int len, char *path) {
	
	unsigned char *seq = (unsigned char *) malloc(len * sizeof(char));
	FILE *file;
	file = fopen("randoutput.txt", "r");

	if (file != NULL) {
		fread(seq, 1, len, file);
	}
	else {
		puts("Can't read file.");
	}
	
	return seq;
}

void uniformity_test (unsigned char *seq, int len, float quantile, float alpha) {

	int bytes[256];
	int n_j = len / 256;
	float stat = 0,
			  hi_def = sqrt(2 * L) * quantile + L;
	unsigned num = 0;

	zeroize_int(bytes, 256);

	for (int i = 0; i < len; i++) {
		bytes[seq[i]]++;
	}

	for (int i = 0; i < 256; i++) {
		stat += (n_j - bytes[i]) * (n_j - bytes[i]) / n_j;
	}

	printf("\t%.4f\t%.4f\t%.4f\t", alpha, hi_def, stat);
	
	if (stat <= hi_def) {
		puts("true");
	}
	else 
		puts("false");

	return;
}

void independency_test (unsigned char *seq, int len, float quantile, float alpha) {

	int bytes_i[256],
			bytes_j[256],
			pairs[256][256];			 
	float stat = 0,
			  hi_def = sqrt(2 * L * L) * quantile + L * L;
	
	for (int i = 0; i < 256; i++) {
		bytes_i[i] = 1;
		bytes_j[i] = 1;
		zeroize_int(pairs[i], 256);
	}
	for (int i = 0; i < len / 2; i++) {
		pairs[seq[2 * i]][seq[2 * i + 1]]++;
		bytes_i[seq[2 * i]]++;
		bytes_j[seq[2 * i + 1]]++;
	}
	
	for (int i = 0; i < 256; i++) {
		for (int j = 0; j < 256; j++) {
			stat += pow(pairs[i][j], 2) / (bytes_i[i] * bytes_j[j]);
		}
	}
	
	stat = (stat - 1) * len / 2;

	printf("\t%.4f\t%.4f\t%.4f\t", alpha, hi_def, stat);
	
	if (stat <= hi_def) {
		puts("true");
	}
	else 
		puts("false");

	return;
}

void homogen_test (unsigned char *seq, int len, float quantile, float alpha) {
	
	int segm = 1000;
	int bytes[256], m = (int) len / segm;
	float stat = 0,
			  hi_def = sqrt(2 * L * L) * quantile + L * (segm - 1);
	
	zeroize_int(bytes, 256);
	for (int i = 0; i < len; i++) {
		bytes[seq[i]]++;
	}

	int *freq_by_segm[256];
	for (int i = 0; i < 256; i++) {
		freq_by_segm[i] = (int *) malloc(segm * sizeof(int));
		zeroize_int(freq_by_segm[i], segm);
	}

	for (int i = 0; i < len; i += m) {
		for (int j = 0; j < m; j++) {
			freq_by_segm[seq[i + j]][i / m]++;
		}
	}

	for (int i = 0; i < 256; i++) {
		for (int j = 0; j < segm; j++) {
			if (bytes[i] == 0)
				continue;
			else {
				stat += pow(freq_by_segm[i][j], 2) / (m * bytes[i]);
			}
		}
	}

	stat = (stat - 1) * segm * m;
	printf("\t%.4f\t%.4f\t%.4f\t", alpha, hi_def, stat);
	
	if (stat <= hi_def) {
		puts("true");
	}
	else 
		puts("false");

	return;
}

void test_sequence (unsigned char *seq, char *seq_name, int len) {

	printf("\nGenerator: %s\n", seq_name);
	
	puts("Uniformity:");
	puts("\talpha\tCB\t\t\tAS");
	uniformity_test(seq, len, alph2, 0.1);
	uniformity_test(seq, len, alph1, 0.05);
	uniformity_test(seq, len, alph0, 0.01);
	
	puts("Independency:");
	puts("\talpha\tCB\t\t\tAS");
	independency_test(seq, len, alph2, 0.1);
	independency_test(seq, len, alph1, 0.05);
	independency_test(seq, len, alph0, 0.01);
	
	puts("Homogenity:");
	puts("\talpha\tCB\t\t\tAS");
	homogen_test(seq, len, alph2, 0.1);
	homogen_test(seq, len, alph1, 0.05);
	homogen_test(seq, len, alph0, 0.01);

	return;
}

void seq_sample (unsigned char *seq, char *string, int len) {

	printf("\n%s : ", string);

	for (int i = 0; i < len; i++)
		for (int j = 0; j < 8; j++)
			printf("%d", (unsigned) (seq[i] >> j) % 2);
	
	puts("");
}

int main () {
	
	//defining ptrs to store generated sequences
	//names are self-explainatory
	unsigned char *seq_built_in,
			 *seq_lehmer_l,
			 *seq_lehmer_h,
			 *seq_L20,
			 *seq_L89,
			 *seq_Geffe,
			 *seq_librarian,
			 *seq_Wolfram,
			 *seq_BM,
			 *seq_BM_bytes,
			 *seq_BBS,
			 *seq_BBS_bytes;
	int len = 1000000,  //one million is concidered enough for our
											//little research	
			len_sample = 20;

	//generating initial state that will be used for all generators
	time_t t;
	srand((unsigned) time(&t));
	int initial_state = rand();

	printf("Seed: %d\n", initial_state);
	//generating sequences
	seq_built_in = built_in_gen(len);
	seq_sample(seq_built_in, "Built-in random", len_sample);
	
	seq_lehmer_h = lehmer(len, initial_state, 'h');
	seq_sample(seq_lehmer_h, "Lehmer_High", len_sample);
	
	seq_lehmer_l = lehmer(len, initial_state, 'l');
	seq_sample(seq_lehmer_l, "Lehmer_Low", len_sample);
	
	seq_L20 = L20(len, initial_state);
	seq_sample(seq_L20, "L20", len_sample);
	
	seq_L89 = L89(len, initial_state);
	seq_sample(seq_L89, "L89", len_sample);
	
	seq_librarian = librarian(200000, "text.txt");
	seq_sample(seq_librarian, "Librarian", len_sample);
	
	seq_Geffe = Geffe(len, initial_state);
	seq_sample(seq_Geffe, "Geffe", len_sample);
	
	seq_Wolfram = Wolfram(len, initial_state);
	seq_sample(seq_Wolfram, "Wolfram", len_sample);
	
	//computing statistical characteristics
	test_sequence(seq_built_in, "built-in", len);
	test_sequence(seq_lehmer_h, "Lehmer high", len);
	test_sequence(seq_lehmer_l, "Lehmer low", len);
	test_sequence(seq_L20, "L20", len);
	test_sequence(seq_L89, "L89", len);
	test_sequence(seq_librarian, "librarian", 200000);
	test_sequence(seq_Geffe, "Geffe", len);
	test_sequence(seq_Wolfram, "Wolfram", len);
	
	return 0;
}
